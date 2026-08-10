#!/usr/bin/env python3
"""Build app-ready, gene-partitioned Parquet files from mouse VIS pseudobulk H5AD.

Input
-----
A pseudobulk AnnData object with:
  * rows = library_label x cluster_alias pseudobulk observations
  * columns = genes
  * X = summed raw counts
  * obs = library/donor/taxonomy metadata, including library_label

Output
------
<output-dir>/
  counts_by_gene/
    gene=<gene_key>/
      data.parquet
  metadata/
    obs_metadata.parquet
    gene_map.parquet
    CPM_scaling.csv
  statistics/
    mouse_gene_summary_statistics.csv
  manifest.json

Each gene has exactly one Parquet file. All gene files use the same sample order.
CPM_scaling.csv has exactly two columns, as requested: library_label and
CPM_scaling_factor. Because each H5AD row is a library x cluster pseudobulk,
library_label may repeat; row order is identical to obs_metadata.parquet and
every gene Parquet file.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
from pathlib import Path
from typing import Iterable, Optional

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse


DEFAULT_INPUT = Path(
    "/data/mouse_vis_cortex_pseudobulk/"
    "Developing_Mouse_Visual_Cortex_library_cluster_pseudobulk.h5ad"
)
DEFAULT_OUTPUT = Path("/results/mouse_vis_cortex_parquet")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the Developing Mouse VIS pseudobulk H5AD into one "
            "Parquet file per gene plus metadata, CPM scaling, and compact "
            "gene statistics."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input pseudobulk H5AD (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output directory (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--gene-block-size",
        type=int,
        default=256,
        help="Number of genes processed together for statistics (default: 256).",
    )
    parser.add_argument(
        "--compression",
        choices=("zstd", "snappy", "gzip", "none"),
        default="zstd",
        help="Parquet compression codec (default: zstd).",
    )
    parser.add_argument(
        "--min-nonzero-n",
        type=int,
        default=3,
        help=(
            "Minimum positive pseudobulk observations required in a group "
            "for nonzero-only statistics (default: 3)."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete an existing output directory before writing.",
    )
    return parser.parse_args()


def locate_input(requested: Path) -> Path:
    if requested.exists():
        return requested

    candidates = list(Path("/data").glob("**/*pseudobulk*.h5ad"))
    if len(candidates) == 1:
        print(f"Input not found at {requested}; using {candidates[0]}")
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(
            f"Input H5AD not found: {requested}. No pseudobulk H5AD was found under /data."
        )
    raise FileNotFoundError(
        f"Input H5AD not found: {requested}. Multiple candidates were found under /data; "
        "rerun with --input followed by the intended path:\n  "
        + "\n  ".join(str(x) for x in candidates)
    )


def ensure_output(output_dir: Path, overwrite: bool) -> tuple[Path, Path, Path]:
    if output_dir.exists():
        if not overwrite:
            raise FileExistsError(
                f"Output directory already exists: {output_dir}. "
                "Use --overwrite to replace it."
            )
        shutil.rmtree(output_dir)

    counts_dir = output_dir / "counts_by_gene"
    metadata_dir = output_dir / "metadata"
    statistics_dir = output_dir / "statistics"
    counts_dir.mkdir(parents=True)
    metadata_dir.mkdir(parents=True)
    statistics_dir.mkdir(parents=True)
    return counts_dir, metadata_dir, statistics_dir


def make_unique_keys(values: Iterable[str]) -> list[str]:
    """Create filesystem-safe, deterministic, unique gene partition keys."""
    seen: dict[str, int] = {}
    keys: list[str] = []
    for raw in values:
        value = str(raw).strip() or "unnamed_gene"
        base = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
        base = base or "unnamed_gene"
        count = seen.get(base, 0) + 1
        seen[base] = count
        keys.append(base if count == 1 else f"{base}__{count}")
    return keys


def choose_column(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    column_set = set(columns)
    for candidate in candidates:
        if candidate in column_set:
            return candidate
    return None


def infer_gene_symbols(adata: ad.AnnData) -> pd.Series:
    for field in ("gene_symbol", "gene_symbols", "symbol", "gene_name"):
        if field in adata.var.columns:
            values = adata.var[field].astype("string")
            fallback = pd.Series(adata.var_names.astype(str), index=adata.var_names)
            return values.fillna(fallback).replace("", pd.NA).fillna(fallback)
    return pd.Series(adata.var_names.astype(str), index=adata.var_names, dtype="string")


def parse_age_value(value: object) -> float:
    """Return an ordering value for numeric, embryonic (E), or postnatal (P) ages.

    Embryonic ages are mapped below zero and postnatal ages at/above zero. The
    exact embryonic offset does not affect youngest/oldest group identification.
    """
    if value is None or pd.isna(value):
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)

    text = str(value).strip().upper()
    if not text:
        return np.nan

    match = re.search(r"([EP])\s*([0-9]+(?:\.[0-9]+)?)", text)
    if match:
        stage, number = match.groups()
        number = float(number)
        return number if stage == "P" else number - 100.0

    match = re.search(r"-?[0-9]+(?:\.[0-9]+)?", text)
    return float(match.group()) if match else np.nan


def safe_log2fc(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.log2((numerator + 1.0) / (denominator + 1.0))


def group_means(
    values: np.ndarray,
    groups: np.ndarray,
    minimum_nonzero_n: int,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Return all-observation and positive-only means by group.

    values has shape observations x genes.
    """
    valid = pd.notna(groups)
    labels = pd.unique(groups[valid].astype(str)).tolist()
    all_means: list[np.ndarray] = []
    nz_means: list[np.ndarray] = []

    for label in labels:
        mask = valid & (groups.astype(str) == label)
        block = values[mask, :]
        all_means.append(np.mean(block, axis=0))

        positive = block > 0
        counts = positive.sum(axis=0)
        sums = np.where(positive, block, 0.0).sum(axis=0)
        out = np.full(block.shape[1], np.nan, dtype=np.float64)
        eligible = counts >= minimum_nonzero_n
        out[eligible] = sums[eligible] / counts[eligible]
        nz_means.append(out)

    return labels, np.vstack(all_means), np.vstack(nz_means)


def max_group_stats(
    labels: list[str],
    means: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return max group label and log2(max / arithmetic mean across groups)."""
    n_genes = means.shape[1]
    names = np.empty(n_genes, dtype=object)
    log2fc = np.full(n_genes, np.nan, dtype=np.float64)

    for j in range(n_genes):
        column = means[:, j]
        valid = np.isfinite(column)
        if not np.any(valid):
            names[j] = None
            continue
        valid_indices = np.flatnonzero(valid)
        max_index = valid_indices[np.argmax(column[valid])]
        max_value = column[max_index]
        mean_value = np.mean(column[valid])
        names[j] = labels[max_index]
        log2fc[j] = safe_log2fc(
            np.asarray([max_value]), np.asarray([mean_value])
        )[0]
    return names, log2fc


def write_gene_file(
    counts_dir: Path,
    gene_key: str,
    sample_ids: np.ndarray,
    counts: np.ndarray,
    compression: Optional[str],
) -> None:
    partition_dir = counts_dir / f"gene={gene_key}"
    partition_dir.mkdir(parents=True, exist_ok=False)
    table = pa.table(
        {
            "sample_id": pa.array(sample_ids, type=pa.string()),
            "summed_counts": pa.array(counts),
        }
    )
    pq.write_table(
        table,
        partition_dir / "data.parquet",
        compression=compression,
        use_dictionary=["sample_id"],
        write_statistics=True,
    )


def main() -> None:
    args = parse_args()
    input_path = locate_input(args.input)
    output_dir = args.output_dir
    compression = None if args.compression == "none" else args.compression

    counts_dir, metadata_dir, statistics_dir = ensure_output(
        output_dir, args.overwrite
    )

    print(f"Reading pseudobulk H5AD: {input_path}")
    adata = ad.read_h5ad(input_path)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"The input H5AD is empty: shape={adata.shape}")
    if "library_label" not in adata.obs.columns:
        raise KeyError("Input H5AD obs must contain 'library_label'.")

    # Retain a stable pseudobulk row identifier and exact row order everywhere.
    sample_ids = adata.obs_names.astype(str).to_numpy()
    if len(np.unique(sample_ids)) != len(sample_ids):
        raise ValueError("AnnData obs_names must be unique.")

    obs_metadata = adata.obs.copy()
    obs_metadata.insert(0, "sample_id", sample_ids)
    obs_metadata.reset_index(drop=True, inplace=True)

    gene_symbols = infer_gene_symbols(adata).astype(str).to_numpy()
    gene_keys = make_unique_keys(adata.var_names.astype(str))
    if len(set(gene_keys)) != len(gene_keys):
        raise RuntimeError("Generated gene keys are unexpectedly non-unique.")

    gene_map = pd.DataFrame(
        {
            "gene_symbol": gene_symbols,
            "gene_column": gene_keys,
            "gene_key": gene_keys,
            "gene_index": np.arange(adata.n_vars, dtype=np.int32),
            "source_var_name": adata.var_names.astype(str),
        }
    )

    # Convert once to CSC for efficient extraction of one complete gene column.
    if sparse.issparse(adata.X):
        X = adata.X.tocsc(copy=False)
    else:
        X = np.asarray(adata.X)

    # Library-size denominator is total raw counts across all genes for each
    # pseudobulk observation (library_label x cluster_alias row).
    row_totals = np.asarray(X.sum(axis=1)).reshape(-1).astype(np.float64)
    cpm_scaling = np.zeros_like(row_totals, dtype=np.float64)
    positive_library_size = row_totals > 0
    cpm_scaling[positive_library_size] = 1_000_000.0 / row_totals[positive_library_size]

    cpm_table = pd.DataFrame(
        {
            "library_label": obs_metadata["library_label"].astype(str).to_numpy(),
            "CPM_scaling_factor": cpm_scaling,
        }
    )

    print("Writing metadata files...")
    obs_metadata.to_parquet(
        metadata_dir / "obs_metadata.parquet",
        index=False,
        compression=compression,
    )
    gene_map.to_parquet(
        metadata_dir / "gene_map.parquet",
        index=False,
        compression=compression,
    )
    cpm_table.to_csv(metadata_dir / "CPM_scaling.csv", index=False)

    age_field = choose_column(
        obs_metadata.columns,
        (
            "age",
            "age_label",
            "donor_age",
            "developmental_age",
            "age_in_days",
            "age_days",
        ),
    )
    region_field = choose_column(
        obs_metadata.columns,
        (
            "region_of_interest_label",
            "Brain Region",
            "brain_region",
            "region_label",
            "structure",
        ),
    )
    cluster_field = choose_column(
        obs_metadata.columns,
        ("cluster_alias", "cluster", "cluster_label"),
    )

    if age_field is None:
        raise KeyError(
            "Could not identify an age field. Available obs columns:\n  "
            + "\n  ".join(obs_metadata.columns)
        )
    if region_field is None:
        raise KeyError(
            "Could not identify a brain-region field. Available obs columns:\n  "
            + "\n  ".join(obs_metadata.columns)
        )
    if cluster_field is None:
        raise KeyError(
            "Could not identify a cluster field. Available obs columns:\n  "
            + "\n  ".join(obs_metadata.columns)
        )

    age_labels = obs_metadata[age_field].astype("string").to_numpy()
    age_values = np.asarray([parse_age_value(x) for x in age_labels], dtype=np.float64)
    finite_age = np.isfinite(age_values)
    if np.unique(age_values[finite_age]).size < 2:
        raise ValueError(
            f"Age field '{age_field}' does not contain at least two parseable ages."
        )
    youngest_value = np.min(age_values[finite_age])
    oldest_value = np.max(age_values[finite_age])
    youngest_mask = finite_age & (age_values == youngest_value)
    oldest_mask = finite_age & (age_values == oldest_value)
    youngest_label = str(pd.unique(age_labels[youngest_mask])[0])
    oldest_label = str(pd.unique(age_labels[oldest_mask])[0])

    region_groups = obs_metadata[region_field].astype("string").to_numpy()
    cluster_groups = obs_metadata[cluster_field].astype("string").to_numpy()

    print(
        f"Detected age='{age_field}' ({youngest_label} to {oldest_label}), "
        f"region='{region_field}', cluster='{cluster_field}'."
    )
    print(
        f"Writing {adata.n_vars:,} single-file gene partitions and statistics "
        f"for {adata.n_obs:,} pseudobulk observations..."
    )

    statistics_blocks: list[pd.DataFrame] = []
    block_size = max(1, int(args.gene_block_size))

    for start in range(0, adata.n_vars, block_size):
        stop = min(start + block_size, adata.n_vars)
        if sparse.issparse(X):
            raw = X[:, start:stop].toarray()
        else:
            raw = np.asarray(X[:, start:stop])
        raw = raw.astype(np.float64, copy=False)
        cpm = raw * cpm_scaling[:, None]

        # Each gene is written exactly once to one data.parquet file.
        for local_j, gene_j in enumerate(range(start, stop)):
            gene_counts = raw[:, local_j]
            if np.all(np.equal(gene_counts, np.floor(gene_counts))):
                gene_counts = gene_counts.astype(np.int64)
            write_gene_file(
                counts_dir=counts_dir,
                gene_key=gene_keys[gene_j],
                sample_ids=sample_ids,
                counts=gene_counts,
                compression=compression,
            )

        overall_mean = np.mean(cpm, axis=0)
        positive = raw > 0
        positive_n = positive.sum(axis=0)
        overall_nonzero_mean = np.full(stop - start, np.nan, dtype=np.float64)
        eligible_overall = positive_n >= args.min_nonzero_n
        overall_nonzero_mean[eligible_overall] = (
            np.where(positive, cpm, 0.0).sum(axis=0)[eligible_overall]
            / positive_n[eligible_overall]
        )

        youngest_mean = np.mean(cpm[youngest_mask, :], axis=0)
        oldest_mean = np.mean(cpm[oldest_mask, :], axis=0)
        age_log2fc = safe_log2fc(oldest_mean, youngest_mean)

        youngest_positive = raw[youngest_mask, :] > 0
        oldest_positive = raw[oldest_mask, :] > 0
        youngest_n = youngest_positive.sum(axis=0)
        oldest_n = oldest_positive.sum(axis=0)
        age_log2fc_nonzero = np.full(stop - start, np.nan, dtype=np.float64)
        eligible_age = (
            (youngest_n >= args.min_nonzero_n)
            & (oldest_n >= args.min_nonzero_n)
        )
        youngest_nz_mean = np.full(stop - start, np.nan, dtype=np.float64)
        oldest_nz_mean = np.full(stop - start, np.nan, dtype=np.float64)
        youngest_nz_mean[eligible_age] = (
            np.where(youngest_positive, cpm[youngest_mask, :], 0.0)
            .sum(axis=0)[eligible_age]
            / youngest_n[eligible_age]
        )
        oldest_nz_mean[eligible_age] = (
            np.where(oldest_positive, cpm[oldest_mask, :], 0.0)
            .sum(axis=0)[eligible_age]
            / oldest_n[eligible_age]
        )
        age_log2fc_nonzero[eligible_age] = safe_log2fc(
            oldest_nz_mean[eligible_age], youngest_nz_mean[eligible_age]
        )

        region_labels, region_means, region_nz_means = group_means(
            cpm, region_groups, args.min_nonzero_n
        )
        max_region, region_log2fc = max_group_stats(region_labels, region_means)
        max_region_nz, region_log2fc_nz = max_group_stats(
            region_labels, region_nz_means
        )

        cluster_labels, cluster_means, cluster_nz_means = group_means(
            cpm, cluster_groups, args.min_nonzero_n
        )
        max_cluster, cluster_log2fc = max_group_stats(
            cluster_labels, cluster_means
        )
        max_cluster_nz, cluster_log2fc_nz = max_group_stats(
            cluster_labels, cluster_nz_means
        )

        statistics_blocks.append(
            pd.DataFrame(
                {
                    "gene_symbol": gene_symbols[start:stop],
                    "gene_column": gene_keys[start:stop],
                    "overall_mean_CPM": overall_mean,
                    "overall_mean_CPM_nonzero": overall_nonzero_mean,
                    "oldest_vs_youngest_log2fc": age_log2fc,
                    "oldest_vs_youngest_log2fc_nonzero": age_log2fc_nonzero,
                    "max_brain_region": max_region,
                    "max_brain_region_nonzero": max_region_nz,
                    "max_brain_region_vs_mean_region_log2fc": region_log2fc,
                    "max_brain_region_vs_mean_region_log2fc_nonzero": region_log2fc_nz,
                    "max_cluster": max_cluster,
                    "max_cluster_nonzero": max_cluster_nz,
                    "max_cluster_vs_mean_cluster_log2fc": cluster_log2fc,
                    "max_cluster_vs_mean_cluster_log2fc_nonzero": cluster_log2fc_nz,
                }
            )
        )

        print(f"  completed genes {start + 1:,}-{stop:,} of {adata.n_vars:,}")

    statistics = pd.concat(statistics_blocks, ignore_index=True)
    statistics.to_csv(
        statistics_dir / "mouse_gene_summary_statistics.csv",
        index=False,
    )

    manifest = {
        "input_h5ad": str(input_path),
        "n_observations": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "gene_files": int(adata.n_vars),
        "gene_file_format": "counts_by_gene/gene=<gene_key>/data.parquet",
        "sample_order": "Identical across obs_metadata.parquet, CPM_scaling.csv, and every gene file",
        "cpm_definition": "summed_counts * (1,000,000 / total raw counts across all genes in the pseudobulk row)",
        "age_field": age_field,
        "youngest_age": youngest_label,
        "oldest_age": oldest_label,
        "region_field": region_field,
        "cluster_field": cluster_field,
        "minimum_nonzero_group_n": int(args.min_nonzero_n),
        "compression": args.compression,
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    # Read-back checks for the files the app will depend on.
    check_obs = pd.read_parquet(metadata_dir / "obs_metadata.parquet")
    check_map = pd.read_parquet(metadata_dir / "gene_map.parquet")
    check_cpm = pd.read_csv(metadata_dir / "CPM_scaling.csv")
    if len(check_obs) != adata.n_obs or len(check_cpm) != adata.n_obs:
        raise RuntimeError("Metadata or CPM row count changed during serialization.")
    if len(check_map) != adata.n_vars:
        raise RuntimeError("Gene-map row count changed during serialization.")
    if not np.array_equal(
        check_obs["library_label"].astype(str).to_numpy(),
        check_cpm["library_label"].astype(str).to_numpy(),
    ):
        raise RuntimeError("CPM_scaling.csv library order does not match metadata.")

    print("\nCompleted successfully.")
    print(f"Output directory: {output_dir}")
    print(f"Gene files: {adata.n_vars:,}")
    print(f"Pseudobulk rows per gene: {adata.n_obs:,}")
    print(f"Statistics: {statistics_dir / 'mouse_gene_summary_statistics.csv'}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
