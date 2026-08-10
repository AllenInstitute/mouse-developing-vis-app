#!/usr/bin/env python3
"""Create a library-by-cluster pseudobulk H5AD for Developing Mouse VIS cortex.

Output
------
./Developing_Mouse_Visual_Cortex_library_cluster_pseudobulk.h5ad

Each observation is one library_label x cluster_alias combination. X contains
summed raw counts by gene. obs contains library, donor, cell, and taxonomy
metadata taken from the first cell in each pseudobulk group, plus
number_of_cells, the number of cells contributing to every value in that row.

The script downloads data through AbcProjectCache. The large source H5AD is
opened in backed mode and aggregated in row chunks rather than loaded into RAM.

Suggested environment
---------------------
pip install "abc_atlas_access[notebooks] @ git+https://github.com/AllenInstitute/abc_atlas_access.git" \
            anndata pandas numpy scipy h5py
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache


DATASET_DIR = "Developing-Mouse-Vis-Cortex-10X"
TAXONOMY_DIR = "Developing-Mouse-Vis-Cortex-taxonomy"
DEFAULT_MANIFEST = "releases/20260331/manifest.json"
DEFAULT_OUTPUT = "Developing_Mouse_Visual_Cortex_library_cluster_pseudobulk.h5ad"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download the Developing Mouse Visual Cortex raw-count H5AD and "
            "create library-by-cluster pseudobulk summed counts."
        )
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("abc_atlas_cache"),
        help="ABC Atlas download cache (default: ./abc_atlas_cache).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT),
        help=f"Output H5AD (default: ./{DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--manifest",
        default=DEFAULT_MANIFEST,
        help=(
            "ABC Atlas manifest to request. If unavailable to the installed "
            "cache, use --manifest latest."
        ),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=4096,
        help="Number of cells read from the source H5AD per chunk.",
    )
    return parser.parse_args()


def load_manifest(cache: AbcProjectCache, requested: str) -> None:
    """Load a requested manifest, or retain the cache default for 'latest'."""
    if requested.lower() == "latest":
        logging.info("Using cache default manifest: %s", cache.current_manifest)
        return

    # Some abc_atlas_access releases require a non-default manifest name to be
    # registered before load_manifest() can see it.
    if requested not in cache.cache.manifest_file_names:
        cache.cache.manifest_file_names.append(requested)
    cache.load_manifest(requested)
    logging.info("Using manifest: %s", cache.current_manifest)


def get_raw_h5ad_path(cache: AbcProjectCache) -> Path:
    """Resolve and download the raw-count expression H5AD."""
    if hasattr(cache, "list_expression_matrix_files"):
        names: Iterable[str] = cache.list_expression_matrix_files(
            directory=DATASET_DIR
        )
    else:
        names = cache.list_data_files(directory=DATASET_DIR)

    names = list(names)
    raw_names = [
        name for name in names
        if str(name).rstrip("/").endswith("/raw") or str(name) == "raw"
    ]
    if len(raw_names) != 1:
        raise RuntimeError(
            f"Expected exactly one raw expression matrix in {DATASET_DIR}; "
            f"found {raw_names!r}. Available entries: {names!r}"
        )

    logging.info("Downloading/resolving raw expression matrix: %s", raw_names[0])
    return Path(
        cache.get_data_path(directory=DATASET_DIR, file_name=raw_names[0])
    )


def load_cell_and_library_metadata(cache: AbcProjectCache) -> pd.DataFrame:
    """Build cell metadata with library, donor, and taxonomy annotations."""
    logging.info("Loading cell, library, donor, and taxonomy metadata")

    cell = cache.get_metadata_dataframe(
        directory=DATASET_DIR, file_name="cell_metadata"
    ).set_index("cell_label", drop=True)
    cell.index = cell.index.astype(str)
    cell.index.name = "cell_label"

    library = cache.get_metadata_dataframe(
        directory=DATASET_DIR, file_name="library"
    ).set_index("library_label", drop=True)
    library.index = library.index.astype(str)

    donor = cache.get_metadata_dataframe(
        directory=DATASET_DIR, file_name="donor"
    ).set_index("donor_label", drop=True)
    donor.index = donor.index.astype(str)

    membership = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIR, file_name="cell_to_cluster_membership"
    ).set_index("cell_label", drop=True)
    membership.index = membership.index.astype(str)

    if "cluster_alias" not in membership.columns:
        raise KeyError("cell_to_cluster_membership lacks cluster_alias")

    # Reproduce the notebook's cluster annotation expansion. The resulting
    # columns normally include class, subclass, cluster, and subcluster terms.
    cluster_term = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIR, file_name="cluster_annotation_term"
    ).set_index("label", drop=True)

    cluster_membership = cache.get_metadata_dataframe(
        directory=TAXONOMY_DIR,
        file_name="cluster_to_cluster_annotation_membership",
    ).set_index("cluster_annotation_term_label", drop=True)

    membership_with_terms = cluster_membership.join(
        cluster_term,
        how="left",
        rsuffix="_annotation_term",
    ).reset_index()

    required_term_columns = {
        "cluster_alias",
        "cluster_annotation_term_set_name",
        "cluster_annotation_term_name",
    }
    missing = required_term_columns.difference(membership_with_terms.columns)
    if missing:
        raise KeyError(
            "Taxonomy annotation tables lack required columns: "
            + ", ".join(sorted(missing))
        )

    cluster_details = (
        membership_with_terms.groupby(
            ["cluster_alias", "cluster_annotation_term_set_name"],
            observed=True,
        )["cluster_annotation_term_name"]
        .first()
        .unstack()
    )

    # Join in the same order used by the official notebooks.
    extended = cell.join(membership[["cluster_alias"]], how="inner")
    extended = extended.join(cluster_details, on="cluster_alias")
    extended = extended.join(library, on="library_label", rsuffix="_library")

    donor_join_column = (
        "donor_label_library"
        if "donor_label_library" in extended.columns
        else "donor_label"
    )
    extended = extended.join(
        donor,
        on=donor_join_column,
        rsuffix="_donor",
    )

    if "library_label" not in extended.columns:
        raise KeyError("Expanded metadata lacks library_label")
    if extended["cluster_alias"].isna().any():
        raise ValueError("Some cells lack cluster_alias assignments")

    extended["library_label"] = extended["library_label"].astype(str)
    extended["cluster_alias"] = extended["cluster_alias"].astype(str)
    return extended


def gene_metadata_for_source(
    cache: AbcProjectCache,
    source: ad.AnnData,
) -> pd.DataFrame:
    """Create var metadata in source-column order, indexed by gene symbol."""
    gene = cache.get_metadata_dataframe(
        directory=DATASET_DIR, file_name="gene"
    )

    if "gene_identifier" not in gene.columns or "gene_symbol" not in gene.columns:
        raise KeyError("gene.csv must contain gene_identifier and gene_symbol")

    gene = gene.drop_duplicates("gene_identifier").set_index("gene_identifier")
    source_ids = pd.Index(source.var_names.astype(str), name="gene_identifier")

    missing = source_ids.difference(gene.index.astype(str))
    if len(missing):
        raise ValueError(
            f"{len(missing)} source gene identifiers are absent from gene.csv; "
            f"examples: {missing[:10].tolist()}"
        )

    gene.index = gene.index.astype(str)
    var = gene.reindex(source_ids).copy()
    symbols = var["gene_symbol"].astype("string")
    symbols = symbols.fillna(pd.Series(source_ids, index=source_ids)).astype(str)

    # AnnData requires unique var_names. Preserve the original symbol and use
    # a stable identifier suffix only when symbols are duplicated.
    symbol_counts = symbols.groupby(symbols).cumcount()
    duplicate_totals = symbols.map(symbols.value_counts())
    unique_symbols = symbols.where(
        duplicate_totals.eq(1),
        symbols + "__" + source_ids.astype(str),
    )

    var.insert(0, "gene_symbol_original", symbols.to_numpy())
    var.index = pd.Index(unique_symbols.to_numpy(), name="gene_key")
    return var


def make_group_metadata(
    aligned_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Define library x cluster pseudobulk groups in deterministic order."""
    group_frame = aligned_metadata[["library_label", "cluster_alias"]].copy()
    group_index = pd.MultiIndex.from_frame(group_frame)
    codes, unique_groups = pd.factorize(group_index, sort=True)

    group_keys = pd.DataFrame(
        unique_groups.tolist(),
        columns=["library_label", "cluster_alias"],
    )
    group_keys["_group_code"] = np.arange(len(group_keys), dtype=np.int64)

    first_metadata = (
        aligned_metadata.assign(_group_code=codes)
        .groupby("_group_code", sort=True, observed=True)
        .first()
    )
    first_metadata = first_metadata.reindex(group_keys["_group_code"])

    # Ensure grouping fields are the first two columns and not duplicated.
    first_metadata = first_metadata.drop(
        columns=["library_label", "cluster_alias"], errors="ignore"
    )
    obs = pd.concat(
        [group_keys.set_index("_group_code"), first_metadata],
        axis=1,
    )

    number_of_cells = np.bincount(codes, minlength=len(obs)).astype(np.int64)
    obs.insert(2, "number_of_cells", number_of_cells)
    obs.index = pd.Index(
        [
            f"{library_label}__cluster_{cluster_alias}"
            for library_label, cluster_alias in zip(
                obs["library_label"], obs["cluster_alias"]
            )
        ],
        name="pseudobulk_id",
    )
    return obs, codes.astype(np.int64, copy=False)


def aggregate_counts(
    source: ad.AnnData,
    group_codes: np.ndarray,
    n_groups: int,
    chunk_size: int,
) -> sparse.csr_matrix:
    """Sum source rows into pseudobulk groups without loading all X into RAM."""
    n_cells, n_genes = source.shape
    if len(group_codes) != n_cells:
        raise ValueError("group_codes length does not match source cell count")

    total = sparse.csr_matrix((n_groups, n_genes), dtype=np.int64)

    for start in range(0, n_cells, chunk_size):
        stop = min(start + chunk_size, n_cells)
        logging.info("Aggregating cells %s-%s of %s", start + 1, stop, n_cells)

        block = source.X[start:stop]
        if sparse.issparse(block):
            block = block.tocsr()
        else:
            block = sparse.csr_matrix(np.asarray(block))

        # Raw counts should be integer-valued. Validate before conversion.
        if block.nnz and not np.allclose(block.data, np.rint(block.data)):
            raise ValueError(
                "The selected source matrix is not integer-valued; ensure the "
                "raw matrix, not log2-normalized data, was selected."
            )
        block.data = np.rint(block.data).astype(np.int64, copy=False)

        local_codes = group_codes[start:stop]
        assignment = sparse.csr_matrix(
            (
                np.ones(stop - start, dtype=np.int64),
                (local_codes, np.arange(stop - start, dtype=np.int64)),
            ),
            shape=(n_groups, stop - start),
        )
        total = total + assignment @ block

    total.sum_duplicates()
    total.eliminate_zeros()
    return total


def sanitize_dataframe_for_h5ad(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert mixed object columns to nullable strings for reliable H5AD I/O."""
    result = frame.copy()
    for column in result.columns:
        if result[column].dtype == object:
            result[column] = result[column].astype("string")
    return result


def main() -> None:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be at least 1")

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    args.output = args.output.resolve()

    cache = AbcProjectCache.from_cache_dir(args.cache_dir)
    load_manifest(cache, args.manifest)

    cell_metadata = load_cell_and_library_metadata(cache)
    raw_h5ad_path = get_raw_h5ad_path(cache)

    logging.info("Opening raw H5AD in backed mode: %s", raw_h5ad_path)
    source = ad.read_h5ad(raw_h5ad_path, backed="r")

    try:
        source_cells = pd.Index(source.obs_names.astype(str), name="cell_label")
        missing_cells = source_cells.difference(cell_metadata.index)
        if len(missing_cells):
            raise ValueError(
                f"{len(missing_cells)} source cells lack expanded metadata; "
                f"examples: {missing_cells[:10].tolist()}"
            )

        # Reindex metadata exactly to source X row order.
        aligned_metadata = cell_metadata.reindex(source_cells)
        obs, group_codes = make_group_metadata(aligned_metadata)
        var = gene_metadata_for_source(cache, source)

        logging.info(
            "Creating %s pseudobulk rows from %s cells across %s genes",
            len(obs), source.n_obs, source.n_vars,
        )
        summed_counts = aggregate_counts(
            source=source,
            group_codes=group_codes,
            n_groups=len(obs),
            chunk_size=args.chunk_size,
        )
    finally:
        # Backed AnnData stores an open h5py handle.
        if getattr(source, "file", None) is not None:
            source.file.close()

    obs = sanitize_dataframe_for_h5ad(obs)
    var = sanitize_dataframe_for_h5ad(var)

    pseudobulk = ad.AnnData(
        X=summed_counts,
        obs=obs,
        var=var,
        uns={
            "dataset": DATASET_DIR,
            "taxonomy": TAXONOMY_DIR,
            "manifest": str(cache.current_manifest),
            "aggregation": "sum of raw counts by library_label and cluster_alias",
            "grouping_columns": ["library_label", "cluster_alias"],
            "count_column": "number_of_cells",
            "source_h5ad": str(raw_h5ad_path),
        },
    )

    logging.info("Writing pseudobulk H5AD: %s", args.output)
    pseudobulk.write_h5ad(args.output, compression="gzip")

    # Final read-back validation catches serialization or schema problems.
    check = ad.read_h5ad(args.output, backed="r")
    try:
        if check.shape != pseudobulk.shape:
            raise RuntimeError(
                f"Read-back shape {check.shape} differs from {pseudobulk.shape}"
            )
        if "number_of_cells" not in check.obs.columns:
            raise RuntimeError("Read-back output lacks number_of_cells")
    finally:
        check.file.close()

    logging.info(
        "Complete: %s pseudobulk rows x %s genes written to %s",
        pseudobulk.n_obs, pseudobulk.n_vars, args.output,
    )


if __name__ == "__main__":
    main()
