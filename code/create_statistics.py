#!/usr/bin/env python3

from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from pybiomart import Dataset


INPUT_H5AD = Path(
    "/data/mouse_vis_cortex_pseudobulk/"
    "Developing_Mouse_Visual_Cortex_library_cluster_pseudobulk.h5ad"
)

OUTPUT_CSV = Path(
    "/results/mouse_gene_summary_statistics.csv"
)

GENE_BLOCK_SIZE = 256
MIN_NONZERO_N = 3


def first_existing(columns, candidates):
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def infer_gene_symbols(adata):
    for field in (
        "gene_symbol",
        "gene_symbols",
        "symbol",
        "gene_name",
    ):
        if field in adata.var.columns:
            values = adata.var[field].astype("string")
            fallback = pd.Series(
                adata.var_names.astype(str),
                index=adata.var_names,
                dtype="string",
            )

            return (
                values
                .replace("", pd.NA)
                .fillna(fallback)
                .astype(str)
                .to_numpy()
            )

    return adata.var_names.astype(str).to_numpy()


def fetch_gene_types(gene_symbols):
    dataset = Dataset(
        name="mmusculus_gene_ensembl",
        host="http://www.ensembl.org",
    )

    biomart = dataset.query(
        attributes=[
            "external_gene_name",
            "gene_biotype",
        ]
    )

    if biomart.shape[1] != 2:
        raise RuntimeError(
            "BioMart did not return the expected gene-symbol and "
            "gene-biotype columns."
        )

    biomart = biomart.iloc[:, :2].copy()
    biomart.columns = [
        "gene_symbol",
        "gene_type",
    ]

    biomart["gene_symbol"] = (
        biomart["gene_symbol"]
        .astype("string")
        .str.strip()
    )
    biomart["gene_type"] = (
        biomart["gene_type"]
        .astype("string")
        .str.strip()
    )

    biomart = biomart.dropna(
        subset=["gene_symbol", "gene_type"]
    )
    biomart = biomart[
        (biomart["gene_symbol"] != "") &
        (biomart["gene_type"] != "")
    ]
    biomart = biomart.drop_duplicates(
        subset=["gene_symbol"],
        keep="first",
    )

    gene_type_lookup = dict(
        zip(
            biomart["gene_symbol"].astype(str),
            biomart["gene_type"].astype(str),
        )
    )

    gene_types = np.asarray(
        [
            gene_type_lookup.get(str(symbol), "unknown")
            for symbol in gene_symbols
        ],
        dtype=object,
    )

    print(
        "Matched BioMart gene_types for "
        f"{np.sum(gene_types != 'unknown'):,} of "
        f"{len(gene_types):,} genes."
    )

    return gene_types


def parse_developmental_age(value):
    if value is None or pd.isna(value):
        return np.nan, None

    text = str(value).strip().upper()

    embryonic_match = re.search(
        r"\bE\s*([0-9]+(?:\.[0-9]+)?)",
        text,
    )

    if embryonic_match:
        return float(embryonic_match.group(1)), "embryonic"

    postnatal_match = re.search(
        r"\bP\s*([0-9]+(?:\.[0-9]+)?)",
        text,
    )

    if postnatal_match:
        return float(postnatal_match.group(1)), "postnatal"

    return np.nan, None


def safe_log2fc(numerator, denominator):
    return np.log2(
        (np.asarray(numerator, dtype=np.float64) + 1.0) /
        (np.asarray(denominator, dtype=np.float64) + 1.0)
    )


def positive_mean(values, minimum_n=MIN_NONZERO_N):
    positive = values > 0
    positive_n = positive.sum(axis=0)

    positive_sum = np.where(
        positive,
        values,
        0.0,
    ).sum(axis=0)

    result = np.full(
        values.shape[1],
        np.nan,
        dtype=np.float64,
    )

    eligible = positive_n >= minimum_n

    result[eligible] = (
        positive_sum[eligible] /
        positive_n[eligible]
    )

    return result


def grouped_means(
    values,
    groups,
    nonzero_only=False,
    minimum_nonzero_n=MIN_NONZERO_N,
):
    groups = np.asarray(groups, dtype=object)

    valid_group = np.asarray(
        [
            value is not None
            and not pd.isna(value)
            and str(value) != ""
            for value in groups
        ],
        dtype=bool,
    )

    labels = pd.unique(
        groups[valid_group].astype(str)
    ).tolist()

    means = []

    for label in labels:
        group_mask = (
            valid_group &
            (groups.astype(str) == label)
        )

        group_values = values[group_mask, :]

        if nonzero_only:
            group_mean = positive_mean(
                group_values,
                minimum_n=minimum_nonzero_n,
            )
        else:
            group_mean = np.mean(
                group_values,
                axis=0,
            )

        means.append(group_mean)

    if not means:
        return [], np.empty(
            (0, values.shape[1]),
            dtype=np.float64,
        )

    return labels, np.vstack(means)


def maximum_group_statistics(labels, means):
    n_genes = means.shape[1]

    maximum_labels = np.empty(
        n_genes,
        dtype=object,
    )

    maximum_labels[:] = None

    log2fc = np.full(
        n_genes,
        np.nan,
        dtype=np.float64,
    )

    for gene_index in range(n_genes):
        gene_means = means[:, gene_index]
        valid = np.isfinite(gene_means)

        if not np.any(valid):
            continue

        valid_indices = np.flatnonzero(valid)

        maximum_index = valid_indices[
            np.argmax(gene_means[valid])
        ]

        maximum_value = gene_means[maximum_index]
        mean_across_groups = np.mean(
            gene_means[valid]
        )

        maximum_labels[gene_index] = labels[
            maximum_index
        ]

        log2fc[gene_index] = safe_log2fc(
            maximum_value,
            mean_across_groups,
        )

    return maximum_labels, log2fc

def weighted_fraction(detected_values, cell_totals):
    cell_totals = np.asarray(cell_totals, dtype=np.float64)
    denominator = np.sum(cell_totals)
    if not np.isfinite(denominator) or denominator <= 0:
        return np.full(detected_values.shape[1], np.nan, dtype=np.float64)
    return np.sum(detected_values, axis=0, dtype=np.float64) / denominator

def grouped_detection_fractions(detected_values, cell_totals, groups):
    groups = np.asarray(groups, dtype=object)
    cell_totals = np.asarray(cell_totals, dtype=np.float64)
    valid_group = np.asarray(
        [
            value is not None
            and not pd.isna(value)
            and str(value) != ""
            for value in groups
        ],
        dtype=bool,
    )
    labels = pd.unique(groups[valid_group].astype(str)).tolist()
    fractions = []
    for label in labels:
        group_mask = valid_group & (groups.astype(str) == label)
        fractions.append(
            weighted_fraction(
                detected_values[group_mask, :],
                cell_totals[group_mask],
            )
        )
    if not fractions:
        return [], np.empty((0, detected_values.shape[1]), dtype=np.float64)
    return labels, np.vstack(fractions)

def maximum_fraction_statistics(labels, fractions):
    n_genes = fractions.shape[1]
    maximum_labels = np.empty(n_genes, dtype=object)
    maximum_labels[:] = None
    fraction_difference = np.full(n_genes, np.nan, dtype=np.float64)
    for gene_index in range(n_genes):
        gene_fractions = fractions[:, gene_index]
        valid = np.isfinite(gene_fractions)
        if not np.any(valid):
            continue
        valid_indices = np.flatnonzero(valid)
        maximum_index = valid_indices[np.argmax(gene_fractions[valid])]
        maximum_labels[gene_index] = labels[maximum_index]
        fraction_difference[gene_index] = (
            gene_fractions[maximum_index] - np.mean(gene_fractions[valid])
        )
    return maximum_labels, fraction_difference

if not INPUT_H5AD.exists():
    candidates = list(
        Path("/data").glob("**/*pseudobulk*.h5ad")
    )

    if len(candidates) == 1:
        INPUT_H5AD = candidates[0]
    else:
        raise FileNotFoundError(
            f"Input H5AD not found: {INPUT_H5AD}"
        )


print(f"Reading: {INPUT_H5AD}")

adata = ad.read_h5ad(INPUT_H5AD)

if adata.n_obs == 0 or adata.n_vars == 0:
    raise ValueError(
        f"Input H5AD is empty: {adata.shape}"
    )


age_field = first_existing(
    adata.obs.columns,
    (
        "donor_age",
        "age",
        "age_label",
        "developmental_age",
        "age_in_days",
        "age_days",
    ),
)

roi_field = first_existing(
    adata.obs.columns,
    (
        "region_of_interest_label",
        "Brain Region",
        "brain_region",
        "region_label",
        "structure",
        "ROI",
    ),
)

subclass_field = first_existing(
    adata.obs.columns,
    (
        "subclass",
        "subclass_label",
        "Subclass",
    ),
)

if age_field is None:
    raise KeyError(
        "Could not identify the developmental-age column."
    )

if roi_field is None:
    raise KeyError(
        "Could not identify the brain-region column."
    )

if subclass_field is None:
    raise KeyError(
        "Could not identify the subclass column."
    )


parsed_ages = [
    parse_developmental_age(value)
    for value in adata.obs[age_field]
]

age_numbers = np.asarray(
    [value[0] for value in parsed_ages],
    dtype=np.float64,
)

age_stages = np.asarray(
    [value[1] for value in parsed_ages],
    dtype=object,
)

embryonic_mask = (
    age_stages == "embryonic"
)

p40_or_older_mask = (
    (age_stages == "postnatal") &
    np.isfinite(age_numbers) &
    (age_numbers >= 40)
)

if embryonic_mask.sum() == 0:
    raise ValueError(
        f"No embryonic ages were found in '{age_field}'."
    )

if p40_or_older_mask.sum() == 0:
    raise ValueError(
        f"No P40-or-older observations were found in '{age_field}'."
    )


gene_symbols = infer_gene_symbols(adata)
gene_types = fetch_gene_types(gene_symbols)

roi_groups = (
    adata.obs[roi_field]
    .astype("string")
    .to_numpy()
)

subclass_groups = (
    adata.obs[subclass_field]
    .astype("string")
    .to_numpy()
)


if sparse.issparse(adata.X):
    expression_matrix = adata.X.tocsc(copy=False)
else:
    expression_matrix = np.asarray(adata.X)

detection_layer = "number_of_cells_expressing"
if detection_layer not in adata.layers:
    raise KeyError(
        f"Input H5AD must contain layers['{detection_layer}']."
    )
if "number_of_cells" not in adata.obs.columns:
    raise KeyError("Input H5AD obs must contain 'number_of_cells'.")

detection_matrix = adata.layers[detection_layer]
if sparse.issparse(detection_matrix):
    detection_matrix = detection_matrix.tocsc(copy=False)
else:
    detection_matrix = np.asarray(detection_matrix)
if detection_matrix.shape != adata.shape:
    raise ValueError("Detection-count layer shape does not match X.")

number_of_cells = pd.to_numeric(
    adata.obs["number_of_cells"],
    errors="coerce",
).to_numpy(dtype=np.float64)
if np.any(~np.isfinite(number_of_cells)) or np.any(number_of_cells <= 0):
    raise ValueError("number_of_cells must contain finite positive values.")

library_totals = np.asarray(
    expression_matrix.sum(axis=1)
).reshape(-1).astype(np.float64)

cpm_scaling_factor = np.zeros(
    adata.n_obs,
    dtype=np.float64,
)

positive_library_total = library_totals > 0

cpm_scaling_factor[positive_library_total] = (
    1_000_000.0 /
    library_totals[positive_library_total]
)


statistics_blocks = []


for block_start in range(
    0,
    adata.n_vars,
    GENE_BLOCK_SIZE,
):
    block_stop = min(
        block_start + GENE_BLOCK_SIZE,
        adata.n_vars,
    )

    if sparse.issparse(expression_matrix):
        raw_counts = (
            expression_matrix[
                :,
                block_start:block_stop,
            ]
            .toarray()
            .astype(np.float64, copy=False)
        )
    else:
        raw_counts = np.asarray(
            expression_matrix[
                :,
                block_start:block_stop,
            ],
            dtype=np.float64,
        )

    if sparse.issparse(detection_matrix):
        detected_counts = (
            detection_matrix[
                :,
                block_start:block_stop,
            ]
            .toarray()
            .astype(np.float64, copy=False)
        )
    else:
        detected_counts = np.asarray(
            detection_matrix[
                :,
                block_start:block_stop,
            ],
            dtype=np.float64,
        )
    if np.any(detected_counts < 0):
        raise ValueError("Detection counts contain negative values.")
    if np.any(detected_counts > number_of_cells[:, np.newaxis]):
        raise ValueError(
            "Detection counts exceed number_of_cells for at least one row."
        )

    cpm = (
        raw_counts *
        cpm_scaling_factor[:, np.newaxis]
    )


    mean_cpm = np.mean(
        cpm,
        axis=0,
    )

    ln_mean_cpm = np.log1p(
        mean_cpm
    )


    embryonic_mean = np.mean(
        cpm[embryonic_mask, :],
        axis=0,
    )

    p40_or_older_mean = np.mean(
        cpm[p40_or_older_mask, :],
        axis=0,
    )

    age_log2fc = safe_log2fc(
        p40_or_older_mean,
        embryonic_mean,
    )


    overall_fraction_expressing = weighted_fraction(
        detected_counts,
        number_of_cells,
    )
    embryonic_fraction_expressing = weighted_fraction(
        detected_counts[embryonic_mask, :],
        number_of_cells[embryonic_mask],
    )
    p40_or_older_fraction_expressing = weighted_fraction(
        detected_counts[p40_or_older_mask, :],
        number_of_cells[p40_or_older_mask],
    )
    age_fraction_expressing_difference = (
        p40_or_older_fraction_expressing -
        embryonic_fraction_expressing
    )
    roi_labels, roi_means = grouped_means(
        cpm,
        roi_groups,
        nonzero_only=False,
    )

    max_roi, roi_log2fc = (
        maximum_group_statistics(
            roi_labels,
            roi_means,
        )
    )


    roi_fraction_labels, roi_fractions = grouped_detection_fractions(
        detected_counts,
        number_of_cells,
        roi_groups,
    )
    max_fraction_roi, fraction_difference_roi = (
        maximum_fraction_statistics(
            roi_fraction_labels,
            roi_fractions,
        )
    )
    subclass_labels, subclass_means = (
        grouped_means(
            cpm,
            subclass_groups,
            nonzero_only=False,
        )
    )

    max_subclass, subclass_log2fc = (
        maximum_group_statistics(
            subclass_labels,
            subclass_means,
        )
    )


    subclass_fraction_labels, subclass_fractions = (
        grouped_detection_fractions(
            detected_counts,
            number_of_cells,
            subclass_groups,
        )
    )
    max_fraction_subclass, fraction_difference_subclass = (
        maximum_fraction_statistics(
            subclass_fraction_labels,
            subclass_fractions,
        )
    )
    statistics_blocks.append(
        pd.DataFrame(
            {
                "gene_symbol": (
                    gene_symbols[
                        block_start:block_stop
                    ]
                ),
                "ln(mean_CPM+1)": (
                    ln_mean_cpm
                ),
                "age_log2fc": (
                    age_log2fc
                ),
                "overall_fraction_expressing": (
                    overall_fraction_expressing
                ),
                "age_fraction_expressing_difference": (
                    age_fraction_expressing_difference
                ),
                "max_ROI": (
                    max_roi
                ),
                "log2fc_ROI": (
                    roi_log2fc
                ),
                "max_fraction_ROI": (
                    max_fraction_roi
                ),
                "fraction_difference_ROI": (
                    fraction_difference_roi
                ),
                "max_subclass": (
                    max_subclass
                ),
                "log2fc_subclass": (
                    subclass_log2fc
                ),
                "max_fraction_subclass": (
                    max_fraction_subclass
                ),
                "fraction_difference_subclass": (
                    fraction_difference_subclass
                ),
                "gene_type": (
                    gene_types[
                        block_start:block_stop
                    ]
                ),
            }
        )
    )

    print(
        "Completed genes "
        f"{block_start + 1:,}-"
        f"{block_stop:,} of "
        f"{adata.n_vars:,}"
    )


statistics = pd.concat(
    statistics_blocks,
    ignore_index=True,
)


expected_columns = [
    "gene_symbol",
    "ln(mean_CPM+1)",
    "age_log2fc",
    "overall_fraction_expressing",
    "age_fraction_expressing_difference",
    "max_ROI",
    "log2fc_ROI",
    "max_fraction_ROI",
    "fraction_difference_ROI",
    "max_subclass",
    "log2fc_subclass",
    "max_fraction_subclass",
    "fraction_difference_subclass",
    "gene_type",
]

statistics = statistics[
    expected_columns
]


if statistics.shape[0] != adata.n_vars:
    raise RuntimeError(
        "The statistics table does not contain one row per gene."
    )

if statistics["gene_symbol"].isna().any():
    raise RuntimeError(
        "The statistics table contains missing gene symbols."
    )


OUTPUT_CSV.parent.mkdir(
    parents=True,
    exist_ok=True,
)

statistics.to_csv(
    OUTPUT_CSV,
    index=False,
)


check = pd.read_csv(OUTPUT_CSV)

if list(check.columns) != expected_columns:
    raise RuntimeError(
        "Output column names changed during serialization."
    )

if check.shape[0] != adata.n_vars:
    raise RuntimeError(
        "Output row count changed during serialization."
    )


print(f"Saved updated statistics table: {OUTPUT_CSV}")
print(f"Rows: {statistics.shape,}")
print(f"Columns: {statistics.shape,}")