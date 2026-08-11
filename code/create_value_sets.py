#!/usr/bin/env python3

from pathlib import Path
from urllib.request import urlopen

import pandas as pd


TAXONOMY_URL = (
    "https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/"
    "metadata/Developing-Mouse-Vis-Cortex-taxonomy/20260331/"
    "cluster_annotation_term.csv"
)

DATASET_VALUE_SETS_URL = (
    "https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/"
    "metadata/Developing-Mouse-Vis-Cortex-10X/20260131/"
    "value_sets.csv"
)

OUTPUT_PATH = Path(
    "/results/mouse_vis_cortex_parquet/metadata/value_sets.csv"
)


def read_public_csv(url):
    with urlopen(url) as response:
        return pd.read_csv(response)


def require_columns(data, required, source_name):
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(
            f"{source_name} is missing required columns: "
            + ", ".join(missing)
        )


print("Reading taxonomy annotation terms...")
taxonomy = read_public_csv(TAXONOMY_URL)

require_columns(
    taxonomy,
    [
        "name",
        "cluster_annotation_term_set_name",
        "color_hex_triplet",
        "term_order",
    ],
    "cluster_annotation_term.csv",
)

taxonomy_value_sets = taxonomy[
    [
        "name",
        "cluster_annotation_term_set_name",
        "color_hex_triplet",
        "term_order",
    ]
].rename(
    columns={
        "name": "label",
        "cluster_annotation_term_set_name": "field",
        "term_order": "order",
    }
)


print("Reading dataset value sets...")
dataset_value_sets = read_public_csv(DATASET_VALUE_SETS_URL)

require_columns(
    dataset_value_sets,
    [
        "label",
        "field",
        "color_hex_triplet",
        "order",
    ],
    "value_sets.csv",
)

dataset_value_sets = dataset_value_sets[
    [
        "label",
        "field",
        "color_hex_triplet",
        "order",
    ]
]


value_sets = pd.concat(
    [
        taxonomy_value_sets,
        dataset_value_sets,
    ],
    ignore_index=True,
)

value_sets["label"] = value_sets["label"].astype("string").str.strip()
value_sets["field"] = value_sets["field"].astype("string").str.strip()
value_sets["color_hex_triplet"] = (
    value_sets["color_hex_triplet"]
    .astype("string")
    .str.strip()
)
value_sets["order"] = pd.to_numeric(
    value_sets["order"],
    errors="coerce",
)

value_sets = value_sets.loc[
    value_sets["label"].notna()
    & value_sets["label"].ne("")
    & value_sets["field"].notna()
    & value_sets["field"].ne("")
].copy()

value_sets = value_sets[
    [
        "label",
        "field",
        "color_hex_triplet",
        "order",
    ]
]

OUTPUT_PATH.parent.mkdir(
    parents=True,
    exist_ok=True,
)

value_sets.to_csv(
    OUTPUT_PATH,
    index=False,
)

check = pd.read_csv(OUTPUT_PATH)

expected_columns = [
    "label",
    "field",
    "color_hex_triplet",
    "order",
]

if list(check.columns) != expected_columns:
    raise RuntimeError(
        "Output columns changed during serialization."
    )

if len(check) != len(value_sets):
    raise RuntimeError(
        "Output row count changed during serialization."
    )

print(f"Saved: {OUTPUT_PATH}")
print(f"Rows: {len(value_sets):,}")
print(f"Fields: {value_sets['field'].nunique():,}")
