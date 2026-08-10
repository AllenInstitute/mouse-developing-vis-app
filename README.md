# mouse-aging-app

Code and reproducible workflows for building a pseudobulk gene-expression resource and interactive viewer for the Developing Mouse Visual Cortex dataset from the Allen Brain Cell Atlas.

## Current workflow

The first stage creates a pseudobulk AnnData (`.h5ad`) file from the raw single-cell count matrix. Expression values are summed for each combination of:

- `library_label`
- `cluster_alias`

Each pseudobulk observation retains representative library, donor, and taxonomy metadata and includes `number_of_cells`, the number of cells contributing to that observation.

## Repository structure

```text
.
├── code/
│   ├── create_dev_mouse_vis_cortex_pseudobulk.py
│   └── run
├── data/
│   └── .gitkeep
├── .gitignore
├── README.md
└── requirements.txt