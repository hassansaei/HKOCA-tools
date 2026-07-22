# HKOCA Level 3 Reannotation

Documentation for the high-resolution reannotation pipeline in `Benchmark/`. This step runs on the scPoli-integrated atlas to refine Level_3 labels via Leiden clustering at resolution 10.0, propagate Level_2 and Level_1 labels, and validate cell type structure with marker genes.


---

## Table of contents

1. [Motivation](#motivation)
2. [Overview](#overview)
3. [Directory layout](#directory-layout)
4. [Workflow](#workflow)
5. [Script reference](#script-reference)
6. [Configuration](#configuration)
7. [Outputs](#outputs)
8. [Usage](#usage)
9. [Requirements](#requirements)
10. [Troubleshooting](#troubleshooting)

---

## Motivation

After scPoli integration, cell type labels from the pre-integration annotation may not fully reflect the integrated embedding structure. High-resolution clustering on `X_emb` identifies micro-clusters whose dominant Level_3 label can differ from the original assignment. Majority-vote relabeling per micro-cluster produces integrated annotation columns (`Level_3_Integrated`, `Level_2_Integrated`, `Level_1_Integrated`) that are consistent with the batch-corrected representation.

Marker gene analysis on a 100k-cell subset provides a QC check that integrated labels remain biologically coherent.

---

## Overview

```
reannotation_config.yaml
        |
        v
reannotate_res10.py --config reannotation_config.yaml
        |
        +-- Load integrated atlas
        +-- Neighbors + Leiden (resolution = 10, use_rep = X_emb)
        +-- Majority-vote Level_3 assignment per micro-cluster
        +-- Propagate Level_2 and Level_1 from Level_3 maps
        +-- Wilcoxon DEGs on 100k subset
        +-- Hierarchical marker matrix plot
        +-- Save reannotated .h5ad
```

**Prerequisite:** scPoli integration must be completed first. Input is `Master_Atlas_scPoli_Level2_Integrated_parameters.h5ad` with `X_emb` in `adata.obsm`.

---

## Directory layout

```
push/
├── reannotation_README.md       # This file
├── reannotation_config.yaml     # Central configuration
└── reannotate_res10.py          # Reannotation script
```

Default output paths point to `parameters/` and `parameters/reannotation/figures/`, matching the original layout.

---

## Workflow

### Step 1: Neighborhood graph and Leiden clustering

Build a kNN graph on the integrated embedding (`X_emb`, 30 neighbors) and run Leiden clustering at resolution 10.0. This produces many high-resolution micro-clusters stored in `leiden_res_10`.

### Step 2: Dominant cell type assignment

For each micro-cluster, assign the mode of `Level_3_latest` as the new Level_3 label. Build Level_3-to-Level_2 and Level_3-to-Level_1 maps from the original annotations (mode per Level_3), then propagate to produce `Level_2_Integrated` and `Level_1_Integrated`.

### Step 3: Marker gene refinement

Subsample 100,000 cells, run Wilcoxon rank-sum tests grouped by `Level_3_Integrated`, extract top-10 markers per cell type, compute a dendrogram, and render a hierarchical matrix plot.

### Step 4: Save

Write the full atlas with new integrated annotation columns to disk.

---

## Script reference

### `reannotate_res10.py`

| Aspect | Detail |
|--------|--------|
| Original | `parameters/reannotation/reannotate_res10.py` |
| Config | `reannotation_config.yaml` (default) or `--config /path/to/config.yaml` |
| Input | Integrated atlas with `X_emb`, `Level_3_latest`, `Level_2`, `Level_1` |
| RAM | Full atlas loaded into memory |

Differences from the original script:

- Paths and parameters moved to YAML config
- `--config` flag support added

All processing logic is unchanged.

---

## Configuration

### Config file: `reannotation_config.yaml`

| Section | Purpose |
|---------|---------|
| `paths` | Input integrated atlas, output h5ad, figure directory |
| `clustering` | Neighbors and Leiden parameters on `X_emb` |
| `annotation` | Source and output metadata column names |
| `marker_refinement` | Subset size, Wilcoxon test settings |
| `plotting` | Matrix plot appearance and filename |

### Full config reference

```yaml
paths:
  input_h5ad: .../Master_Atlas_scPoli_Level2_Integrated_parameters.h5ad
  output_h5ad: .../Master_Atlas_scPoli_Integrated_Reannotated.h5ad
  figure_dir: .../reannotation/figures

clustering:
  use_rep: X_emb
  n_neighbors: 30
  random_state: 42
  resolution: 10.0
  leiden_key: leiden_res_10

annotation:
  source_level3: Level_3_latest
  source_level2: Level_2
  source_level1: Level_1
  output_level3: Level_3_Integrated
  output_level2: Level_2_Integrated
  output_level1: Level_1_Integrated

marker_refinement:
  n_subset: 100000
  random_state: 42
  groupby: Level_3_Integrated
  method: wilcoxon
  n_top_genes: 10
  use_raw: false
  pts: true

plotting:
  dpi: 300
  facecolor: white
  cmap: Blues
  standard_scale: var
  title: Top 10 Marker Genes Hierarchical Clustering
  matrixplot_filename: Hierarchical_Marker_Genes_Refinement.png
```

### Example: write output to a new file

Do not overwrite existing source `.h5ad` files. Use a distinct output filename:

```yaml
paths:
  output_h5ad: /path/to/Master_Atlas_scPoli_Integrated_Reannotated_v2.h5ad
```

---

## Outputs

| File | Location | Description |
|------|----------|-------------|
| `Master_Atlas_scPoli_Integrated_Reannotated.h5ad` | `parameters/` (default) | Full atlas with integrated annotation columns |
| `Hierarchical_Marker_Genes_Refinement.png` | `parameters/reannotation/figures/` (default) | Top-10 marker genes per cell type with dendrogram |

### New observation columns

| Column | Description |
|--------|-------------|
| `leiden_res_10` | High-resolution Leiden cluster ID |
| `Level_3_Integrated` | Majority-vote Level_3 label per micro-cluster |
| `Level_2_Integrated` | Level_2 propagated from Level_3 map |
| `Level_1_Integrated` | Level_1 propagated from Level_3 map |

---

## Usage

```bash
cd scripts_benchmark/push

python reannotate_res10.py --config reannotation_config.yaml
```

With default config (omit `--config`):

```bash
python reannotate_res10.py
```

---

## Requirements

### Python

| Package | Required |
|---------|----------|
| scanpy | Yes |
| anndata | Yes |
| pandas | Yes |
| numpy | Yes |
| matplotlib | Yes |
| pyyaml | Yes |

### Hardware

| Resource | Recommendation |
|----------|----------------|
| GPU | Not required |
| RAM | 300+ GB (full atlas loaded into RAM) |

---

## Troubleshooting

**`X_emb` not found**

Run scPoli integration first. The script expects the integrated embedding in `adata.obsm['X_emb']`.

**Out of memory on load**

The full atlas is loaded into RAM. Run on a high-memory node or reduce atlas size upstream.

**High relabel rate**

Expected at resolution 10.0. Review `Hierarchical_Marker_Genes_Refinement.png` to confirm marker coherence across integrated labels.

**Missing source columns**

Ensure `Level_3_latest`, `Level_2`, and `Level_1` exist in `adata.obs` before running.

**Wilcoxon test fails on subset**

Confirm the subset retains multiple cells per `Level_3_Integrated` category. Very rare types may need a larger `n_subset` or manual review.
