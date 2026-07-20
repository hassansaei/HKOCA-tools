# Snapseed Annotation

Hierarchical cell-type annotation of QC-filtered, harmonized kidney organoid single-cell RNA-seq datasets using [Snapseed](https://github.com/theislab/snapseed), with curated marker hierarchies, optional manual overrides, and downstream cross-dataset reporting.

This directory contains the primary analysis notebook, a portable configuration file, and the marker gene hierarchy used by Snapseed.

---

## Contents

| File | Description |
|------|-------------|
| [`snapseed_annotation.ipynb`](./snapseed_annotation.ipynb) | End-to-end annotation notebook (no stored outputs; suitable for version control) |
| [`snapseed_annotation_config.yaml`](./snapseed_annotation_config.yaml) | Pipeline configuration and parameters |
| [`snapseed_markers_v4.yaml`](./snapseed_markers_v4.yaml) | Hierarchical Snapseed marker gene dictionary (Level 1 to Level 3) |
| [`README.md`](./README.md) | This documentation |

Runtime outputs (not committed by default):

```
├── clustered objects/     # intermediate clustered .h5ad
├── annotated objects/     # final annotated .h5ad
└── summary_figures/       # figures and CSV tables
```

---

## Overview

The workflow performs **per-dataset** annotation of kidney organoid scRNA-seq objects that have already been quality-controlled and gene-harmonized. For each study:

1. Load a QC-filtered AnnData object (`.h5ad`).
2. Normalize / log-transform when needed; recover unscaled expression if a scaled matrix is detected.
3. Select highly variable genes, compute PCA, neighbors, Leiden clustering, and UMAP.
4. Rank marker genes (Wilcoxon) and classify clusters with **Snapseed** against the hierarchical marker dictionary.
5. Optionally apply **manual label overrides** and QC-driven cluster removal.
6. Persist clustered and annotated objects.
7. Aggregate all annotated datasets into summary figures, protocol-balance plots, and adult/fetal pseudobulk correlations.

Annotation levels written to `adata.obs`:

| Column | Meaning |
|--------|---------|
| `Level_1` | Broad lineage (for example Nephron, Stromal_cells, Off_target_cells) |
| `Level_2` | Intermediate subtype (for example Proximal_tubule, Stromal_progenitors) |
| `Level_3` | Finest Snapseed assignment when available |
| `Level_3_latest` | Deepest non-missing label with fallback to Level 2, then Level 1, then cluster ID |

---

## What is Snapseed?

[Snapseed](https://github.com/devsystemslab/snapseed) is a marker-based annotation tool for single-cell RNA-seq data. It assigns labels to **predefined groups of cells** (typically Leiden clusters) by matching their expression profiles against **manually curated marker gene sets**. It was developed to annotate large organoid atlases quickly, including GPU-accelerated scoring for high-throughput use.

Unlike reference-mapping methods (which transfer labels from an annotated atlas), Snapseed does not require a reference embedding. Knowledge enters through the marker dictionary: experts define which genes characterize each cell type, and Snapseed scores every cluster against those definitions.

### Flat vs hierarchical annotation

Snapseed supports two modes:

| Mode | Function | Use when |
|------|----------|----------|
| Flat | `annotate` | Cell types sit at one level with no nesting |
| Hierarchical | `annotate_hierarchy` | Cell types form a tree (lineage → subtype → fine state) |

This notebook uses **`annotate_hierarchy`**, because kidney organoid biology is naturally nested (for example Nephron → Podocytes → Mature_podocytes).

### How hierarchical annotation works

1. **Marker tree.** Markers are supplied as a nested YAML dictionary (`snapseed_markers_v4.yaml`). Each node lists `marker_genes` and optional `subtypes`.
2. **Cluster-level scoring.** For each Leiden cluster, Snapseed evaluates how well the cluster expresses the markers of each candidate cell type at the current tree level. Scoring favors genes that are enriched in the cluster relative to other cells.
3. **Top-down assignment.** At Level 1, each cluster is assigned the best-matching broad lineage (Nephron, Stromal_cells, Off_target_cells, and so on).
4. **Descent into subtypes.** Given that Level 1 call, Snapseed only considers children under that node for Level 2, then Level 3. A cluster labeled Nephron is scored against Nephron subtypes (NPC, podocyte, tubule, ...), not against stromal or endothelial subtypes.
5. **Outputs.** The call returns:
   - **assignments**: best label per cluster at each depth
   - **metrics**: confidence / enrichment scores used for the call

Clusters that lack a confident finer label may stop early in the tree. This notebook then builds `Level_3_latest` by taking the deepest available assignment, falling back to Level 2, Level 1, or the cluster ID.

### What Snapseed does not do

- It does **not** invent cell types absent from the marker YAML.
- It does **not** replace expert review: borderline clusters, doublets, and off-target contaminants often need manual overrides.
- It does **not** recompute clustering; labels are attached to an existing `groupby` column (here, Leiden clusters).

In this pipeline, Scanpy prepares clusters and differential expression; Snapseed maps those clusters onto the curated kidney organoid hierarchy; optional curation corrects known edge cases.

---

## Requirements

### Software

| Package | Role |
|---------|------|
| Python ≥ 3.10 | Runtime |
| [scanpy](https://scanpy.readthedocs.io/) | Preprocessing, clustering, DE |
| [anndata](https://anndata.readthedocs.io/) | Data containers |
| [snapseed](https://github.com/theislab/snapseed) | Hierarchical marker-based annotation |
| pandas, numpy, scipy | Numerics and tables |
| matplotlib, plotly | Visualization |
| PyYAML | Config and marker loading |

Optional (downstream correlation / protocol scatter cells):

- `scipy.stats.spearmanr`
- Adult and fetal reference `.h5ad` files (configured in `snapseed_annotation_config.yaml`)

### Hardware

Per-dataset cells load one object at a time. Memory requirements scale with the largest organoid dataset (typically multi-GB RAM for full PCA + Leiden + UMAP). Downstream summary cells load all annotated objects; reserve sufficient memory for that stage.

---

## How to run

1. Install the dependencies listed above and open [`snapseed_annotation.ipynb`](./snapseed_annotation.ipynb).
2. Run the setup cell first (loads `snapseed_annotation_config.yaml` and `snapseed_markers_v4.yaml`).
3. Run the pipeline definition cell (`process_and_annotate_dataset`).
4. Run per-dataset annotation cells (one study per cell). Adjust Leiden resolution and `curated_overrides` as needed.
5. Optionally run the cross-dataset summary, protocol composition, and adult/fetal correlation sections.

---

## Configuration (`snapseed_annotation_config.yaml`)

Shared pipeline hyperparameters live in this file so the notebook remains portable across environments.

### Schema

```yaml
hkoca_tools_root: null   # optional; null = auto-detect repository root

paths:
  qc_input_dir: ...
  marker_yaml: ...            # snapseed_markers_v4.yaml
  annotated_output_dir: ...
  clustered_output_dir: ...
  summary_figures_dir: ...
  adult_ref_h5ad: ...
  fetal_ref_h5ad: ...
  adult_markers_dir: ...
  fetal_markers_dir: ...

parameters:
  seed: 42
  scanpy_verbosity: 0
  default_cluster_key: leiden_clusters
  hvg_n_top_genes: 3000
  pca_n_comps: 50
  neighbors_n_neighbors: 20
  neighbors_n_pcs: 30
  scale_max_value: 10
  normalize_target_sum: 10000
  dpi: 300
```

### Parameters used by the pipeline

| Key | Default | Used for |
|-----|---------|----------|
| `seed` | `42` | `random`, NumPy, Scanpy PCA / Leiden / UMAP |
| `hvg_n_top_genes` | `3000` | Highly variable gene selection |
| `pca_n_comps` | `50` | PCA on HVG-scaled matrix |
| `neighbors_n_neighbors` | `20` | kNN graph |
| `neighbors_n_pcs` | `30` | PCs retained for neighbors |
| `scale_max_value` | `10` | Gene scaling clip for PCA |
| `normalize_target_sum` | `10000` | Library-size normalization when raw counts are detected |
| `dpi` | `300` | Figure and savefig resolution |

---

## Marker hierarchy (`snapseed_markers_v4.yaml`)

Snapseed annotation is driven by this YAML dictionary. Structure matches Snapseed's hierarchical input: each node has `marker_genes` and optional nested `subtypes`.

### Top-level lineages (Level 1)

| Lineage | Role in kidney organoids |
|---------|---------------------------|
| **Nephron** | NPC to podocyte / tubule lineages |
| **Ureteric_Epithelium** | Ureteric bud / collecting duct |
| **Stromal_cells** | Stromal progenitors, fibroblasts, mesangium, pericytes |
| **Endothelial_cells** | Early and late endothelium |
| **Off_target_cells** | Non-renal contaminants (neural, muscle, glial, endoderm) |

### Example nesting

```yaml
Nephron:
  marker_genes: [NPHS1, NPHS2, EPCAM, ...]
  subtypes:
    Podocytes:
      marker_genes: [NPHS1, NPHS2, MAFB, ...]
      subtypes:
        Precursor_podocytes:
          marker_genes: [OLFM3, MAFB, CCN2, NPHS1]
        Mature_podocytes:
          marker_genes: [PTPRO, CLIC5, DDN, ...]
```

Editing markers:

1. Update gene symbols in [`snapseed_markers_v4.yaml`](./snapseed_markers_v4.yaml).
2. Re-run the setup cell to reload `MARKER_HIERARCHY`.
3. Re-run annotation for affected datasets.

Keep gene symbols consistent with the harmonized feature names in the input `.h5ad` objects (typically HGNC symbols after gene harmonization).

---

## Notebook structure

Sections follow the order of execution.

### 1. Setup

Loads `snapseed_annotation_config.yaml`, sets seeds and figure DPI, imports libraries, loads `snapseed_markers_v4.yaml` into `MARKER_HIERARCHY`, and defines helpers:

- `display_ranked_markers` - top DE genes scored as  
  \(\mathrm{score} = \mathrm{logFC} \times (\mathrm{pct\_in})^2 \times (1 - \mathrm{pct\_out})\)
- `_plot_umap_with_labels` / `_plot_annotation_umaps_with_labels` - cluster and annotation UMAPs

### 2. Pipeline (`process_and_annotate_dataset`)

Core function signature (conceptually):

```python
process_and_annotate_dataset(
    file_path,
    marker_dict,
    resolution=1.0,
    cluster_key="leiden_clusters",
    manual_annotations=None,
    output_dir=OUTPUT_DIR,
    clustered_dir=CLUSTERED_DIR,
)
```

**Processing steps**

1. Read `.h5ad`; set gene names from `var["features"]` when present.
2. Guard against pre-scaled matrices (negative values): restore from `.raw` or a counts/raw/normalized layer when available.
3. Normalize + `log1p` if max expression suggests raw counts; otherwise skip.
4. HVG selection (config `hvg_n_top_genes`).
5. Scale HVGs on a copy (preserve unscaled `.X` for DE / Snapseed), PCA, neighbors, Leiden, UMAP.
6. Sanitize reserved AnnData column names (`_index`) that break h5ad writes.
7. Write clustered object.
8. Wilcoxon DE; Snapseed `annotate_hierarchy`.
9. Apply optional `manual_annotations` (cluster ID to Level_1 / Level_2 / Level_3).
10. Map assignments to `obs`, build `Level_3_latest`, plot annotation UMAPs.
11. Write annotated object (`*_annotated.h5ad`).

### 3. Per-dataset annotation

One code cell per study. Typical pattern:

```python
dataset_name = "<study>_harmonized_singlets_filtered.h5ad"
target_file = str(BASE_DIR / dataset_name)
chosen_resolution = 0.4
dynamic_cluster_key = f"leiden_res_{chosen_resolution}"
curated_overrides = None  # or {"3": {"Level_1": "...", "Level_2": "...", "Level_3": "..."}}

adata_... = process_and_annotate_dataset(
    file_path=target_file,
    marker_dict=MARKER_HIERARCHY,
    resolution=chosen_resolution,
    cluster_key=dynamic_cluster_key,
    manual_annotations=curated_overrides,
)
```

**QC-aware cells** (outlier cluster removal without recomputing neighbors/UMAP) are marked with markdown notes in the notebook. After subsetting a bad cluster, Snapseed is re-run and manual overrides are reapplied on object-typed `obs` columns before converting back to categorical.

### 4. Cross-dataset summary

Loads all `*_annotated.h5ad` files and produces:

- Level 1 composition bar charts
- Level 2 consistency / detection matrices
- Summary statistics panels
- CSV tables (composition %, cell counts, detection rates, dataset metadata)

### 5. Protocol composition

Uses `diff_protocol` (or equivalent) metadata to compare nephron vs stromal balance across differentiation protocols (for example Morizane, Takasato, hybrids).

### 6. Adult / fetal pseudobulk correlation

Builds mean-expression pseudobulk profiles for organoid Level 3 labels and adult/fetal reference cell types, then computes Spearman correlations on shared genes (optionally restricted to highly variable genes). Reference files are set in the config.

---

## Inputs and outputs

### Inputs

| Source | Description |
|--------|-------------|
| QC-filtered `.h5ad` files | Harmonized single-dataset objects |
| `snapseed_markers_v4.yaml` | Marker hierarchy |
| Adult / fetal references | Optional; required only for the correlation section |

Expected gene metadata: `adata.var["features"]` holding gene symbols when present. Expression may be counts or already log-normalized; the pipeline branches accordingly.

### Outputs

| Artifact | Description |
|----------|-------------|
| `*_clustered.h5ad` | Intermediate clustered objects |
| `*_annotated.h5ad` | Final annotated objects |
| `*_cleaned_annotated.h5ad` | Annotated objects after QC cluster removal (selected studies) |
| Figures (PDF/PNG) | Summary and protocol plots |
| Tables (CSV) | Composition, counts, detection, metadata |

Annotated objects store Snapseed results in:

- `adata.uns["snapseed_assignments"]`
- `adata.uns["snapseed_metrics"]`
- `adata.obs["Level_1" | "Level_2" | "Level_3" | "Level_3_latest"]`

---

## Reproducibility

- Global `SEED` (default `42`) is applied to Python `random`, NumPy, and Scanpy (`pca`, `leiden`, `umap`).
- Leiden resolution is set **per dataset** via `chosen_resolution` and recorded in the cluster key name (`leiden_res_<resolution>`).
- Manual overrides and QC cluster IDs are explicit in each dataset cell so curation decisions remain auditable in the notebook history.
- Figure DPI defaults to `300` via the config parameter `dpi`.

---

## Design notes

**Why a separate config file?**  
Configuration stays outside the notebook logic, so hyperparameters and environment settings can change without editing analysis code.

**Why keep markers in YAML?**  
The hierarchy is shared, reviewable, and editable without touching Python. Snapseed consumes the same nested structure.

**Why intermediate clustered objects?**  
Clustering is expensive. Saving clustered AnnData allows re-annotation or marker experimentation without recomputing PCA / neighbors / UMAP when only labels change (QC-removal cells already preserve the embedding).

**Manual overrides**  
Snapseed assignments are automated. Curated dictionaries correct known biology (for example off-target neural clusters, stroma vs nephron borderline labels) after expert review of markers and UMAPs.

---

## Troubleshooting

| Symptom | Likely cause | Action |
|---------|--------------|--------|
| Config / marker file not found | Wrong working directory, or repository root not detected | Open the notebook from this folder; set `hkoca_tools_root` in the config if needed |
| `FileNotFoundError` for a dataset | Filename mismatch in the QC input directory | Align `dataset_name` with the available `.h5ad` files |
| Negative values / scaled matrix errors | Object was scaled without `.raw` or layers | Ensure QC export keeps unscaled counts in `.X`, `.raw`, or a named layer |
| Empty Snapseed levels | Markers missing from `var_names` | Check harmonized gene symbols vs `snapseed_markers_v4.yaml` |
| Summary cell finds no files | Annotation cells not run yet | Re-run annotation cells first |
| Categorical assignment errors after overrides | Labels edited while still categorical | Pipeline / QC cells cast to `object` before overrides, then back to `category` |


