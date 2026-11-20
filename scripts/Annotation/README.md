# Annotation pipeline

Annotation workflow to turn QC-filtered Seurat objects into a Snapseed-annotated kidney organoid atlas:

**RDS → h5ad → concatenate → Leiden (res 1.0) → Snapseed v5**

This folder ships the runnable scripts, their scheduler, configs, and the marker hierarchy. Large intermediate `.h5ad` files are written to paths defined in `config/` (edit those for your machine).

---

## What is Snapseed?

[Snapseed](https://github.com/) is a **marker-based hierarchical cell-type annotation** tool for single-cell data. Given:

1. an AnnData object with **cluster labels** (here: Leiden `leiden_res_1.0`), and  
2. a **YAML tree of cell types** with marker genes at each level,

it scores every cluster against each node in the hierarchy (default method: **AUROC** on marker expression) and assigns the best-matching labels.

In this project we use **markers v5** (`markers/snapseed_markers_v5.yaml`), a kidney-organoid hierarchy with up to three levels, for example:

| Level | Example |
|-------|---------|
| **Level_1** | Nephron, Stromal_cells, Endothelial_cells, Off_target_cells, … |
| **Level_2** | Podocytes, Proximal_tubule, S1/S2/S3, Neural_progenitor_cells, … |
| **Level_3** | Precursor_podocytes, Mature_podocytes, Dev_PT, PT, … |

`annotate_snapseed_v5.py` wraps `snapseed.annotate_hierarchy`, writes labels onto `adata.obs` (`Level_1`, `Level_2`, `Level_3`, `Level_3_latest`), and **flags clusters** whose overall annotation score is **&lt; 0.6** for manual review.

Clustering (normalize → HVG → PCA → neighbors → Leiden → UMAP) is included in the same script when `leiden_res_1.0` is missing, so you do not need a separate 600-line clustering job.

---

## Directory layout

```
annotation_pipeline/
├── README.md
├── scripts/          # pipeline code + scheduler
├── config/           # paths and parameters (YAML)
└── markers/          # snapseed_markers_v5.yaml
```

---

## Pipeline steps

| # | Step | Script | Role |
|---|------|--------|------|
| 1 | **RDS → h5ad** | `scripts/rds_to_h5ad.r` (invoked by the scheduler) | Convert each QC-filtered Seurat `.rds` to `.h5ad` via SeuratDisk (`SaveH5Seurat` → `Convert`). Assay5 RNA is cast to Assay before save. Optionally syncs h5ads into the concat input folder. |
| 2 | **Concatenate** | `scripts/concatenate_filtered_singlets.py` | Merge filtered singlet h5ads into one atlas AnnData (shared gene space, standardized `sample_id` / metadata columns). |
| 3 | **Cluster** | `scripts/annotate_snapseed_v5.py --cluster-only` | Leiden resolution **1.0** on the concatenated object; writes a clustered checkpoint h5ad. |
| 4 | **Annotate** | `scripts/annotate_snapseed_v5.py` | Run Snapseed v5 on the clustered object; write annotated h5ad, score tables, flags, UMAP, and a markdown report. |

Afterwards, we discovered a low quality cluster, which was removed 

---

## Scheduler

`scripts/run_annotation_pipeline.sh` runs the four steps in order, reading paths from `config/*.yaml`. Logs go under `logs/`.

```bash
cd annotation_pipeline   # or Organization/annotation_pipeline in this repo

# Full chain
bash scripts/run_annotation_pipeline.sh

# Preview commands without running
bash scripts/run_annotation_pipeline.sh --dry-run

# Resume from a step (skips earlier ones)
bash scripts/run_annotation_pipeline.sh --from concat
bash scripts/run_annotation_pipeline.sh --from cluster
bash scripts/run_annotation_pipeline.sh --from annotate

# Run a single step
bash scripts/run_annotation_pipeline.sh --only concat
bash scripts/run_annotation_pipeline.sh --only annotate
```

Environment overrides (optional):

- `PY` — Python binary (default: atlas conda env)
- `RSCRIPT` — Rscript binary
- `OMP_NUM_THREADS` / `MKL_NUM_THREADS` — thread caps (default 8)

Standalone annotate (no full chain):

```bash
bash scripts/run_annotate_snapseed_v5.sh
# or SLURM:
# sbatch scripts/sbatch_annotate_snapseed_v5.sh
```

---

## Config and markers

Edit these before running on a new machine:

| File | Purpose |
|------|---------|
| `config/rds_to_h5ad_parameters.yaml` | Input RDS directory, output h5ad directory |
| `config/concatenation_parameters.yaml` | Singlet input dir, concatenated output h5ad path |
| `config/snapseed_annotation_parameters.yaml` | Concat/clustered/annotated paths, Leiden resolution, score flag threshold (0.6) |
| `markers/snapseed_markers_v5.yaml` | Hierarchical marker definitions for Snapseed |

---

## Main outputs (annotate step)

Written under the `output_dir` from `snapseed_annotation_parameters.yaml` (default: `…/annotation_snapseed_v5/`):

| Output | Description |
|--------|-------------|
| `clustered_leiden_res1.0.h5ad` | Checkpoint after Leiden (cluster-only step) |
| `annotated_snapseed_v5_leiden_res1.0.h5ad` | Final object with Level_* labels and confidence flags |
| `snapseed_assignments.csv` | Cluster → Level_1 / Level_2 / Level_3 |
| `snapseed_cluster_annotation_scores.csv` | Per-cluster scores and flags |
| `snapseed_flagged_clusters_score_lt_0.6.csv` | Clusters with overall score &lt; 0.6 |
| `cells_per_Level_*.csv` | Cell counts per label |
| `umap_snapseed_annotation.png` | UMAP coloured by cluster / levels / flags |
| `annotation_report.md` | Short human-readable summary |

On the annotated object, useful `obs` columns include: `leiden_res_1.0`, `Level_1`, `Level_2`, `Level_3`, `Level_3_latest`, `snapseed_overall_score`, `snapseed_low_confidence`.

---

## Dependencies

- **R**: Seurat, SeuratDisk (RDS → h5ad)
- **Python**: scanpy, anndata, snapseed, pandas, numpy, matplotlib, PyYAML  
  Typical env used here: conda env `atlas` (`PY` / `RSCRIPT` as in the scheduler)

---

## Notes

- Leiden **res 1.0** here is the **pre-integration** clustering used for Snapseed. Later atlas work may use a higher resolution (e.g. Leiden res 2.0 after scPoli) for reannotation; that is a separate step, not this scheduler.
- Marker genes missing from the AnnData are skipped by Snapseed and listed in the annotate logs / `marker_genes_missing.csv` when present.
- Paths in `config/` currently point at lab absolute paths under `/data-master/...`; change them when cloning elsewhere.
