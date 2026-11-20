# Annotation

Workflow to turn QC-filtered single-cell objects into a Snapseed-annotated kidney organoid atlas:

**(optional) RDS → h5ad → concatenate → Leiden (res 1.0) → Snapseed v5**

This folder (`annotation/`) ships the runnable scripts, scheduler, configs, and marker hierarchy. Large intermediate `.h5ad` files are written to paths in `config/` (edit those for your machine).

---

## What is Snapseed?

Snapseed is a **marker-based hierarchical cell-type annotation** tool for single-cell data. Given:

1. an AnnData object with **cluster labels** (here: Leiden `leiden_res_1.0`), and  
2. a **YAML tree of cell types** with marker genes at each level,

it scores every cluster against each node in the hierarchy (default method: **AUROC** on marker expression) and assigns the best-matching labels.

This project uses **markers v5** (`markers/snapseed_markers_v5.yaml`), a kidney-organoid hierarchy with up to three levels:

| Level | Example |
|-------|---------|
| **Level_1** | Nephron, Stromal_cells, Endothelial_cells, Off_target_cells, … |
| **Level_2** | Podocytes, Proximal_tubule, S1/S2/S3, Neural_progenitor_cells, … |
| **Level_3** | Precursor_podocytes, Mature_podocytes, Dev_PT, PT, … |

`annotate_snapseed_v5.py` wraps `snapseed.annotate_hierarchy`, writes labels onto `adata.obs` (`Level_1`, `Level_2`, `Level_3`, `Level_3_latest`), and **flags clusters** whose overall annotation score is **&lt; 0.6** for manual review.

Clustering (normalize → HVG → PCA → neighbors → Leiden → UMAP) runs inside the same script when `leiden_res_1.0` is missing.

---

## Directory layout

```
annotation/
├── README.md
├── scripts/          # pipeline code + scheduler
├── config/           # paths and parameters (YAML)
└── markers/          # snapseed_markers_v5.yaml
```

---

## Pipeline steps

| # | Step | Script | Role |
|---|------|--------|------|
| 1 | **RDS → h5ad** *(optional)* | `scripts/rds_to_h5ad.r` via scheduler | Convert Seurat `.rds` → `.h5ad` (SeuratDisk). **Skipped** if filtered singlet `.h5ad` already exist (e.g. QC wrote them). Runs only when RDS are present and no matching h5ad are found. |
| 2 | **Concatenate** | `scripts/concatenate_filtered_singlets.py` | Merge filtered singlet h5ads into one atlas AnnData (shared gene space, standardized metadata). |
| 3 | **Cluster** | `scripts/annotate_snapseed_v5.py --cluster-only` | Leiden resolution **1.0**; write clustered checkpoint h5ad. |
| 4 | **Annotate** | `scripts/annotate_snapseed_v5.py` | Snapseed v5; write annotated h5ad, scores, flags, UMAP, report. |

### Optional conversion logic

The scheduler looks for `*_harmonized_singlets_filtered.h5ad` in:

- `config/concatenation_parameters.yaml` → `input_dir`, and  
- `config/rds_to_h5ad_parameters.yaml` → `output_h5ad_dir`

| Situation | Behaviour |
|-----------|-----------|
| Matching **h5ad** found | Skip conversion (QC already exported h5ad). Sync into concat `input_dir` if needed. |
| No h5ad, but **RDS** in `input_rds_dir` | Run RDS → h5ad, then sync into concat `input_dir`. |
| Neither | Error and stop. |

---

## Scheduler

`scripts/run_annotation.sh` runs the chain above, reading paths from `config/*.yaml`. Logs go under `logs/`.

```bash
cd annotation

# Full chain (conversion auto-skipped when h5ad already exist)
bash scripts/run_annotation.sh

# Preview commands without running
bash scripts/run_annotation.sh --dry-run

# Resume from a step
bash scripts/run_annotation.sh --from concat
bash scripts/run_annotation.sh --from cluster
bash scripts/run_annotation.sh --from annotate

# Run a single step
bash scripts/run_annotation.sh --only concat
bash scripts/run_annotation.sh --only annotate
```

Environment overrides (optional):

- `PY` — Python binary (default: atlas conda env)
- `RSCRIPT` — Rscript binary (needed only if conversion runs)
- `OMP_NUM_THREADS` / `MKL_NUM_THREADS` — thread caps (default 8)

Standalone annotate:

```bash
bash scripts/run_annotate_snapseed_v5.sh
# or: sbatch scripts/sbatch_annotate_snapseed_v5.sh
```

---

## Config and markers

| File | Purpose |
|------|---------|
| `config/rds_to_h5ad_parameters.yaml` | RDS / convert-output dirs (used only if conversion runs) |
| `config/concatenation_parameters.yaml` | Singlet h5ad input dir, concatenated output path |
| `config/snapseed_annotation_parameters.yaml` | Clustered/annotated paths, Leiden resolution, score flag threshold (0.6) |
| `markers/snapseed_markers_v5.yaml` | Hierarchical markers for Snapseed |

---

## Main outputs (annotate step)

Under `output_dir` from `snapseed_annotation_parameters.yaml`:

| Output | Description |
|--------|-------------|
| `clustered_leiden_res1.0.h5ad` | Checkpoint after Leiden |
| `annotated_snapseed_v5_leiden_res1.0.h5ad` | Final object with Level_* labels and flags |
| `snapseed_assignments.csv` | Cluster → Level_1 / Level_2 / Level_3 |
| `snapseed_cluster_annotation_scores.csv` | Per-cluster scores and flags |
| `snapseed_flagged_clusters_score_lt_0.6.csv` | Clusters with overall score &lt; 0.6 |
| `cells_per_Level_*.csv` | Cell counts per label |
| `umap_snapseed_annotation.png` | UMAP by cluster / levels / flags |
| `annotation_report.md` | Short summary |

Useful `obs` columns: `leiden_res_1.0`, `Level_1`, `Level_2`, `Level_3`, `Level_3_latest`, `snapseed_overall_score`, `snapseed_low_confidence`.

---

## Dependencies

- **R** (only if converting RDS): Seurat, SeuratDisk  
- **Python**: scanpy, anndata, snapseed, pandas, numpy, matplotlib, PyYAML  

---

## Notes

- Leiden **res 1.0** is the **pre-integration** clustering used for Snapseed. Later reannotation may use higher resolution (e.g. after scPoli); that is outside this scheduler.
- Missing marker genes are skipped by Snapseed and listed in annotate logs / `marker_genes_missing.csv` when written.
- Edit absolute paths in `config/` when cloning to another machine.
