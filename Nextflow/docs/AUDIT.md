# Module and environment audit

Nextflow only orchestrates. Filtering thresholds, markers, integration parameters, annotation rules, and projection settings stay in `hkoca/` / `conda/`.

Package version at audit: `0.1.0-alpha.101` (`hkoca/_version.py`).

## Stage order (scientific)

| Nextflow stage | Existing CLI / entry | Notes |
|---|---|---|
| CellBender | `hkoca cellbender h5\|mtx` | Optional; skipped with `--skip-cellbender` / `run_cellbender=false` |
| Harmonization | `python -m hkoca.qc_filter.harmonize` (also step 1 of `hkoca qc-filter`) | Gene-space harmonize + Seurat RDS (`--to-rds`) |
| QC | `hkoca.qc_filter.qc_runner` / `QC_scdbl_Combined.R` (step 2 of `hkoca qc-filter`) | scDblFinder + MAD/QC floors + optional H5AD export |
| Annotation | `hkoca annotation run` | Snapseed + Leiden (0.4, 0.6, 1.0) |
| Integration | `hkoca integration prep` then `hkoca integration run` | SCT prep, then Harmony / RPCA / CCA (+ optional scIB) |
| Projection | `hkoca projection map` | scPoli surgery onto frozen atlas UMAP (GPU) |

Original `hkoca pipeline` stages are `cellbender -> qc_filter -> annotation -> integration -> [projection]`, where `qc_filter` always runs harmonize then QC. Nextflow splits that into two processes for clearer dependency tracking, still calling the same modules.

## Environments (`conda/`)

| Env name (`name:`) | YAML | Used by |
|---|---|---|
| `hkoca_cellbender` | `environment_cellbender.yaml` | CellBender (PyTorch cu117 + cellbender 0.3.2) |
| `hkoca_harmonize` | `environment_harmonize.yaml` | Harmonization, annotation, pipeline driver; Snapseed from GitHub |
| `sc_qc_pipeline` | `environment_qc.yaml` | QC Rscript (scDblFinder, Seurat, anndata) |
| `hkoca_integration` | `environment_integration.yaml` | Integration Seurat R (+ reticulate) |
| `hkoca_projection` | `environment_projection.yaml` | scPoli / scArches + PyTorch cu121 |
| `hkoca_annotation` | `environment_annotation.yaml` | Optional annotation-only; prefer Snapseed inside `hkoca_harmonize` |

Cross-env wiring (unchanged): `hkoca/conda_env.py` resolves sibling envs; QC uses `HKOCA_QC_ENV` / `HKOCA_RSCRIPT`; projection uses `HKOCA_ROOT` + `hkoca_projection` Python.

## Per-stage findings

### CellBender

- **Input:** CellRanger `raw_feature_bc_matrix.h5` (h5 mode) or MTX directory (mtx mode) under each sample directory.
- **Command:** `cellbender remove-background --input ... --output {sample_id}_filtered.h5` with packaged defaults (epochs=100, droplets=35000, lr=5e-5, threads=24, cuda=true).
- **Config:** `hkoca/cellbender/cellbender.config`.
- **Output:** `{sample_dir}/{sample_id}_filtered.h5`; CellBender also emits a cell-called `{sample_id}_filtered_filtered.h5`, which the pipeline uses as harmonize `data_path` after CellBender.
- **Env:** `hkoca_cellbender`.

### Harmonization

- **Input:** Metadata CSV (`data_path`, `sample_id`, `study`, ...), GRCh38 GTF, optional transgenes.
- **Command:** Harmonize CLI with `--to-rds` (forced by `qc-filter`).
- **Config:** `hkoca/qc_filter/harmonize/harmonize.config`.
- **Output:** `{output}/rds/{study}_harmonized.rds` (and intermediate `raw/` / `harmonized/` h5ad).
- **Env:** `hkoca_harmonize`.

### QC

- **Input:** Harmonized RDS under harmonize `output_root` matching `_harmonized\.rds$`.
- **Command:** `QC_scdbl_Combined.R --run all` via `sc_qc_pipeline` Rscript.
- **Config:** `hkoca/qc_filter/r/qc_config.dcf` (floors: `qc_min_features=200`, scDblFinder batch=`sample_id`, H5AD conversion on).
- **Output:** `qc_filtered_rds/{study}*_filtered.rds`, `h5ad_converted/{study}*_filtered.h5ad`, `doublet_filtered_rds/`, `Summary/`.
- **Env:** `sc_qc_pipeline` (R); invoked from harmonize env via env resolution.

### Annotation

- **Input:** QC-filtered `.h5ad`; marker YAML (default `hkoca/config/snapseed_markers_v4.yaml`).
- **Command:** `hkoca annotation run --input ... --output-dir ... --markers ...`.
- **Config:** `hkoca/config/annotation.config.yaml` (seed=42, HVG=3000, PCA=50, resolutions 0.4/0.6/1.0).
- **Output:** `annotated_obj/*_annotated.h5ad`, `clustered_obj/`, `plots/`.
- **Env:** `hkoca_harmonize` (+ Snapseed).

### Integration

- **Input:** QC-filtered Seurat RDS + optional annotated h5ad for label transfer.
- **Command:** `hkoca integration prep` then `hkoca integration run --methods harmony,rpca,cca`.
- **Config:** `hkoca/integration/r/integration.config.dcf` (SCT n_features=2500, silhouette seed=42, min resolution 0.4, `R_MAX_VSIZE=150Gb`).
- **Output:** `prep/sct_prepared.rds`, `objects/integrated_{harmony,rpca,cca}.rds`, plots/tables/benchmark.
- **Env:** `hkoca_integration` (Rscript); annotation Python via reticulate/`HKOCA_ANNOTATION_PYTHON`.

### Projection

- **Input:** Integration prep `sct_prepared.rds` (RNA counts), atlas h5ad (`Level_*_Integrated`, `X_umap_scpoli`), scPoli reference model dir.
- **Command:** `hkoca projection map --query ... --atlas ... --model-dir ... --output-dir ...`.
- **Config:** `hkoca/projection/projection.config.yaml` (seed=42, n_epochs=50, pretrain=40, eta=10).
- **Output:** `projected_obj/sct_prepared_projected.h5ad`, plots, tables, surgery model.
- **Env:** `hkoca_projection` (GPU); `HKOCA_ROOT` must point at this clone root.

## Non-goals

- Do not change QC thresholds, marker hierarchies, SCT/integration settings, or scPoli hyperparameters.
- Do not overwrite analysis modules under `hkoca/`.
- Do not mutate unrelated named conda envs; Nextflow should create caches from the YAML files.
