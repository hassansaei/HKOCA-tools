# File dependency map

Each stage is declared so Nextflow can track inputs and outputs. Commands call existing `hkoca` modules; parameters stay in packaged configs unless overridden by CLI/params.

## 1. CellBender

| Field | Spec |
|---|---|
| **Inputs** | Per-sample CellRanger raw H5 or MTX dir; optional `cellbender.config`; sample_info columns `sample_id`, `sample_dir` / `data_path`, `run_cellbender` |
| **Command** | `hkoca cellbender h5` or `hkoca cellbender mtx` (same flags as today) |
| **Environment** | `conda/environment_cellbender.yaml` (`hkoca_cellbender`) |
| **Outputs** | `{sample_id}_filtered.h5` (+ CellBender cell-called `*_filtered_filtered.h5` used as harmonize `data_path`) |
| **Resources** | GPU preferred (`cuda`); default 24 CPU threads; high GPU memory for large droplet counts |
| **Validation** | Non-empty cell-called H5; sample_id present in resume set |

## 2. Harmonization

| Field | Spec |
|---|---|
| **Inputs** | Harmonize metadata CSV (`data_path` must point at matrices after CellBender rewrite); GTF; optional transgenes; `harmonize.config` |
| **Command** | `python -m hkoca.qc_filter.harmonize --csv ... --gtf ... --output ... --to-rds` |
| **Environment** | `conda/environment_harmonize.yaml` (`hkoca_harmonize`) |
| **Outputs** | `{output}/rds/{study}_harmonized.rds` (plus intermediate h5ad under `harmonized/` / `raw/`) |
| **Resources** | Multi-CPU; memory scales with genes x cells across samples in a study |
| **Validation** | Non-empty `_harmonized.rds` per study in sample CSV |

## 3. QC

| Field | Spec |
|---|---|
| **Inputs** | Harmonize output root (recursive discovery of `*_harmonized.rds`); `qc_config.dcf`; `--run all\|doublet\|qc\|h5ad` |
| **Command** | Rscript `QC_scdbl_Combined.R` via `hkoca.qc_filter.qc_runner` (or equivalent argv) |
| **Environment** | `conda/environment_qc.yaml` (`sc_qc_pipeline`) |
| **Outputs** | `qc_filtered_rds/*_filtered.rds`, `h5ad_converted/*_filtered.h5ad`, `doublet_filtered_rds/`, `Summary/`, logs |
| **Resources** | Multi-CPU; R heap grows with object size |
| **Validation** | Non-empty filtered RDS + H5AD per study (when H5AD conversion enabled) |

## 4. Annotation

| Field | Spec |
|---|---|
| **Inputs** | QC-filtered H5AD; marker YAML (`snapseed_markers_v4.yaml`); optional `annotation.config.yaml` |
| **Command** | `hkoca annotation run --input <filtered.h5ad> --output-dir ... --markers ...` |
| **Environment** | `conda/environment_harmonize.yaml` (Snapseed) |
| **Outputs** | `annotated_obj/{stem}_annotated.h5ad`, `clustered_obj/`, `plots/` |
| **Resources** | Multi-CPU; memory for HVG/PCA/Leiden + Snapseed |
| **Validation** | Non-empty `*_annotated.h5ad` with Level_* columns for configured resolutions |

## 5. Integration

| Field | Spec |
|---|---|
| **Inputs** | QC-filtered RDS; annotated H5AD (label transfer); `integration.config.dcf`; methods list |
| **Command** | `hkoca integration prep` then `hkoca integration run` |
| **Environment** | `conda/environment_integration.yaml` (`hkoca_integration`); reticulate may use harmonize Python |
| **Outputs** | `prep/sct_prepared.rds`, `objects/integrated_{method}.rds`, method UMAPs, tables, optional `benchmark/` |
| **Resources** | High memory (`R_MAX_VSIZE=150Gb` default); multi-CPU |
| **Validation** | Non-empty `sct_prepared.rds` and each requested `integrated_*.rds` |

## 6. Projection

| Field | Spec |
|---|---|
| **Inputs** | `prep/sct_prepared.rds`; atlas H5AD; scPoli `model-dir` (`model_params.pt`, `attr.pkl`, `var_names.csv`); `projection.config.yaml` |
| **Command** | `hkoca projection map --query ... --atlas ... --model-dir ... --output-dir ...` |
| **Environment** | `conda/environment_projection.yaml` (`hkoca_projection`); `HKOCA_ROOT=<clone root>` |
| **Outputs** | `projected_obj/sct_prepared_projected.h5ad`, `plots/`, `tables/`, `models/scpoli_query_surgery/` |
| **Resources** | GPU node; VRAM depends on query size and surgery epochs |
| **Validation** | Non-empty projected H5AD with `Level_3_pred` / `X_scpoli` (or equivalent keys from runner) |

## Real dependency graph (not simplified)

```
sample_info.csv + raw matrices
        |
        v
   [CellBender] ----optional skip----\
        |                             |
        v                             v
  rewritten meta CSV -----------> [Harmonization] --GTF-->
                                        |
                                        v
                                      [QC]
                                   /         \
                                  v           v
                         filtered.rds    filtered.h5ad
                                  \           /
                                   \         v
                                    \   [Annotation]
                                     \       |
                                      \      v
                                       [Integration prep+run]
                                              |
                                              |  + atlas.h5ad + model-dir
                                              v
                                         [Projection]
```

Projection does **not** consume integrated Harmony/RPCA/CCA objects by default; it consumes the SCT prep RDS (`sct_prepared.rds`) plus external atlas and model references.
