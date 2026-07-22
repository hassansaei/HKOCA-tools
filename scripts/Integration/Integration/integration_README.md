# HKOCA Integration Tools

Documentation for the integration workflows in `scripts_benchmark/push/`. This directory contains two independent, config-driven pipelines:

1. **Integration benchmark** — compare twelve integration methods with scIB metrics across repeated subsamples.
2. **scPoli full-atlas integration** — train scPoli on the complete organoid atlas and produce an integrated reference object with QC figures.

Both pipelines read parameters from YAML config files and accept a `--config` flag. Original scripts outside `push/` are not modified.

---

## Table of contents

1. [Directory layout](#directory-layout)
2. [Pipeline 1: Integration benchmark](#pipeline-1-integration-benchmark)
3. [Pipeline 2: scPoli full-atlas integration](#pipeline-2-scpoli-full-atlas-integration)
4. [Configuration reference](#configuration-reference)
5. [Requirements](#requirements)
6. [Troubleshooting](#troubleshooting)

---

## Directory layout

```
push/
├── integration_README.md              # This file
├── README.md                          # Detailed benchmark documentation

# Benchmark pipeline
├── benchmark_config.yaml
├── workflow_scheduler
├── python_integration_methods.py
├── R_integration_methods.R
├── run_metrics_updated.py
├── Visualisation_results.py
├── checkpoints/                       # Created at runtime
└── logs/                              # Created at runtime

# scPoli full-atlas integration
├── scpoli_integration_config.yaml
└── run_integration_scpoli_parameters.py
```

---

## Pipeline 1: Integration benchmark

### Purpose

Systematically compare integration methods on the HKOCA organoid atlas using the [scIB framework](https://github.com/theislab/scib). Each iteration subsamples ten non-overlapping `sample_id` values, runs all methods, and records biological conservation and batch correction metrics.

### Workflow

```
benchmark_config.yaml
        |
        v
workflow_scheduler --config benchmark_config.yaml --iter N
        |
        +-- python_integration_methods.py     (8 Python methods + MTX export)
        +-- R_integration_methods.R           (CCA, RPCA, CSS Pearson, CSS Spearman)
        +-- run_metrics_updated.py            (scIB metrics for all 12 methods)

(after all iterations)

Visualisation_results.py                  (aggregate tables and figures)
```

### Methods evaluated

| Method | Language | Description |
|--------|----------|-------------|
| `pca` | Python | Unintegrated PCA baseline |
| `harmony` | Python | Harmony batch correction |
| `scanorama` | Python | Scanorama MNN integration |
| `combat` | Python | ComBat linear correction |
| `bbknn` | Python | Batch-balanced kNN |
| `scvi` | Python | scVI VAE |
| `scanvi` | Python | scANVI semi-supervised VAE |
| `scpoli` | Python | scPoli conditional VAE |
| `cca` | R | CCA integration (Seurat) |
| `rpca` | R | RPCA integration (Seurat) |
| `css_pearson` | R | Cluster Similarity Spectrum, Pearson |
| `css_spearman` | R | Cluster Similarity Spectrum, Spearman |

### Usage

```bash
cd scripts_benchmark/push

# Run one iteration
./workflow_scheduler --config benchmark_config.yaml --iter 1

# Run all iterations
for i in $(seq 1 10); do
  ./workflow_scheduler --config benchmark_config.yaml --iter "${i}"
done

# Aggregate results
python Visualisation_results.py --config benchmark_config.yaml
```

### Key outputs

| Location | File | Description |
|----------|------|-------------|
| `checkpoints/subsample_iter_N/` | `checkpoint_<method>.h5ad` | Per-method embedding checkpoints |
| `checkpoints/subsample_iter_N/` | `scib_raw_metrics_full.csv` | scIB metrics for one iteration |
| `checkpoints/` | `scib_mean_metrics.csv` | Mean metrics across iterations |
| `checkpoints/` | `MASTER_SCORECARD.png` | Summary scorecard figure |

See [README.md](README.md) for full benchmark documentation including config sections, metric definitions, and troubleshooting.

---

## Pipeline 2: scPoli full-atlas integration

### Purpose

Train scPoli on the **full** unintegrated organoid atlas (not a subsample) to produce a production-ready integrated reference. The script:

1. Loads and sanitizes the master atlas
2. Selects highly variable genes
3. Computes a PCA baseline for before/after comparison
4. Trains scPoli with study and sample_id as batch keys and Level_2 as the cell type key
5. Saves the integrated `.h5ad` object
6. Renders four before/after UMAP comparison figures

This pipeline is separate from the benchmark subsampling workflow. It is intended for generating the final integrated atlas object used downstream.

### Workflow

```
scpoli_integration_config.yaml
        |
        v
run_integration_scpoli_parameters.py --config scpoli_integration_config.yaml
        |
        +-- Load and sanitize master atlas
        +-- HVG selection (4,000 genes, batch_key=study)
        +-- Baseline PCA + UMAP
        +-- scPoli training
        +-- Integrated embedding + UMAP
        +-- Save .h5ad
        +-- Render comparison figures
```

### Script: `run_integration_scpoli_parameters.py`

The push copy is logic-identical to `parameters/run_integration_scpoli_parameters.py`. The only differences are:

- Hardcoded paths and parameters moved to `scpoli_integration_config.yaml`
- Emojis and unnecessary comments removed
- `--config` flag support added

All processing steps, model architecture, training parameters, and figure generation match the original script.

### Usage

```bash
cd scripts_benchmark/push

python run_integration_scpoli_parameters.py --config scpoli_integration_config.yaml
```

### Key outputs

| File | Description |
|------|-------------|
| `Master_Atlas_scPoli_Level2_Integrated_parameters.h5ad` | Integrated atlas with PCA and scPoli embeddings |
| `Comparison_Sample_Alignment.png` | Before/after UMAP colored by sample_id |
| `Comparison_Study_Alignment.png` | Before/after UMAP colored by study |
| `Comparison_Level2_Biology.png` | Before/after UMAP colored by Level_2 |
| `Comparison_Level1_Biology.png` | Before/after UMAP colored by Level_1 |

Default output directory: `parameters/` (configurable via `paths.output_dir`).

### Config: `scpoli_integration_config.yaml`

| Section | Purpose |
|---------|---------|
| `paths` | Input atlas, output directory, output filename |
| `metadata` | Observation columns to sanitize for missing values |
| `feature_selection` | HVG count, flavor, batch key |
| `pca` | PCA/UMAP parameters for baseline representation |
| `scpoli` | Model architecture (hidden layers, latent dim, condition keys) |
| `training` | Epochs, batch size, early stopping |
| `plotting` | Figure sizes, DPI, Level_1 color palette |

---

## Configuration reference

Both pipelines use the same `--config` pattern:

```bash
python <script>.py --config /path/to/config.yaml
```

If `--config` is omitted, each script defaults to its config file in the same directory.

### Benchmark config (`benchmark_config.yaml`)

Controls subsampling, all twelve integration methods, scIB metric weights, and visualization settings. Key sections: `paths`, `iteration`, `metadata`, `subsample`, `integration`, `r_integration`, `metrics`, `visualization`.

### scPoli config (`scpoli_integration_config.yaml`)

Controls the full-atlas scPoli training run. Key sections: `paths`, `metadata`, `feature_selection`, `pca`, `scpoli`, `training`, `plotting`.

### Example: change benchmark iteration count

```yaml
# benchmark_config.yaml
iteration:
  n_iterations: 5
  samples_per_iteration: 10

visualization:
  n_iterations: 5
```

### Example: change scPoli output location

```yaml
# scpoli_integration_config.yaml
paths:
  output_dir: /path/to/my/output
  output_h5ad: My_Integrated_Atlas.h5ad
```

---

## Requirements

### Python

| Package | Benchmark | scPoli |
|---------|-----------|--------|
| scanpy, anndata | Yes | Yes |
| scib | Yes | No |
| harmonypy, scanorama, bbknn | Yes | No |
| scvi-tools, scarches | Yes | Yes |
| pyyaml | Yes | Yes |
| matplotlib, seaborn | Yes | Yes |

### R

| Package | Benchmark | scPoli |
|---------|-----------|--------|
| Seurat | Yes | No |
| simspec | Yes | No |
| yaml | Yes | No |

### Hardware

| Pipeline | GPU | RAM |
|----------|-----|-----|
| Benchmark (per iteration) | Required for scVI/scANVI/scPoli | 128 GB recommended |
| scPoli full-atlas | Required | 300+ GB recommended |

---

## Troubleshooting

### Benchmark

**`STRATIFICATION ERROR`**: Atlas has fewer unique `sample_id` values than `n_iterations * samples_per_iteration`. Reduce iteration count or samples per iteration in `benchmark_config.yaml`.

**Missing MTX after Step 1**: Check logs for Python integration failures. Verify `paths.master_atlas` is readable.

**All R integrations failed**: Confirm Seurat and simspec are installed. Verify Step 1 produced `heoca_counts.mtx`.

### scPoli full-atlas

**anndata compatibility warnings at startup**: Expected on older anndata versions. The script includes compatibility patches for scarches/scvi-tools.

**Out of memory during training**: Reduce `feature_selection.n_top_genes` or `training.batch_size` in `scpoli_integration_config.yaml`.

**Missing cells after sanitization**: Cells with null values in `study`, `sample_id`, `Level_2`, or `Level_1` are removed. Check upstream annotation quality if cell loss is excessive.
