# HKOCA Integration Benchmark

A reproducible benchmarking pipeline for comparing single-cell integration methods on the HKOCA organoid atlas. The workflow subsamples the atlas into non-overlapping sample batches, runs twelve integration methods across repeated iterations, evaluates each method with [scIB](https://github.com/theislab/scib) metrics, and aggregates results into publication-ready summary tables and figures.

This package lives under `scripts_benchmark/push/` and is fully driven by a single configuration file (`benchmark_config.yaml`). All paths, hyperparameters, and method lists are defined there.

---

## Table of contents

1. [Motivation](#motivation)
2. [Overview](#overview)
3. [Directory layout](#directory-layout)
4. [Workflow](#workflow)
5. [Integration methods](#integration-methods)
6. [Script reference](#script-reference)
7. [Configuration](#configuration)
8. [Outputs](#outputs)
9. [Usage](#usage)
10. [Requirements](#requirements)
11. [Design decisions](#design-decisions)
12. [Troubleshooting](#troubleshooting)

---

## Motivation

Batch effects between organoid samples can dominate downstream analysis if integration is poor. Choosing an integration method for the HKOCA atlas requires a systematic comparison that balances:

- **Biological conservation**: cell types from the same lineage remain separable after integration.
- **Batch correction**: technical variation between samples is reduced without erasing biology.

This pipeline implements the [scIB benchmarking framework](https://github.com/theislab/scib) on repeated subsamples of the atlas. Each iteration uses a disjoint set of ten `sample_id` values, so performance estimates reflect variability across sample compositions rather than a single lucky split.

---

## Overview

```
benchmark_config.yaml
        |
        v
  workflow_scheduler  (one iteration per invocation)
        |
        +-- Step 1: python_integration_methods.py
        |           Subsample atlas, run 8 Python integrations, export MTX for R
        |
        +-- Step 2: R_integration_methods.R
        |           CCA/RPCA + CSS (Pearson/Spearman) embeddings
        |
        +-- Step 3: run_metrics_updated.py
                    Import R embeddings, compute scIB metrics per method

  (after all iterations complete)

  Visualisation_results.py
        Aggregate metrics, Wilcoxon tests, heatmaps, scorecard
```

Each invocation of `workflow_scheduler` runs all three steps sequentially for a single iteration. Submit multiple iterations by calling the scheduler repeatedly with different `--iter` values, or wrap it in your own cluster submission script.

---

## Directory layout

```
Benchmark/
├── benchmark_config.yaml          # Central configuration (paths, params)
├── workflow_scheduler             # Workflow entry point; orchestrates all three steps
├── python_integration_methods.py  # Step 1: Python integration methods
├── R_integration_methods.R        # Step 2: R integration methods
├── run_metrics_updated.py         # Step 3: scIB metric evaluation
├── Visualisation_results.py       # Post-hoc aggregation and plotting
├── README.md                      # This file
├── checkpoints/                   # Per-iteration outputs (created at runtime)
│   └── subsample_iter_N/
│       ├── checkpoint_<method>.h5ad
│       ├── heoca_counts.mtx
│       ├── heoca_metadata.csv
│       ├── heoca_genes.csv
│       ├── cca_embeddings.csv
│       ├── rpca_embeddings.csv
│       ├── css_pearson_embeddings.csv
│       ├── css_spearman_embeddings.csv
│       └── scib_raw_metrics_full.csv
└── logs/                          # Optional log directory (created at runtime)
    ├── bench_iter_1.log
    └── bench_iter_1.err
```

---

## Workflow

### Phase 0: Configuration

Edit `benchmark_config.yaml` before launching. At minimum, verify:

- `paths.master_atlas`: path to the unintegrated `.h5ad` input
- `paths.checkpoints_dir` and `paths.log_dir`: writable output directories
- `paths.conda_profile` and `paths.conda_env`: Python/R environment for GPU integrations
- `iteration.n_iterations` and `iteration.samples_per_iteration`: subsampling design

### Phase 1: Subsampling and Python integrations

`python_integration_methods.py` loads the master atlas, shuffles all `sample_id` values with a fixed seed, and selects a non-overlapping window of samples for the current array task. It recalculates batch-aware highly variable genes (HVGs), then runs eight Python methods and writes lightweight checkpoint files.

It also exports a Matrix Market (`.mtx`) file plus metadata and gene lists for the R step.

### Phase 2: R integrations

`R_integration_methods.R` reads the MTX export, runs CCA integration, RPCA integration, and Cluster Similarity Spectrum (CSS) with Pearson and Spearman correlation. Each method writes a cell-level embedding CSV to the iteration checkpoint directory.

### Phase 3: scIB metrics

`run_metrics_updated.py` converts R embedding CSVs into `.h5ad` checkpoints, then computes scIB metrics for all twelve methods. Metrics are written to `scib_raw_metrics_full.csv` per iteration.

### Phase 4: Visualization (manual, after all iterations)

`Visualisation_results.py` aggregates per-iteration CSVs, computes mean and standard deviation across iterations, runs paired Wilcoxon tests, and saves heatmaps and a master scorecard figure.

---

## Integration methods

| Method | Language | Description |
|--------|----------|-------------|
| `pca` | Python | Unintegrated PCA baseline |
| `harmony` | Python | Harmony batch correction on PCA |
| `scanorama` | Python | Scanorama mutual nearest-neighbor integration |
| `combat` | Python | ComBat linear batch correction |
| `bbknn` | Python | Batch-balanced kNN graph (PCA embedding for metrics) |
| `scvi` | Python | scVI variational autoencoder |
| `scanvi` | Python | scANVI semi-supervised VAE (warm-started from scVI) |
| `scpoli` | Python | scPoli conditional VAE |
| `cca` | R | CCA integration (Seurat) |
| `rpca` | R | RPCA integration (Seurat) |
| `css_pearson` | R | Cluster Similarity Spectrum, Pearson |
| `css_spearman` | R | Cluster Similarity Spectrum, Spearman |

---

## Script reference

### `workflow_scheduler`

Shell script that runs the full pipeline for one iteration. Parses `--config` and optional `--iter`, loads paths from the YAML file, activates the conda environment, and executes Steps 1 through 3 in order with validation gates between steps.

**Arguments**

| Flag | Required | Description |
|------|----------|-------------|
| `--config` | No (defaults to `./benchmark_config.yaml`) | Path to the YAML configuration file |
| `--iter` | No (defaults to `iteration.default_iter_id` in config) | Iteration index for subsample selection and checkpoint naming |

---

### `python_integration_methods.py`

Runs subsampling, HVG selection, and all Python integration methods for the current iteration (set via `--iter` or the environment variable defined in `iteration.slurm_env_var`).

**Key operations**

1. Load master atlas from `paths.master_atlas`
2. Select `samples_per_iteration` unique `sample_id` values (non-overlapping across iterations)
3. Recalculate HVGs (`subsample.n_top_genes`, batch-aware)
4. Export `heoca_counts.mtx`, `heoca_metadata.csv`, `heoca_genes.csv`
5. Run PCA, Harmony, Scanorama, ComBat, BBKNN, scVI, scANVI, scPoli
6. Save `checkpoint_<method>.h5ad` per method

**Checkpoint format**

Each checkpoint stores the embedding (`obsm["X_emb"]`), kNN graph (`obsp`), and observation metadata. The PCA baseline additionally retains raw counts in `.X` so scIB can compute the Principal Component Regression (PCR) batch metric.

**Arguments**

```
python python_integration_methods.py --config benchmark_config.yaml
```

---

### `R_integration_methods.R`

Runs Seurat and CSS integrations on the MTX export from Step 1.

**Key operations**

1. Load counts matrix and metadata from the iteration checkpoint directory
2. CCA: anchor finding and integrated PCA embedding
3. RPCA: anchor finding and integrated PCA embedding
4. CSS Pearson and CSS Spearman embeddings on a PCA-reduced master object

**Arguments**

```
Rscript R_integration_methods.R --config benchmark_config.yaml
```

---

### `run_metrics_updated.py`

Evaluates all twelve methods with scIB metrics for the current iteration.

**Key operations**

1. Import R embedding CSVs and write `checkpoint_<method>.h5ad` files
2. For each method: optimal Leiden resolution, then compute biological and batch metrics
3. Compute weighted overall score: `0.6 * bio_mean + 0.4 * batch_mean`
4. Write `scib_raw_metrics_full.csv`

**scIB metrics computed**

| Category | Metrics |
|----------|---------|
| Biological | ARI, NMI, ASW_celltype, iso_label_F1, iso_label_ASW, cLISI |
| Batch | ASW_batch, graph_connectivity, iLISI, kBET, PCR |
| Summary | `overall` (weighted composite) |

**Arguments**

```
python run_metrics_updated.py --config benchmark_config.yaml
```

---

### `Visualisation_results.py`

Aggregates results across all completed iterations. Run after all iterations have finished.

**Key operations**

1. Load `scib_raw_metrics_full.csv` from each `subsample_iter_N/`
2. Write `scib_mean_metrics.csv` and `scib_std_metrics.csv`
3. Wilcoxon signed-rank test: winner vs PCA baseline and winner vs runner-up
4. Save `scib_mean_sd_heatmap.png`, `scib_overall_errorbars.png`, `MASTER_SCORECARD.png`

**Arguments**

```
python Visualisation_results.py --config benchmark_config.yaml
```

---

## Configuration

All tunable parameters live in `benchmark_config.yaml`. Every script accepts `--config <path>`; if omitted, the default is `benchmark_config.yaml` in the same directory as the script.

### Main sections

| Section | Purpose |
|---------|---------|
| `paths` | Input atlas, output directories, conda environment |
| `iteration` | Number of iterations, samples per iteration, subsample seed, iteration env var |
| `metadata` | AnnData column names (`batch_key`, `label_key`, `counts_layer`) |
| `subsample` | HVG selection parameters |
| `integration` | PCA, Harmony, Scanorama, scVI, scANVI, scPoli hyperparameters |
| `r_integration` | Seurat dims, CSS label tag, embedding filenames |
| `export` | MTX export filenames |
| `metrics` | scIB method list, metric weights, LISI/kBET/PCR settings |
| `visualization` | Figure filenames, DPI, layout |

### Example: change iteration count

```yaml
iteration:
  n_iterations: 5
  samples_per_iteration: 10
```

And `visualization.n_iterations`:

```yaml
visualization:
  n_iterations: 5
```

### Example: point to a different atlas

```yaml
paths:
  master_atlas: /path/to/your/unintegrated_atlas.h5ad
```

---

## Outputs

### Per-iteration (`checkpoints/subsample_iter_N/`)

| File | Description |
|------|-------------|
| `checkpoint_<method>.h5ad` | Embedding + graph checkpoint per integration method |
| `heoca_counts.mtx` | Raw counts matrix for R (genes x cells) |
| `heoca_metadata.csv` | Cell metadata |
| `heoca_genes.csv` | Gene names |
| `cca_embeddings.csv` | CCA cell embeddings |
| `rpca_embeddings.csv` | RPCA cell embeddings |
| `css_pearson_embeddings.csv` | CSS Pearson embeddings |
| `css_spearman_embeddings.csv` | CSS Spearman embeddings |
| `scib_raw_metrics_full.csv` | scIB metrics for all methods, one iteration |

### Aggregated (`checkpoints/`)

| File | Description |
|------|-------------|
| `scib_mean_metrics.csv` | Mean metrics across iterations |
| `scib_std_metrics.csv` | Standard deviation across iterations |
| `scib_mean_sd_heatmap.png` | Side-by-side mean and SD heatmaps |
| `scib_overall_errorbars.png` | Overall score bar chart with error bars |
| `MASTER_SCORECARD.png` | Full scorecard: bio, batch, category, and overall scores |

### Logs (`logs/`)

| File | Description |
|------|-------------|
| `bench_iter_N.log` | Stdout for iteration N |
| `bench_iter_N.err` | Stderr for iteration N |

---

## Usage

### Full benchmark (recommended)

```bash
cd /path/to/scripts_benchmark/push

# 1. Edit configuration
vim benchmark_config.yaml

# 2. Run all iterations
for i in $(seq 1 10); do
  ./workflow_scheduler --config benchmark_config.yaml --iter "${i}"
done

# 3. Aggregate and plot
python Visualisation_results.py --config benchmark_config.yaml
```

### Run a single iteration

```bash
./workflow_scheduler --config benchmark_config.yaml --iter 3
```

---

## Requirements

### Python environment (`paths.conda_env`, default: `atlas`)

| Package | Used by |
|---------|---------|
| scanpy, anndata | All Python scripts |
| harmonypy | Harmony |
| scanorama | Scanorama |
| bbknn | BBKNN |
| scvi-tools, torch | scVI, scANVI |
| scarches | scPoli |
| scib | Metrics |
| pyyaml | Config loading |
| matplotlib, seaborn, scipy, pandas, numpy | Visualization |

### R packages

| Package | Used by |
|---------|---------|
| Seurat | CCA and RPCA integrations |
| Matrix | MTX reading |
| yaml | Config loading |
| simspec | CSS integration |

### Hardware

| Resource | Requirement |
|----------|-------------|
| GPU | Required for scVI, scANVI, scPoli (Step 1) |
| RAM | 128 GB recommended per iteration |
| CPUs | 12 recommended |

---

## Design decisions

**Non-overlapping subsampling.** Each iteration draws a distinct window of ten `sample_id` values from a globally shuffled list (seed 42). This requires at least `n_iterations * samples_per_iteration` unique samples in the atlas (default: 100).

**Lightweight checkpoints.** Only the PCA baseline retains a gene expression matrix. All other methods store embeddings and graphs only, keeping disk usage manageable across twelve methods and ten iterations.

**scIB as the evaluation standard.** scIB provides a balanced, published metric suite used widely in the single-cell integration literature. The weighted overall score (60% biological, 40% batch) follows the scIB paper defaults.

**Single config file.** All scripts read from `benchmark_config.yaml` via `--config`, avoiding hardcoded paths and making it straightforward to version-control parameter sets for different benchmark runs.

**No differential expression step.** DE concordance was removed from the pipeline. It required full expression matrices for every method and failed on lightweight checkpoints. The benchmark now focuses exclusively on embedding-based scIB metrics.

---

## Troubleshooting

**`STRATIFICATION ERROR: Your atlas only contains N unique samples`**

The atlas has fewer than `n_iterations * samples_per_iteration` unique `sample_id` values. Reduce `iteration.n_iterations` or `iteration.samples_per_iteration` in the config.

**Step 1 completes but `heoca_counts.mtx` is missing**

Check `logs/bench_iter_N.err` for Python errors during subsampling or MTX export. Verify `paths.master_atlas` exists and is readable.

**All R integrations failed (Step 2)**

Confirm the conda/R environment has Seurat and simspec installed. Check that Step 1 wrote all three MTX export files.

**scPoli fails but other methods succeed**

scPoli is wrapped in a try/except block and logs the error without aborting the pipeline. Inspect the log for the specific failure (memory, package version, etc.).

**`No scIB CSV files found` in visualization**

Not all iterations have completed. Run `Visualisation_results.py` only after every `subsample_iter_N/scib_raw_metrics_full.csv` exists.

**Config not found**

Pass an absolute path: `--config /full/path/to/benchmark_config.yaml`

---

## Citation

If you use this pipeline in a publication, cite the scIB framework:

> Luecken, M.D. et al. Benchmarking atlas-level data integration in single-cell genomics. *Nature Methods* (2022).

Method-specific citations apply for individual integration tools (Harmony, scVI, Seurat, and others).
