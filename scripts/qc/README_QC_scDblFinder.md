# Single-Cell QC Pipeline

**`qc_pipeline.R`** — CSV-driven quality control filtering followed by batch-aware doublet detection for single-cell RNA-seq data. Both stages are controlled by a single configuration file and share a unified log.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Repository layout](#3-repository-layout)
4. [Quick start](#4-quick-start)
5. [Configuration file reference](#5-configuration-file-reference)
6. [CLI reference](#6-cli-reference)
7. [Stage 1 — QC Filtering](#7-stage-1--qc-filtering)
8. [Stage 2 — Doublet Detection](#8-stage-2--doublet-detection)
9. [Outputs](#9-outputs)
10. [Threshold method reference](#10-threshold-method-reference)
11. [Parameter-Aware Resume Logic](#11-parameter-aware-resume-logic)
12. [Common workflows](#12-common-workflows)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. Overview

The pipeline runs in two sequential stages:

```
Raw harmonized RDS files
      │
      ▼  Stage 1 — QC Filtering
      │  • Per-dataset CSV-driven thresholds (nFeature, nCount, %mito)
      │  • 16 QC PNGs per dataset (violin + scatter, raw + filtered)
      │  • Before/After audit PDF
      │  • Full 7-panel dashboard (ridge plots, inflation ratios, ...)
      │  • Filtered RDS files  ──────────────────────────────────────┐
      │                                                              │
      ▼  Stage 2 — Doublet Detection                                │
      │  • Reads filtered RDS files ◄──────────────────────────────┘
      │  • Ghost-cell filter (removes CellBender low-count artifacts)
      │  • Normalisation → PCA → UMAP
      │  • scDblFinder (batch-aware when multi-sample metadata present)
      │  • Per-sample configurable expected doublet rate (DBR)
      │  • Singlet-only RDS + full-call RDS
      │  • 4-panel audit PDF per sample
      └  • Doublet summary CSV
```

Both stages read from the **same** `qc_config.dcf` file and write to the **same** log file.

---

## 2. Requirements

### Package Management & Runtime Environment

This pipeline uses Conda to cleanly manage dependencies (R ≥ 4.2, Seurat, scDblFinder, and critical system libraries for robust PDF/PNG rendering like Cairo and HarfBuzz). 

The script strictly prevents in-script installation to maintain cluster security and absolute reproducibility.

### Step 1 — Create and activate the Conda Environment

An `environment.yml` file is provided. You must build and activate this environment before running the R script.

```bash
# Create the environment (only needed the first time)
conda env create -f environment.yml

# Activate the environment
conda activate hkoca_qc_env
```

Full package list managed by Conda:

| Package | Source | Role |
|---|---|---|
| `r-base` | conda-forge | R runtime |
| `cairo`, `harfbuzz`, `fribidi` | conda-forge | System dependencies for plotting |
| `r-seurat` | conda-forge | Seurat object handling |
| `bioconductor-singlecellexperiment` | bioconda | SCE conversion for scDblFinder |
| `bioconductor-scdblfinder` | bioconda | Doublet detection |
| `r-ggplot2` | conda-forge | All plotting |
| `r-patchwork` | conda-forge | Multi-panel plot composition |
| `r-cowplot` | conda-forge | Plot layout helpers |
| `r-ggridges` | conda-forge | Ridge plots in dashboard |
| `r-reshape2` | conda-forge | Data reshaping for dashboard |
| `r-mass` | conda-forge | 2-D kernel density estimation |
| `r-scales` | conda-forge | Axis/colour helpers |
| `r-yaml` | conda-forge | Config parsing |
| `r-dplyr` | conda-forge | Data manipulation |
| `r-ggrastr` | conda-forge | Raster rendering for >50k-cell scatters |

---

## 3. Repository layout

```
scripts/qc/
├── qc_scDblFinder_combined.R          # Main pipeline script
├── README_QC_scDblFinder.md           # this file
├── qc_config.dcf          # Configuration file (edit this)
├── QC_Decisions.csv
│     
├── environment.yml 
└── qc_results/            # Output root (set via output_dir in config)
    ├── logs/
    │   └── qc_pipeline.log
    ├── Summary/           # QC audit outputs
    │   ├── QC_Before_After_Report.pdf
    │   ├── summary_qc_full_dashboard.png
    │   ├── qc_summary_detailed.csv
    │   └── <dataset_name>/
    │       ├── Violin_qc.raw.nFeature_RNA.png
    │       ├── Plot_qc.raw.nFeature_RNA.png
    │       └── ...        # 16 PNGs per dataset
    ├── qc_filtered_rds/   # Stage 1 output (input for Stage 2)
    │   └── <name>_filtered.rds
    └── scDblFinder_results/   # Stage 2 output
        ├── Doublet_Summary/
        │   ├── Doublet_Audit_Report_<timestamp>.pdf
        │   └── doublet_summary_<timestamp>.csv
        ├── Doublet_filtered_rds/
        │── <name>_singlets.rds
        └── <name>_with_doublet_calls.rds
```

---

## 4. Quick start

### Step 1 — Initialize the Environment

```bash
conda env create -f environment.yml
conda activate hkoca_qc_env
```

### Step 2 — Edit the config file

Open `qc_config.dcf` and set at minimum:

```
rds_dir:        /path/to/your/raw_rds
decisions_csv:  /path/to/your/QC_Decisions.csv
output_dir:     /path/to/your/output
```

### Step 3 — Run the full pipeline

```bash
Rscript qc_pipeline.R --config qc_config.dcf
```

This runs both stages in sequence: QC filtering first, then doublet detection on the filtered output.

---

## 5. Configuration file reference

The config file uses **DCF format** (Debian Control File): one `key: value` pair per line, no quotes, `#` for comments.

```
key: value
```

All keys are described below, grouped by section.

---

### Section 1 — Input paths

| Key | Required | Description |
|---|---|---|
| `rds_dir` | Yes | Directory containing the raw Seurat RDS files to process |
| `decisions_csv` | Yes | Path to the CSV file with per-dataset QC threshold methods (see [Section 7](#7-stage-1--qc-filtering) for required columns) |

---

### Section 2 — Output paths

| Key | Default | Description |
|---|---|---|
| `output_dir` | `postqc7` | Root directory for all pipeline outputs |
| `summary_subdir` | `Summary` | Subdirectory under `output_dir` for QC plots and CSVs |
| `filtered_subdir` | `qc_filtered_rds` | Subdirectory for QC-filtered RDS files (Stage 1 → Stage 2 handoff) |
| `doublet_subdir` | `doublet_filtered_rds` | Subdirectory for doublet-detection outputs |
| `log_file` | `logs/qc_pipeline.log` | Path to the shared log file (relative paths are resolved relative to `output_dir`) |

---

### Section 3 — File discovery

| Key | Default | Description |
|---|---|---|
| `rds_pattern` | `\.rds$` | Regex pattern used to find RDS files in `rds_dir` |
| `recursive_discovery` | `FALSE` | Set `TRUE` to search `rds_dir` recursively |

---

### Section 4 — Dataset exclusions

| Key | Default | Description |
|---|---|---|
| `skip_noisy_datasets` | `TRUE` | Whether to skip datasets listed in `noisy_datasets` |
| `noisy_datasets` | *(hardcoded list)* | Comma-separated list of RDS basenames to exclude. If left blank, the pipeline uses its internal hardcoded exclusion list |

---

### Section 5 — scDblFinder parameters

| Key | Default | Description |
|---|---|---|
| `dbl_batch_col` | `sample_id` | Metadata column that identifies individual samples within a merged object. When multiple unique values are present, scDblFinder runs independently per sample (batch-aware mode) |
| `dbl_min_count` | `100` | Minimum UMI count to keep a cell before running scDblFinder. Cells below this threshold are removed to prevent size-factor estimation failures |
| `dbl_min_feature` | `50` | Minimum number of detected genes. Same purpose as `dbl_min_count` |
| `dbl_umap_dims` | `20` | Number of PCA dimensions used when computing the UMAP (used only for visualisation in the audit PDF) |
| `dbl_default_platform` | `10x` | Default sequencing platform applied to any sample not listed in Section 6. For `10x`, scDblFinder auto-estimates the doublet rate from cell count when no explicit DBR is given |

---

### Section 6 — Per-sample platform and doublet rate (DBR) overrides

Per-sample settings are encoded directly in the config file using a dotted-key convention:

```
sample.<SAMPLE_NAME>.platform:  <platform>
sample.<SAMPLE_NAME>.dbr:       <numeric>
```

The `<SAMPLE_NAME>` must exactly match the RDS filename after stripping the suffix:

| RDS filename | Sample name used in config |
|---|---|
| `my_sample_filtered.rds` | `my_sample` |
| `my_sample_harmonized_filtered.rds` | `my_sample` |
| `my_sample.rds` | `my_sample` |

**Platform rules:**

| Platform | DBR key | Behaviour |
|---|---|---|
| `10x` | Optional | If omitted, scDblFinder auto-estimates from cell count |
| `10x` | Provided | Uses your value directly (must be in `(0, 0.30]`) |
| Any other (e.g. `dropseq`, `indrops`) | **Required** | Pipeline stops with an error if DBR is missing |

**Example config entries:**

```
# 10x sample with a high cell load — explicit DBR recommended
sample.kidney_organoid_rep1.platform: 10x
sample.kidney_organoid_rep1.dbr: 0.08

# 10x sample — auto-estimate (most common case, no dbr needed)
sample.kidney_organoid_rep2.platform: 10x

# Drop-seq sample — DBR required
sample.dropseq_run_A.platform: dropseq
sample.dropseq_run_A.dbr: 0.05
```

Samples not listed in this section inherit `dbl_default_platform` and use auto-estimation.

---

## 6. CLI reference

```
Rscript qc_pipeline.R [--config PATH] [--stage STAGE] [OPTIONS]
```

### Flags

| Flag | Default | Description |
|---|---|---|
| `--config PATH` | `qc_config.dcf` (next to script) | Path to the configuration file |
| `--stage STR` | `all` | Which stage(s) to run: `all`, `qc`, or `doublet` |
| `--rds_dir PATH` | *(from config)* | Override raw RDS input directory |
| `--decisions_csv PATH` | *(from config)* | Override QC decisions CSV |
| `--output_dir PATH` | *(from config)* | Override output base directory |
| `--summary_subdir NAME` | *(from config)* | Override QC summary subdirectory name |
| `--filtered_subdir NAME` | *(from config)* | Override filtered RDS subdirectory name |
| `--doublet_subdir NAME` | *(from config)* | Override doublet output subdirectory name |
| `--log_file PATH` | *(from config)* | Override log file path |
| `--recursive_discovery` | *(from config)* | Enable recursive RDS file search |
| `--rds_pattern REGEX` | *(from config)* | Override RDS filename regex |
| `--force_overwrite` | `FALSE` | Force reprocessing if outputs already exist |
| `--help` | — | Print usage and exit |

**CLI flags always override the config file.** This allows you to use a shared base config and customise individual runs without editing the file.

### Examples

```bash
# Full pipeline with default config
Rscript qc_pipeline.R --config qc_config.dcf

# QC only (skip doublet detection)
Rscript qc_pipeline.R --config qc_config.dcf --stage qc

# Doublet detection only (on already-filtered RDS)
Rscript qc_pipeline.R --config qc_config.dcf --stage doublet

# Override output directory at runtime
Rscript qc_pipeline.R --config qc_config.dcf --output_dir /scratch/my_run

# Point to a different RDS directory without editing the config
Rscript qc_pipeline.R --config qc_config.dcf --rds_dir /data/new_batch

# Run in background and redirect output
nohup Rscript qc_pipeline.R --config qc_config.dcf > run.log 2>&1 &
```

---

## 7. Stage 1 — QC Filtering

### Input

- Raw Seurat RDS files from `rds_dir`
- `QC_Decisions.csv` from `decisions_csv`

### QC_Decisions.csv format

The CSV must contain one row per dataset and these columns:

| Column | Description |
|---|---|
| `Dataset_Name` | RDS filename (basename, with or without `.rds`) |
| `Lower_Feature_Method` | Method for the lower nFeature_RNA threshold |
| `Upper_Feature_Method` | Method for the upper nFeature_RNA threshold |
| `Lower_Count_Method` | Method for the lower nCount_RNA threshold |
| `Upper_Count_Method` | Method for the upper nCount_RNA threshold |
| `Upper_Mito_Method` | Method for the upper percent.mito threshold |

Method values accepted in each column — see [Section 10](#10-threshold-method-reference) for full details.

### What Stage 1 does, per dataset

1. **Loads** the Seurat object and slims it (`DietSeurat`)
2. **Computes** QC metrics if absent: `nCount_RNA`, `nFeature_RNA`, `percent.mito`, `percent.ribo`
3. **Removes** barcodes with zero counts (empty droplets)
4. **Computes thresholds** by applying the algorithm named in each CSV column to the metric distribution
5. **Filters** cells outside the threshold windows
6. **Generates 16 PNGs** per dataset (4 features × 2 plot types × raw/filtered):
   - `Violin_qc.raw.<feature>.png` — violin with cutoff lines
   - `Plot_qc.raw.<feature>.png` — scatter coloured by density with cutoff lines
   - `Violin_qc.filtered.<feature>.png` — violin after filtering (no cutoff lines)
   - `Plot_qc.filtered.<feature>.png` — scatter after filtering
7. **Appends two pages** to the audit PDF (before and after density plots)
8. **Saves** the filtered object as `<name>_filtered.rds` in `filtered_subdir`
9. **Logs** an inflation ratio diagnostic — if the median nCount increases by more than 1.5× after filtering, a warning is emitted (this indicates over-aggressive low-end cutoffs)

### Dashboard

After all datasets are processed, a 7-panel `summary_qc_full_dashboard.png` is produced:

| Panel | Content |
|---|---|
| P1 | Cell retention bar chart (before vs after) |
| P2 | Filtering intensity (% cells removed) with 35 % warning line |
| P3 | Median nCount inflation ratio with thresholds at 1.0×, 1.4×, 1.5× |
| P4 | Median nCount_RNA before vs after (line plot) |
| P5 | Median nFeature_RNA before vs after (line plot) |
| P6 | Median %mito and %ribo before vs after (facet bar) |
| P7a | Per-dataset log1p(nCount) density with cutoff vlines |
| P7b | Per-dataset log1p(nFeature) density with cutoff vlines |
| P8a | Ridge plot: log1p(nCount) across all datasets |
| P8b | Ridge plot: log1p(nFeature) across all datasets |

---

## 8. Stage 2 — Doublet Detection

### Input

- QC-filtered RDS files from `filtered_subdir` (produced by Stage 1, or pre-existing if running `--stage doublet` directly)

### What Stage 2 does, per sample

1. **Ghost-cell filter** — removes cells with `nCount_RNA ≤ dbl_min_count` or `nFeature_RNA ≤ dbl_min_feature`. These residual low-quality cells (common after CellBender) cause size-factor estimation to fail inside scDblFinder.
2. **S4SXP fix** — coerces `meta.data` to a base `data.frame`, resolving a Seurat v5 S4 slot incompatibility with `as.SingleCellExperiment()`.
3. **Normalisation** — `NormalizeData → FindVariableFeatures → ScaleData → RunPCA → RunUMAP` (UMAP is used only for the audit plot).
4. **Doublet rate resolution** — looks up the sample in `DBL_SAMPLE_CFG` (populated from `qc_config.dcf`, Section 6). If no override is found, falls back to `dbl_default_platform` and auto-estimation.
5. **scDblFinder** — runs in either:
   - **Batch mode**: when `dbl_batch_col` is present in metadata with more than one unique value, each sub-sample is processed independently. This prevents the spurious ~30 % doublet rate that occurs when a merged multi-sample object is treated as a single pool.
   - **Single mode**: standard run on the whole object.
6. **Result mapping** — `scDblFinder.class` (singlet/doublet) and `scDblFinder.score` are added back to the Seurat metadata.
7. **Audit plots** — a 4-panel figure is appended to the PDF per sample: UMAP coloured by class, complexity scatter, feature violin, score density.
8. **Saves** two RDS files:
   - `<name>_with_doublet_calls.rds` — full object with doublet annotations
   - `<name>_singlets.rds` — singlets only (ready for downstream analysis)

---

## 9. Outputs

### Stage 1 outputs

| Path | Description |
|---|---|
| `<output_dir>/<summary_subdir>/QC_Before_After_Report.pdf` | Density plots before and after filtering for every dataset |
| `<output_dir>/<summary_subdir>/summary_qc_full_dashboard.png` | 7-panel composite dashboard |
| `<output_dir>/<summary_subdir>/qc_summary_detailed.csv` | Per-dataset stats table (40+ columns: cell counts, medians, cutoffs, inflation ratio) |
| `<output_dir>/<summary_subdir>/<dataset>/Violin_qc.*.png` | 8 violin PNGs per dataset |
| `<output_dir>/<summary_subdir>/<dataset>/Plot_qc.*.png` | 8 scatter PNGs per dataset |
| `<output_dir>/<filtered_subdir>/<name>_filtered.rds` | QC-filtered Seurat object (input for Stage 2) |

### Stage 2 outputs

| Path | Description |
|---|---|
| `<output_dir>/<doublet_subdir>/Summary/Doublet_Audit_Report_<ts>.pdf` | 4-panel audit plots per sample |
| `<output_dir>/<doublet_subdir>/Summary/doublet_summary_<ts>.csv` | Per-sample doublet counts and rates |
| `<output_dir>/<doublet_subdir>/<name>_with_doublet_calls.rds` | Full object with `scDblFinder.class` and `scDblFinder.score` in metadata |
| `<output_dir>/<doublet_subdir>/<name>_singlets.rds` | Singlets-only Seurat object |

### Shared

| Path | Description |
|---|---|
| `<output_dir>/logs/qc_pipeline.log` | Timestamped log covering both stages |

---

## 10. Threshold method reference

These strings are used in `QC_Decisions.csv` column values:

| Method string | Algorithm | Notes |
|---|---|---|
| `lower_tri` | Triangle threshold (lower tail) | Finds the bin at maximum perpendicular distance from the line connecting the histogram peak to the left endpoint |
| `upper_tri` | Triangle threshold (upper tail) | Same, but peak to right endpoint |
| `renyi` | Rényi entropy binarisation | Maximises the sum of Rényi entropies of the two histogram partitions |
| `knee` | Kneedle algorithm | Fits a curve to the ranked values and finds the point of maximum curvature |
| `mad3` | Median ± 3 × MAD | For `upper_*` columns: median + 3×MAD; for `lower_*`: median − 3×MAD |
| `none` | No threshold | Sets the cutoff to `Inf` (upper) or `-Inf` (lower) — effectively no filter |
| `manual_<value>` | Fixed numeric | Sets the cutoff to `<value>` directly, e.g. `manual_200` |

---

## 11. Parameter-Aware Resume Logic

The pipeline uses a **Smart Metadata Checkpoint** system to enable seamless and granular re-runs without requiring the user to manually track parameters or frequently use the `--force_overwrite` flag.

### How it works:
1. When a dataset finishes processing in Stage 1, the pipeline saves the filtered RDS object alongside a tiny, hidden `.json` sidecar file (e.g., `<output_dir>/qc_filtered_rds/.<dataset>_filtered.rds.meta.json`).
2. This JSON file contains a cryptographic hash of the exact threshold logics (from `QC_Decisions.csv`) used to produce that dataset.
3. **If you restart the pipeline:** it cross-references the current CSV configurations with the `.json` sidecars.
   - **Match**: If the dataset was already processed and the configuration has not changed, it is gracefully skipped (saving time while still injecting its metrics into the global dashboard).
   - **Mismatch**: If you tweaked even a single threshold for that sample, the pipeline detects the change and automatically re-runs *only* that dataset. The rest are skipped.

This ensures your workflow remains completely safe, idempotent, and hands-free for iterative parameter tuning. (Note: Using `--force_overwrite` completely bypasses this logic and re-runs everything).

---

## 12. Common workflows

### Run everything from scratch

```bash
Rscript qc_pipeline.R --config qc_config.dcf
```

### Re-run doublet detection only (QC already done)

```bash
Rscript qc_pipeline.R --config qc_config.dcf --stage doublet
```

### Test on a single dataset without touching your main output

```bash
Rscript qc_pipeline.R \
  --config qc_config.dcf \
  --rds_dir /path/to/single_sample \
  --output_dir /tmp/test_run
```

### Add a per-sample DBR override without re-running QC

1. Add the entry to `qc_config.dcf`:
   ```
   sample.my_sample.platform: dropseq
   sample.my_sample.dbr: 0.06
   ```
2. Run Stage 2 only:
   ```bash
   Rscript qc_pipeline.R --config qc_config.dcf --stage doublet
   ```

### Run with a different config for a second cohort

```bash
Rscript qc_pipeline.R --config cohort2_config.dcf
```

### Submit to a SLURM cluster

```bash
sbatch --job-name=qc_pipeline \
       --mem=64G \
       --cpus-per-task=4 \
       --time=04:00:00 \
       --wrap="Rscript qc_pipeline.R --config qc_config.dcf"
```

---

## 13. Troubleshooting

### `No RDS files found`

Check that `rds_dir` points to the correct directory and that files match `rds_pattern`. Try:

```bash
ls /your/rds_dir/*.rds
```

If files are in subdirectories, set `recursive_discovery: TRUE` in the config.

---

### `QC decisions CSV not found`

Verify the `decisions_csv` path. The CSV must exist before Stage 1 runs. Check that `Dataset_Name` values in the CSV match the RDS basenames in `rds_dir`.

---

### `No filtered RDS files found` (Stage 2)

Stage 2 reads from `filtered_subdir`. Either:
- Run `--stage all` or `--stage qc` first to produce the filtered files, or
- Make sure `filtered_subdir` in the config matches the directory where your pre-filtered RDS files are stored.

---

### scDblFinder size-factor warnings / very high doublet rates (~30 %)

Two separate causes:

1. **Ghost cells**: lower `dbl_min_count` and `dbl_min_feature` thresholds are not strict enough. Try increasing them (e.g. `dbl_min_count: 200`).
2. **Multi-sample objects processed as one pool**: ensure `dbl_batch_col` matches the metadata column that identifies individual samples. Confirm the column is present with `unique(obj$sample_id)` in R.

---

### `Non-10x platform requires sample.<n>.dbr`

You listed a sample with a non-10x platform but did not provide a `dbr` key. Add it to `qc_config.dcf`:

```
sample.<name>.platform: dropseq
sample.<name>.dbr: 0.05
```

---

### `Invalid DBR` error

The `dbr` value must be strictly between 0 and 0.30. A value of `0.08` means 8 % expected doublets. Values above 30 % are biologically implausible and rejected.

---

### Inflation ratio warning (`⚠ OVER-FILTERING`)

The median nCount_RNA increased by more than 1.5× after QC, meaning the lower nCount cutoff may be too aggressive and is removing genuine low-expression cells rather than empty droplets. Review the `lower_tri` or `knee` threshold for that dataset in the QC PNGs and consider using `manual_<value>` or `none` for the lower bound.

---

### Dashboard ridge plots missing

The ridge plot section requires the raw RDS files to still be accessible at the path in `rds_dir` after filtering (it reloads them briefly for the distribution overlay). If `rds_dir` is not available at dashboard-generation time, the pipeline falls back to a 3-panel dashboard and logs a warning.
