# Single-Cell QC Pipeline

**`QC_scdbl_Combined.R`** is the main script for automated doublet detection and quality control filtering of single-cell RNA-seq data. It handles everything from raw Seurat RDS files to clean, filtered outputs ready for downstream analysis. Both stages are driven by a single config file and write to a shared log, so you always know exactly what happened and why.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Repository Layout](#3-repository-layout)
4. [Quick Start](#4-quick-start)
5. [Configuration File Reference](#5-configuration-file-reference)
6. [CLI Reference](#6-cli-reference)
7. [Stage 1 - Doublet Detection](#7-stage-1---doublet-detection)
8. [Stage 2 - QC Filtering](#8-stage-2---qc-filtering)
9. [Integrated Summary](#9-integrated-summary)
10. [Outputs](#10-outputs)
11. [Threshold Method Reference](#11-threshold-method-reference)
12. [Automatic sc_protocol Detection and DBR Resolution](#12-automatic-sc_protocol-detection-and-dbr-resolution)
13. [Parameter-Aware Resume Logic](#13-parameter-aware-resume-logic)
14. [Common Workflows](#14-common-workflows)
15. [Troubleshooting](#15-troubleshooting)

---

## 1. Overview

The pipeline runs two stages in a fixed order, followed by an integrated summary:

```
Raw harmonized RDS files
      |
      v  Stage 1 - Doublet Detection (scDblFinder)
      |  - Auto-reads sc_protocol metadata from each RDS
      |  - Parses platform and chemistry (10x v2/v3/v3.1/v4/5',
      |    Drop-seq, PIPseq)
      |  - Ghost-cell filter (removes CellBender low-count artifacts)
      |  - Normalisation -> PCA -> UMAP
      |  - scDblFinder (batch-aware when multi-sample metadata present)
      |  - Per-sample DBR: auto-estimated from the appropriate 10x table,
      |    or looked up from platform defaults, or set manually in config
      |  - Singlet-only RDS + full-call RDS
      |  - 4-panel audit PDF per sample
      |  - Doublet summary CSV
      |
      v  Stage 2 - QC Filtering
      |  - Reads singlet-only RDS from Stage 1
      |  - Per-dataset config-driven thresholds (nFeature, nCount, %mito)
      |  - Thresholds defined inline in the config file (qc_decisions_table)
      |  - 16 QC PNGs per dataset (violin + scatter, raw + filtered)
      |  - Before/After audit PDF
      |  - Full 7-panel dashboard (ridge plots, inflation ratios, ...)
      |  - Filtered RDS files
      |
      v  Integrated Summary
         - 9-panel end-to-end dashboard (PNG + PDF)
         - Per-sample CSV tracking cell counts across every stage
         - Waterfall bar chart, funnel plot, density overlays
```

Doublet detection always runs first. This is intentional: removing doublets before QC filtering prevents them from inflating the QC metric distributions, which would throw off automatic thresholds.

Both stages read from the **same** `qc_config.dcf` file and write to the **same** log file.

---

## 2. Requirements

### Package Management and Runtime Environment

This pipeline uses Conda to manage all dependencies cleanly (R >= 4.2, Seurat, scDblFinder, and system libraries needed for PDF/PNG rendering such as Cairo and HarfBuzz).

The script does not install packages at runtime. This is a deliberate choice for cluster security and reproducibility.

### Step 1 - Create and activate the Conda environment

An `environment.yml` file is provided under `scripts/qc/`. Build and activate it before running:

```bash
# From the repository root, build the environment (first time only)
conda env create -f scripts/qc/environment.yml

# Activate it
conda activate hkoca_qc_env
```

### Full package list

| Package | Source | What it does |
|---|---|---|
| `r-base` | conda-forge | R runtime |
| `cairo`, `harfbuzz`, `fribidi` | conda-forge | System libraries for plotting |
| `r-seurat` | conda-forge | Seurat object handling |
| `bioconductor-singlecellexperiment` | bioconda | SCE conversion for scDblFinder |
| `bioconductor-scdblfinder` | bioconda | Doublet detection |
| `r-ggplot2` | conda-forge | All plotting |
| `r-patchwork` | conda-forge | Multi-panel plot composition |
| `r-cowplot` | conda-forge | Plot layout helpers |
| `r-ggridges` | conda-forge | Ridge plots in dashboards |
| `r-ggrepel` | conda-forge | Non-overlapping labels in integrated scatter plots |
| `r-reshape2` | conda-forge | Data reshaping for dashboard |
| `r-mass` | conda-forge | 2-D kernel density estimation |
| `r-scales` | conda-forge | Axis and colour helpers |
| `r-yaml` | conda-forge | Config parsing |
| `r-dplyr` | conda-forge | Data manipulation |
| `r-ggrastr` | conda-forge | Raster rendering for large scatters (>50k cells) |
| `r-jsonlite` | conda-forge | Parameter-aware resume checkpointing |
| `r-digest` | conda-forge | Cryptographic hashing for checkpoint comparison |

---

## 3. Repository Layout

```
scripts/test/qc/final/
|-- QC_scdbl_Combined.R                # Main pipeline script
|-- README_QC_scDblFinder.md           # This file
|-- qc_config.dcf                      # Configuration file (edit this)
|-- my_job.sh                          # Example SLURM submission script
|-- raw_rds/                           # Symlinked or actual input RDS files
+-- qc_results_<date>/                 # Output root (set via output_dir)
    |-- logs/
   |   +-- qc_pipeline.log
    |-- doublet_filtered_rds/          # Stage 1 output -> Stage 2 input
    |   |-- Summary/
    |   |   |-- Doublet_Audit_Report_<timestamp>.pdf
    |   |   +-- doublet_summary_<timestamp>.csv
   |   |   +-- raw_cell_counts_cache.csv
    |   |-- <name>_singlets.rds
    |   +-- <name>_with_doublet_calls.rds
    |-- Summary/                       # Stage 2 output
    |   |-- QC_Before_After_Report.pdf
    |   |-- summary_qc_full_dashboard.png
    |   |-- qc_summary_detailed.csv
    |   +-- <dataset_name>/
    |       |-- Violin_qc.raw.nFeature_RNA.png
    |       |-- Plot_qc.raw.nFeature_RNA.png
    |       +-- ...                    # 16 PNGs per dataset
    |-- qc_filtered_rds/               # Stage 2 filtered output
    |   +-- <name>_filtered.rds
   |-- integrated_qc_doublet_summary.csv    # Integrated summary (or under integrated_summary_subdir/)
   |-- integrated_summary_dashboard.png     # Integrated dashboard (or under integrated_summary_subdir/)
   +-- integrated_summary_dashboard.pdf
```

---

## 4. Quick Start

### Step 1 - Set up the environment

```bash
conda env create -f scripts/qc/environment.yml
conda activate hkoca_qc_env
```

### Step 2 - Edit the config file

Open `qc_config.dcf` and point it at your data:

```
rds_dir:    /path/to/your/raw_rds
output_dir: /path/to/your/output
```

The QC threshold decisions are embedded directly in the config file using the `qc_decisions_table` key (see [Section 5](#5-configuration-file-reference) for the format).

### Step 3 - Run the full pipeline

```bash
Rscript QC_scdbl_Combined.R --config qc_config.dcf
```

This runs doublet detection first, then QC filtering on the singlet-only output, and finishes with the integrated summary. That's it.

---

## 5. Configuration File Reference

The config uses DCF format (Debian Control File): one `key: value` pair per line, no quotes needed, comments start with `#`.

```
key: value
```

---

### Section 1 - Input Paths

| Key | Required | Description |
|---|---|---|
| `rds_dir` | Yes | Directory containing the raw Seurat RDS files to process |

---

### Section 2 - QC Threshold Decisions (embedded table)

Instead of a separate CSV, threshold decisions live directly in the config file using `qc_decisions_table`. The first line is the header, and each subsequent indented line is one dataset:

```
qc_decisions_table: Dataset_Name,Lower_Feature_Method,Lower_Count_Method,Upper_Feature_Method,Upper_Count_Method,Upper_Mito_Method
 dataset_A_harmonized.rds,lower_tri,lower_tri,mad3,mad3,renyi
 dataset_B_harmonized.rds,manual_1000,manual_1800,mad3,mad3,renyi
```

**Required columns:**

| Column | Description |
|---|---|
| `Dataset_Name` | RDS filename (basename, with or without `.rds`) |
| `Lower_Feature_Method` | Method for the lower nFeature_RNA threshold |
| `Lower_Count_Method` | Method for the lower nCount_RNA threshold |
| `Upper_Feature_Method` | Method for the upper nFeature_RNA threshold |
| `Upper_Count_Method` | Method for the upper nCount_RNA threshold |
| `Upper_Mito_Method` | Method for the upper percent.mito threshold |

See [Section 11](#11-threshold-method-reference) for all accepted method strings.

---

### Section 3 - Output Paths

| Key | Default | Description |
|---|---|---|
| `output_dir` | (required) | Root directory for all pipeline outputs |
| `summary_subdir` | `Summary` | Subdirectory for QC plots and CSVs |
| `filtered_subdir` | `qc_filtered_rds` | Subdirectory for QC-filtered RDS files |
| `doublet_subdir` | `doublet_filtered_rds` | Subdirectory for doublet outputs (this is also where Stage 2 reads its input) |
| `integrated_summary_subdir` | (output_dir root) | Subdirectory for integrated summary files; if blank, they go directly under `output_dir` |
| `log_file` | `logs/qc_pipeline.log` | Path to the shared log file (relative paths resolve from `output_dir`) |
| `reverse_mode` | `TRUE` | Controls integrated-summary ordering. Default reverse mode reports raw -> ghost filter -> doublet removal -> QC filtering |

---

### Section 4 - File Discovery

| Key | Default | Description |
|---|---|---|
| `rds_pattern` | `\.rds$` | Regex pattern for finding RDS files in `rds_dir` |
| `recursive_discovery` | `FALSE` | Set to `TRUE` to search `rds_dir` recursively |

---

### Section 5 - Dataset Exclusions

| Key | Default | Description |
|---|---|---|
| `skip_noisy_datasets` | `TRUE` | Whether to skip datasets listed in `noisy_datasets` |
| `noisy_datasets` | (empty) | Comma-separated list of RDS basenames to exclude |

---

### Section 6 - scDblFinder Parameters

| Key | Default | Description |
|---|---|---|
| `dbl_batch_col` | `sample_id` | Metadata column identifying individual samples. When multiple unique values are present, scDblFinder runs independently per sample (batch-aware mode) |
| `dbl_min_count` | `100` | Minimum UMI count to keep a cell before doublet detection. Cells below this are removed so they don't cause size-factor estimation failures |
| `dbl_min_feature` | `50` | Minimum number of detected genes (same purpose as `dbl_min_count`) |
| `dbl_umap_dims` | `20` | Number of PCA dimensions for the UMAP (used only in the audit plots) |
| `dbl_default_platform` | `10x` | Default sequencing platform when `sc_protocol` metadata is missing and there's no config override |
| `dbl_default_chemistry` | `v3` | Default 10x chemistry when not determinable from `sc_protocol` |

---

### Section 7 - Platform-Specific DBR Defaults (non-10x)

For platforms without a 10x-style multiplet rate table, you can set a fixed doublet rate:

```
platform_dbr.<platform_name>: <fraction>
```

**Pre-configured defaults:**

| Key | Value | Description |
|---|---|---|
| `platform_dbr.dropseq` | `0.049` | Default DBR for Drop-seq samples |
| `platform_dbr.pipseq` | `0.07` | Default DBR for PIPseq samples |

You can add more platforms or override existing values:

```
platform_dbr.indrops: 0.06
```

---

### Section 8 - Per-Sample Platform and DBR Overrides

Per-sample settings use a dotted-key convention. These always override auto-detected values from `sc_protocol`:

```
sample.<SAMPLE_NAME>.platform:  <platform>
sample.<SAMPLE_NAME>.chemistry: <chemistry>
sample.<SAMPLE_NAME>.dbr:       <numeric>
```

The `<SAMPLE_NAME>` must match the RDS filename after stripping suffixes (see below).

**Example:**

```
# 10x sample with explicit DBR (skips table estimation)
sample.kidney_organoid_rep1.platform: 10x
sample.kidney_organoid_rep1.chemistry: v4
sample.kidney_organoid_rep1.dbr: 0.08

# Drop-seq sample with explicit DBR (overrides the platform_dbr.dropseq default)
sample.dropseq_run_A.platform: dropseq
sample.dropseq_run_A.dbr: 0.05
```

Most samples don't need per-sample overrides at all. The pipeline auto-detects platform and chemistry from `sc_protocol` metadata (see [Section 12](#12-automatic-sc_protocol-detection-and-dbr-resolution)).

### Sample-Name Normalization

Before looking up per-sample config, the script normalizes names by stripping these suffixes (in order):

- `.rds`
- `_harmonized_filtered`
- `_filtered`
- `_with_doublet_calls`
- `_singlets`

So all of these map to the same config entry:

| Name seen by the pipeline | Config lookup key |
|---|---|
| `my_sample.rds` | `my_sample` |
| `my_sample_filtered.rds` | `my_sample` |
| `my_sample_harmonized_filtered.rds` | `my_sample` |
| `my_sample_with_doublet_calls.rds` | `my_sample` |
| `my_sample_singlets.rds` | `my_sample` |

This is why config entries keep working even when Stage 1 and Stage 2 refer to files with different suffixes.

### Dots in Sample Names

The config parser uses this pattern:

```
^sample\.(.+)\.(platform|chemistry|dbr)$
```

The sample name is matched greedily and the field is anchored at the end, so sample names containing dots (like `11H_neg-j38.denoised`) are handled correctly.

---

## 6. CLI Reference

```
Rscript QC_scdbl_Combined.R [--config PATH] [--stage STAGE] [OPTIONS]
```

### Flags

| Flag | Default | Description |
|---|---|---|
| `--config PATH` | `qc_config.dcf` (next to script) | Path to the configuration file |
| `--stage STR` | `all` | Which stage(s) to run: `all`, `qc`, or `doublet` |
| `--rds_dir PATH` | (from config) | Override raw RDS input directory |
| `--output_dir PATH` | (from config) | Override output base directory |
| `--summary_subdir NAME` | (from config) | Override QC summary subdirectory name |
| `--filtered_subdir NAME` | (from config) | Override filtered RDS subdirectory name |
| `--doublet_subdir NAME` | (from config) | Override doublet output subdirectory name |
| `--log_file PATH` | (from config) | Override log file path |
| `--recursive_discovery` | (from config) | Enable recursive RDS file search |
| `--rds_pattern REGEX` | (from config) | Override RDS filename regex |
| `--force_overwrite` | `FALSE` | Force reprocessing of everything (bypasses resume checkpoints) |
| `--help` | | Print usage and exit |

CLI flags always override the config file. This lets you use a shared base config and customize individual runs without editing the file.

### Examples

```bash
# Full pipeline with default config (doublet -> QC -> integrated summary)
Rscript QC_scdbl_Combined.R --config qc_config.dcf

# Doublet detection only (Stage 1)
Rscript QC_scdbl_Combined.R --config qc_config.dcf --stage doublet

# QC filtering only (Stage 2, requires doublet-filtered RDS from Stage 1)
Rscript QC_scdbl_Combined.R --config qc_config.dcf --stage qc

# Override output directory at runtime
Rscript QC_scdbl_Combined.R --config qc_config.dcf --output_dir /scratch/my_run

# Point to a different RDS directory without editing the config
Rscript QC_scdbl_Combined.R --config qc_config.dcf --rds_dir /data/new_batch

# Run in background and redirect output
nohup Rscript QC_scdbl_Combined.R --config qc_config.dcf > run.log 2>&1 &
```

---

## 7. Stage 1 - Doublet Detection

### Input

- Raw Seurat RDS files from `rds_dir`

### What happens for each sample

1. **Load** the Seurat object with `readRDS()`.
2. **Ghost-cell filter** - removes cells with `nCount_RNA <= dbl_min_count` or `nFeature_RNA <= dbl_min_feature`. These leftover low-quality cells (common after CellBender) can cause size-factor estimation to blow up inside scDblFinder.
3. **S4/SXP fix** - coerces `meta.data` to a plain `data.frame` to work around a Seurat v5 S4 slot incompatibility with `as.SingleCellExperiment()`.
4. **Platform and chemistry auto-detection** - reads `sc_protocol` from the object's metadata and parses it into a normalised platform and chemistry (see [Section 12](#12-automatic-sc_protocol-detection-and-dbr-resolution)). Config overrides win if specified.
5. **Doublet rate (DBR) resolution** - resolves the expected doublet rate through a priority chain:
   - User override - explicit `sample.<name>.dbr` in config
   - 10x table estimation - for 10x platforms, piecewise-linear interpolation on the published multiplet rate table, matched to the right chemistry (v2/v3, v4, or 5')
   - Platform default - for non-10x platforms, looked up from `platform_dbr.*` config keys
   - Fallback - if nothing else works, let scDblFinder auto-estimate
6. **Normalisation** - `NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA -> RunUMAP` (UMAP is only used for the audit plot, not for filtering).
7. **scDblFinder** runs in one of two modes:
   - **Batch mode**: when `dbl_batch_col` has more than one unique value, each sub-sample is processed independently. This prevents the spurious ~30% doublet rate you'd get by treating a merged multi-sample object as a single pool.
   - **Single mode**: standard run on the whole object.
8. **Result mapping** - `scDblFinder.class` (singlet/doublet) and `scDblFinder.score` are added to the Seurat metadata.
9. **Audit plots** - a 4-panel figure is appended to the audit PDF: UMAP coloured by class, complexity scatter, feature violin, score density.
10. **Saves** two RDS files per sample:
    - `<name>_with_doublet_calls.rds` - full object with doublet annotations
    - `<name>_singlets.rds` - singlets only (this is what Stage 2 reads)

### 10x Multiplet Rate Tables

The pipeline ships with three built-in tables from 10x Genomics user guides:

**v2/v3 chemistry** (also covers v3.1) - [source PDF, page 18](https://cdn.10xgenomics.com/image/upload/v1722285481/support-documents/CG000315_ChromiumNextGEMSingleCell3__GeneExpression_v3.1_DualIndex__RevF.pdf):

| Loaded cells | Recovered | Multiplet rate |
|---:|---:|---:|
| 800 | 500 | 0.4% |
| 1,600 | 1,000 | 0.8% |
| 3,200 | 2,000 | 1.6% |
| 4,800 | 3,000 | 2.3% |
| 6,400 | 4,000 | 3.1% |
| 8,000 | 5,000 | 3.9% |
| 9,600 | 6,000 | 4.6% |
| 11,200 | 7,000 | 5.4% |
| 12,800 | 8,000 | 6.1% |
| 14,400 | 9,000 | 6.9% |
| 16,000 | 10,000 | 7.6% |

**v4 chemistry** (Single Cell 3' Reagent Kits v4, CG000731 Rev B - roughly 0.58x the v3 rate) - [source PDF, page 19](https://cdn.10xgenomics.com/image/upload/v1725314293/support-documents/CG000731_ChromiumGEM-X_SingleCell3v4_UserGuide_RevB.pdf):

| Loaded cells | Recovered | Multiplet rate |
|---:|---:|---:|
| 725 | 500 | 0.2% |
| 2,900 | 2,000 | 0.8% |
| 5,800 | 4,000 | 1.6% |
| 14,500 | 10,000 | 4.0% |
| 29,000 | 20,000 | 8.0% |

**5' chemistry** (Single Cell 5' Reagent Kit v3, CG000733 Rev A - same curve as v4) - [source PDF, page 21](https://cdn.10xgenomics.com/image/upload/v1710231087/support-documents/CG000733_ChromiumGEM-X_SingleCell5_ReagentKitsv3_UserGuide_RevA.pdf):

| Loaded cells | Recovered | Multiplet rate |
|---:|---:|---:|
| 725 | 500 | 0.2% |
| 2,900 | 2,000 | 0.8% |
| 5,800 | 4,000 | 1.6% |
| 14,500 | 10,000 | 4.0% |
| 29,000 | 20,000 | 8.0% |

For cell counts within the table range, the DBR is computed by piecewise-linear interpolation. For cell counts beyond the table range, the edge segment slope is extrapolated, clamped to [0, 0.30].

---

## 8. Stage 2 - QC Filtering

### Input

- Singlet-only RDS files from `doublet_subdir` (produced by Stage 1, or placed there beforehand if you run `--stage qc` directly).
- Threshold decisions from the `qc_decisions_table` key in the config.

`--stage qc` still reads from `doublet_subdir`; it does not read directly from `rds_dir`.

### What happens for each dataset

1. **Load** the Seurat object, rebuild it from scratch to avoid Seurat v5 update conflicts.
2. **Compute** QC metrics if missing: `nCount_RNA`, `nFeature_RNA`, `percent.mito`, `percent.ribo`.
3. **Remove** barcodes with zero counts (empty droplets that slipped through).
4. **Compute thresholds** by applying the algorithm named in each config column to the metric distribution.
5. **Filter** cells falling outside the threshold windows.
6. **Generate 16 PNGs** per dataset (4 features x 2 plot types x raw/filtered):
   - `Violin_qc.raw.<feature>.png` - violin with cutoff lines
   - `Plot_qc.raw.<feature>.png` - density scatter with cutoff lines
   - `Violin_qc.filtered.<feature>.png` - violin after filtering
   - `Plot_qc.filtered.<feature>.png` - scatter after filtering
7. **Append two pages** to the audit PDF (before and after density plots). This is the primary document for assessing how well filtering went: comparing the BEFORE and AFTER distributions of nFeature_RNA, nCount_RNA, and percent.mito side by side makes it immediately obvious whether the thresholds landed in the right places or need adjustment.
8. **Save** the filtered object as `<name>_filtered.rds` in `filtered_subdir`.
9. **Log** an inflation ratio diagnostic. If the median nCount jumps by more than 1.5x after filtering, a warning is emitted. This is a sign that the lower cutoff might be too aggressive.

### Dashboard

After all datasets are processed, a multi-panel `summary_qc_full_dashboard.png` is produced:

| Panel | Content |
|---|---|
| P1 | Cell retention bar chart (before vs after) |
| P2 | Filtering intensity (% cells removed) with a 35% warning line |
| P3 | Median nCount inflation ratio with thresholds at 1.0x, 1.4x, 1.5x |
| P4 | Median nCount_RNA before vs after (line plot) |
| P5 | Median nFeature_RNA before vs after (line plot) |
| P6 | Median %mito and %ribo before vs after (facet bar) |
| P7a | Per-dataset log1p(nCount) density with cutoff lines |
| P7b | Per-dataset log1p(nFeature) density with cutoff lines |
| P8a | Ridge plot: log1p(nCount) across all datasets |
| P8b | Ridge plot: log1p(nFeature) across all datasets |

---

## 9. Integrated Summary

When both stages have run (or their output files already exist on disk), the pipeline produces an end-to-end summary that tracks every cell through the full filtering funnel.

### Output Files

| File | Description |
|---|---|
| `integrated_qc_doublet_summary.csv` | One row per sample with cell counts and loss percentages at each stage |
| `integrated_summary_dashboard.png` | Composite integrated summary figure |
| `integrated_summary_dashboard.pdf` | Same figure in PDF format |

If `integrated_summary_subdir` is set in the config, these files are written there; otherwise they are written directly under `output_dir`.

### CSV Columns

| Column | Description |
|---|---|
| `sample` | Sample identifier (canonicalized) |
| `cells_raw` | Total cells in the raw RDS |
| `cells_after_ghost` | Cells remaining after the ghost-cell filter |
| `dbr_used` | Doublet rate used by scDblFinder |
| `doublets_found` | Number of cells classified as doublets |
| `pct_doublets` | Percentage classified as doublets |
| `cells_after_doublet_removal` | Cells remaining after doublet removal |
| `cells_after_qc` | Cells remaining after QC filtering |
| `pct_qc_removal` | % removed by QC (relative to post-doublet count) |
| `total_cells_removal` | Total cells removed across both stages |
| `pct_total_cells_removal` | Overall % of raw cells removed |

### Dashboard Panels

| Panel | Content |
|---|---|
| 1 | Reverse-mode waterfall chart (raw -> post-ghost -> post-doublet -> post-QC) |
| 2 | Stacked % loss breakdown (doublet vs QC) |
| 3 | Cumulative funnel line plot across the same reverse sequence |
| 4 | log1p(nCount_RNA) density overlay using raw, post-doublet, and post-QC files |
| 5 | log1p(nFeature_RNA) density overlay |
| 6 | % mito overlay (raw vs post-QC) |
| 7 | Doublet score distribution per sample (ridge plot) |
| 8 | Doublet % vs QC removal % scatter |
| 9 | Total cells lost per stage (pie/donut) |

By default the script builds the integrated summary in `reverse_mode`, which matches the current execution order: doublet detection first, then QC. If `reverse_mode` is set to `FALSE`, the integrated summary falls back to the legacy raw -> QC -> doublet interpretation.

---

## 10. Outputs

### Stage 1 - Doublet Detection

| Path | Description |
|---|---|
| `<output_dir>/<doublet_subdir>/Summary/Doublet_Audit_Report_<ts>.pdf` | 4-panel audit plots per sample |
| `<output_dir>/<doublet_subdir>/Summary/doublet_summary_<ts>.csv` | Per-sample doublet counts and rates |
| `<output_dir>/<doublet_subdir>/<name>_with_doublet_calls.rds` | Full object with `scDblFinder.class` and `scDblFinder.score` |
| `<output_dir>/<doublet_subdir>/<name>_singlets.rds` | Singlets-only Seurat object (input for Stage 2) |

### Stage 2 - QC Filtering

| Path | Description |
|---|---|
| `<output_dir>/<summary_subdir>/QC_Before_After_Report.pdf` | Density plots before and after filtering |
| `<output_dir>/<summary_subdir>/summary_qc_full_dashboard.png` | 7-panel composite dashboard |
| `<output_dir>/<summary_subdir>/qc_summary_detailed.csv` | Per-dataset stats table (40+ columns: cell counts, medians, cutoffs, inflation ratio) |
| `<output_dir>/<summary_subdir>/<dataset>/Violin_qc.*.png` | 8 violin PNGs per dataset |
| `<output_dir>/<summary_subdir>/<dataset>/Plot_qc.*.png` | 8 scatter PNGs per dataset |
| `<output_dir>/<filtered_subdir>/<name>_filtered.rds` | Final QC-filtered Seurat object |

### Integrated Summary

| Path | Description |
|---|---|
| `<output_dir>/<integrated_summary_subdir>/integrated_qc_doublet_summary.csv` | Per-sample end-to-end cell counts (or directly under `output_dir` if the subdir is blank) |
| `<output_dir>/<integrated_summary_subdir>/integrated_summary_dashboard.png` | Integrated composite dashboard |
| `<output_dir>/<integrated_summary_subdir>/integrated_summary_dashboard.pdf` | Same dashboard in PDF |

### Shared

| Path | Description |
|---|---|
| `<output_dir>/logs/qc_pipeline.log` | Timestamped log covering all stages |

---

## 11. Threshold Method Reference

These strings go in the `qc_decisions_table` column values:

| Method string | Algorithm | Notes |
|---|---|---|
| `lower_tri` | Triangle threshold (lower tail) | Finds the histogram bin at maximum perpendicular distance from the line connecting the peak to the left endpoint |
| `upper_tri` | Triangle threshold (upper tail) | Same idea, but peak to right endpoint |
| `renyi` | Renyi entropy binarisation | Maximises the sum of Renyi entropies of the two histogram partitions |
| `knee` | Kneedle algorithm | Fits a curve to ranked values and finds the point of maximum curvature |
| `mad3` | Median +/- 3 x MAD | For upper columns: median + 3xMAD; for lower columns: median - 3xMAD |
| `none` | No threshold | Sets the cutoff to `Inf` (upper) or `-Inf` (lower), effectively no filtering |
| `manual_<value>` | Fixed numeric | Sets the cutoff directly, e.g. `manual_200` means the threshold is 200 |

### A note on choosing thresholds

`lower_tri` works well for most datasets because it automatically adapts to the data distribution. If you know the biology well and want precise control, use `manual_<value>`. Be careful with very aggressive lower bounds since most single-cell data has a median of 2,000-5,000 genes per cell, and you could accidentally filter out everything.

---

## 12. Automatic sc_protocol Detection and DBR Resolution

The pipeline automatically reads the `sc_protocol` column from each RDS object's metadata to determine the sequencing platform and chemistry. For most datasets, you won't need to configure `sample.<name>.platform` or `sample.<name>.chemistry` at all.

### How it works

When processing each sample in Stage 1:

1. **Read metadata** - the pipeline checks for `obj@meta.data$sc_protocol`. If the column exists and has a single unique non-NA value, it's used. If there are multiple values, the first one is used (with a warning).
2. **Parse the protocol string** - it's normalised into a platform and chemistry:

| `sc_protocol` value | Parsed platform | Parsed chemistry |
|---|---|---|
| `10x_3_v3` | `10x` | `v3` |
| `10x_3_v3.1` | `10x` | `v3` (v3.1 shares the v3 multiplet table) |
| `10x_3_v2` / `10x_v2` | `10x` | `v2` |
| `10x_3_v1` | `10x` | `v2` (v1 is closest to the v2 table) |
| `10x_v4` | `10x` | `v4` |
| `10X_5` | `10x` | `5p` (dedicated 5' multiplet rate table) |
| `Dropseq` / `Drop-seq` | `dropseq` | N/A |
| `PIPseq_V_T10_3` | `pipseq` | N/A |

3. **Apply config override** - if a `sample.<name>.platform` or `sample.<name>.chemistry` key exists in the config, it wins over the auto-detected value.
4. **Resolve DBR** - the final priority chain:

```
 +-------------------------------------------------------------+
 |  1. User override (sample.<name>.dbr in config)              |
 |  2. 10x table estimation (for 10x platforms)                 |
 |     -- OR --                                                 |
 |  2. Platform default (platform_dbr.<name> in config)         |
 |     for non-10x platforms (Drop-seq -> 0.049, PIPseq -> 0.07)|
 |  3. scDblFinder auto-estimation (last resort)                |
 +-------------------------------------------------------------+
```

### DBR source labels in logs

The log records which resolution path was used, so you can always trace why a particular rate was chosen:

| `DBR_Source` | Meaning |
|---|---|
| `SC_PROTOCOL_AUTO` | Platform/chemistry parsed from `sc_protocol` metadata |
| `TENX_V3_TABLE_FROM_FILTERED_CELLS` | DBR estimated from the 10x v3 multiplet table |
| `TENX_V4_TABLE_FROM_FILTERED_CELLS` | DBR estimated from the 10x v4 multiplet table |
| `TENX_5P_TABLE_FROM_FILTERED_CELLS` | DBR estimated from the 10x 5' multiplet table |
| `TENX_V*_TABLE_PER_BATCH_FILTERED_CELLS` | For merged 10x inputs, DBR was estimated separately for each internal batch |
| `PLATFORM_DEFAULT_DROPSEQ` | Used `platform_dbr.dropseq` from config (0.049) |
| `PLATFORM_DEFAULT_PIPSEQ` | Used `platform_dbr.pipseq` from config (0.07) |
| `USER_OVERRIDE_10X` | Explicit `sample.<name>.dbr` for a 10x sample |
| `USER_OVERRIDE` | Explicit `sample.<name>.dbr` for a non-10x sample |
| `AUTO_10X_FALLBACK` | All other methods failed, scDblFinder auto-estimates |

The summary CSV can also contain run-state labels such as `SKIPPED_RESUME` for reused outputs.

### When sc_protocol is missing

If the `sc_protocol` column isn't in a sample's metadata, the pipeline falls back to:
1. Per-sample config overrides (`sample.<name>.platform` / `.chemistry`)
2. Global defaults (`dbl_default_platform` and `dbl_default_chemistry`)

A warning is logged so you know it happened.

---

## 13. Parameter-Aware Resume Logic

The pipeline uses a smart checkpoint system so you can re-run it without reprocessing samples that haven't changed.

### How it works

1. When a dataset finishes in Stage 1 or Stage 2, the script saves a tiny hidden `.json` sidecar file next to the stage output RDS.
2. This JSON file contains a cryptographic hash (via `digest`) of the exact parameters used to produce that output.
3. On the next run, the pipeline compares the current config against each sidecar:
   - **Match**: the dataset was already processed with the same thresholds, so it's skipped. Its metrics are still picked up for the summary dashboard.
   - **Mismatch**: you changed even one threshold for that sample, so it gets reprocessed automatically. Everything else is still skipped.

For older doublet outputs that predate sidecars, the script also supports a legacy file-existence fallback so previous results can still be reused.

This makes iterative threshold tuning painless. You tweak one sample's config line, re-run, and only that sample gets reprocessed.

Using `--force_overwrite` bypasses this entirely and reprocesses everything.

---

## 14. Common Workflows

### Run everything from scratch

```bash
Rscript QC_scdbl_Combined.R --config qc_config.dcf
```

### Run doublet detection only (Stage 1)

```bash
Rscript QC_scdbl_Combined.R --config qc_config.dcf --stage doublet
```

### Re-run QC only after tweaking thresholds (doublets already done)

```bash
Rscript QC_scdbl_Combined.R --config qc_config.dcf --stage qc
```

The resume system will only reprocess samples whose thresholds changed. Make sure the singlet RDS files already exist in `doublet_subdir`, because QC-only mode reads from that location.

### Test on a single dataset without touching your main output

```bash
Rscript QC_scdbl_Combined.R \
  --config qc_config.dcf \
  --rds_dir /path/to/single_sample \
  --output_dir /tmp/test_run
```

### Add a platform-wide DBR for a new platform

Add this line to `qc_config.dcf`:

```
platform_dbr.indrops: 0.06
```

No per-sample overrides needed. The pipeline auto-detects the platform from `sc_protocol` and applies this DBR.

### Override a specific sample's DBR

```
sample.my_special_sample.dbr: 0.10
```

Then re-run with `--stage doublet` or `--stage all`.

### Submit to a SLURM cluster

```bash
sbatch --job-name=qc_pipeline \
       --mem=64G \
       --cpus-per-task=4 \
       --time=04:00:00 \
       --wrap="Rscript QC_scdbl_Combined.R --config qc_config.dcf"
```

---

## 15. Troubleshooting

### "No RDS files found"

Make sure `rds_dir` points to the right directory and that files match `rds_pattern`. A quick check:

```bash
ls /your/rds_dir/*.rds
```

If files are in subdirectories, set `recursive_discovery: TRUE`.

---

### "No filtered RDS files found" (Stage 2)

Stage 2 reads from `doublet_subdir` (the singlet-only output of Stage 1). Either:
- Run `--stage all` or `--stage doublet` first to produce the singlet files, or
- Check that `doublet_subdir` in the config matches where your pre-filtered RDS files actually are.

---

### "No cells found" error

This usually means the lower thresholds in `qc_decisions_table` are too aggressive. For example, `manual_12000` as a lower feature threshold requires over 12,000 genes per cell, but most scRNA-seq cells only have 2,000-5,000 genes, so every cell gets filtered out.

Fix: change the lower thresholds to something reasonable like `lower_tri` (auto-adaptive) or `manual_1200`.

---

### "sc_protocol column not found" warnings

The pipeline still works without `sc_protocol` - it falls back to global defaults (`dbl_default_platform`, `dbl_default_chemistry`). But for the most accurate DBR estimation, add the `sc_protocol` column to your Seurat objects during harmonisation. Accepted values:

```
10x_3_v3, 10x_3_v3.1, 10x_3_v2, 10x_3_v1, 10x_v4, 10X_5,
Dropseq, Drop-seq, PIPseq_V_T10_3
```

---

### scDblFinder warnings or very high doublet rates (~30%)

Two common causes:

1. **Ghost cells**: the `dbl_min_count` and `dbl_min_feature` thresholds aren't strict enough. Try bumping them up (e.g. `dbl_min_count: 200`).
2. **Multi-sample objects processed as one pool**: make sure `dbl_batch_col` matches the metadata column that identifies individual samples. Verify with `unique(obj$sample_id)` in R.

---

### "Invalid DBR" error

The `dbr` value must be between 0 and 0.30. A value of `0.08` means 8% expected doublets. Anything above 30% is biologically implausible and gets rejected.

---

### Non-10x platform with no DBR available

If the pipeline encounters a platform it doesn't recognise, it falls back to 10x auto-estimation and logs a warning. To fix:

1. Add a platform default: `platform_dbr.myplatform: 0.06`
2. Or provide a per-sample override: `sample.my_sample.dbr: 0.06`

---

### Inflation ratio warning ("OVER-FILTERING")

The median nCount_RNA increased by more than 1.5x after QC. This means the lower nCount cutoff is probably too aggressive - it's removing genuine low-expression cells instead of just empty droplets. Check the QC PNGs for that dataset, and consider using `lower_tri`, `manual_<value>` with a lower number, or `none` for the lower bound.

---

### Dashboard ridge plots missing

The ridge plot section needs access to the raw RDS files at `rds_dir` after filtering (it briefly reloads them for the distribution overlay). If `rds_dir` is unavailable at dashboard time, the pipeline falls back to a simpler dashboard and logs a warning.
