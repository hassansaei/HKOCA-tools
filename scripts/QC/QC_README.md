# Single-Cell RNA-Seq QC, Doublet Detection, and H5AD Conversion Pipeline

## 1. Project Title & Overview

### 1.1 Executive Summary

This repository contains a single, unified R pipeline for processing raw single-cell RNA sequencing (scRNA-seq) data stored as Seurat `.rds` objects, cleaning it through doublet detection and quality-control filtering, and exporting the result directly to the `.h5ad` (AnnData) format for downstream use in the Python single-cell ecosystem (Scanpy, scVI, etc.).

| Component | Role |
|---|---|
| **`QC_scdbl_Combined.R`** | A configuration-driven, multi-stage analytical engine. It performs **doublet detection** (via `scDblFinder`), **quality control filtering** (nFeature/nCount/percent-mitochondrial thresholding), and — as of this revision — **native `.h5ad` conversion**, all from a single script and a single run. |
| **`qc_config.dcf`** | The single DCF configuration file that drives every stage of the pipeline: input/output paths, doublet tuning, per-sample platform/chemistry/DBR overrides, QC filter criteria, and H5AD conversion settings. |
| **`environment_qc.yaml`** | A Conda environment specification providing the exact R, Python, CRAN, Bioconductor, and Python-package versions required to run the pipeline end-to-end, including the Python `anndata` package consumed via `reticulate`. |

**The problem this pipeline solves:** Raw single-cell count matrices are contaminated by technical artifacts that, if left unaddressed, corrupt downstream biological inference. These artifacts include empty droplets, ambient RNA "ghost cells," dying or stressed cells (high mitochondrial content), and droplets containing two or more cells (doublets), which is common in droplet-based platforms such as 10x Genomics. Manually identifying clean, sample-specific thresholds for hundreds of datasets is labor-intensive, error-prone, and difficult to audit. Compounding this, single-cell labs typically work in a mixed R/Python environment, and previously this pipeline required a **second, separate R script** (`rds_to_h5ad.r`, now retired) invoking the external `SeuratDisk`/`hdf5r` toolchain to bridge that gap.

**What changed in this revision:** The RDS→h5ad conversion logic has been **fully absorbed into `QC_scdbl_Combined.R`** as a new, optional pipeline stage (Section 6.5, "H5AD Conversion"). It no longer depends on `SeuratDisk` or `hdf5r`; instead it constructs an `AnnData` object natively in R via the `anndata` and `reticulate` packages, sourcing counts and cell metadata directly from the same Seurat object already in memory during the pipeline run. Conversion is now driven by the **same DCF configuration file** as the rest of the pipeline (`run_h5ad_conversion` / `h5ad_output_subdir` keys), eliminating the manual synchronization step that the old two-script design required. A dedicated `environment_qc.yaml` is now provided to reproducibly provision every dependency this consolidated pipeline needs, including the Python-side `anndata` and `umap-learn` packages that `reticulate` calls into from R.

---

## 2. Architecture & Interaction

### 2.1 High-Level Workflow

The pipeline is now a **single script, single process, single configuration file** design. All three stages execute sequentially within one `Rscript` invocation (subject to the `--stage` flag and the `run_h5ad_conversion` config key described below).

```
                         ┌─────────────────────────────┐
                         │   Raw Seurat .rds files      │
                         │   (one file per sample)      │
                         └───────────────┬─────────────┘
                                         │
                                         ▼
                ┌────────────────────────────────────────────┐
                │              QC_scdbl_Combined.R             │
                │        (single script, single config)        │
                │                                              │
                │  STAGE 1 — Doublet Detection (scDblFinder)  │
                │   • Ghost-cell pre-filter                    │
                │   • Platform/chemistry-aware DBR estimation  │
                │   • Ambiguous-zone (KDE overlap) detection   │
                │   • Outputs: *_singlets.rds, audit PDF, CSV  │
                │                       │                       │
                │                       ▼                       │
                │  STAGE 2 — QC Filtering                      │
                │   • Noisy-dataset exclusion (skip list)      │
                │   • nFeature / nCount / %mito thresholding   │
                │   • 5 selectable threshold algorithms        │
                │   • Outputs: *_filtered.rds, audit PDF, CSV  │
                │                       │                       │
                │                       ▼                       │
                │  STAGE 3 — Integrated Summary                │
                │   • Combined funnel CSV/PNG/PDF dashboards   │
                └────────────────────────┬─────────────────────┘
                                         │
                                         ▼
                ┌────────────────────────────────────────────┐
                │  STAGE 4 (Optional) — H5AD Conversion        │
                │   • Enabled via run_h5ad_conversion: TRUE    │
                │   • Native anndata/reticulate construction   │
                │   • No intermediate h5Seurat file            │
                │   • Outputs: *.h5ad                          │
                │                                              │
                └────────────────────────┬─────────────────────┘
                                         │
                                         ▼
                         Python / Scanpy / AnnData
                          downstream analysis
```

### 2.2 Retirement of the Standalone Conversion Script

Previous revisions of this pipeline shipped a second script, `rds_to_h5ad.r`, which was invoked separately after `QC_scdbl_Combined.R` completed, and which relied on `SeuratDisk::SaveH5Seurat()` / `SeuratDisk::Convert()` plus a manual `hdf5r`-based metadata parity check. **This script is now retired and superseded** by the integrated Section 6.5 stage described below.

> **Migration Note:** If your environment or downstream tooling still references `rds_to_h5ad.r` or `INPUT_DIR`/`OUTPUT_DIR` constants from that script, those integration points should be removed. The new `run_h5ad_conversion` and `h5ad_output_subdir` keys in `qc_config.dcf` now fully control h5ad output, and no second script invocation or manual path synchronization is required.

### 2.3 Internal Architecture of `QC_scdbl_Combined.R`

The script is organized into numbered sections, executing in this order:

| Section | Name | Responsibility |
|---|---|---|
| 0 | Package Loading | Verifies all required R packages — now including `anndata` and `reticulate` — are installed; fails fast with a clear, actionable error directing the operator to build the environment via `conda env create -f environment_qc.yaml`. |
| 1 | Shared Runtime Utilities | CLI parsing, DCF config reading, type coercion, path resolution, logging, RDS file discovery, 10x Genomics multiplet-rate lookup tables, parameter-hash-based skip/resume logic, three-tier QC filter resolution (`--filters` CSV → inline `qc_decisions_table` → auto-generated MAD3 fallback), and small-sample safety guards for `scDblFinder`. |
| 2 | Configuration & Path Setup | Resolves the active configuration (file + CLI overrides), determines which stage(s) to run, resolves the new H5AD output directory and flag (`RUN_H5AD`, `H5AD_DIR`), and creates the output directory tree. |
| 3 | Threshold Algorithm Implementations | Five interchangeable statistical methods for deriving QC cutoffs (`lower_tri`, `upper_tri`, `renyi`, `knee`, `mad3`, `none`, `manual_<value>`). |
| 4 | QC Plot Builders | Reusable `ggplot2` builder functions for density-colored scatterplots, violin plots (global and per-sample-faceted), and before/after distribution overlays. |
| 5 | Stage 2 — QC Filtering | The main QC loop: resolves QC filter criteria, optionally excludes datasets listed under `noisy_datasets`, loads each raw/doublet-filtered RDS, computes thresholds, filters cells, generates per-sample plot sets, and writes a detailed summary CSV plus a combined dashboard PDF/PNG. **Now includes native Seurat v5 Compatibility** (see callout below). |
| 6 | Stage 1 — Doublet Detection | The main doublet loop: ghost-cell pre-filtering, platform/chemistry-aware DBR resolution, `scDblFinder` execution with graceful fallbacks, ambiguous-zone classification, PCA/UMAP visualization, and an audit PDF/summary CSV. **Now includes Dynamic Mixed Chemistry Batching** (see callout below). |
| **6.5** | **Stage 3 — H5AD Conversion (new)** | Optional, config-gated stage. Iterates the appropriate upstream output directory, reads each `.rds` file, casts any Seurat v5 `Assay5` to a legacy `Assay`, extracts the raw counts matrix (transposed to Cells × Genes) and cell metadata, constructs an `AnnData` object natively via the `anndata` R package, and writes it directly to `.h5ad`. Skips files whose `.h5ad` output already exists unless `force_overwrite` is set. |
| 7 | Pipeline Complete | Final console/log summary of stage outcomes and output file locations. |
| 8 | Integrated Summary | Merges the QC and doublet summary CSVs into a single per-sample funnel table (`raw → ghost-filtered → doublet-filtered → QC-filtered`) with four dashboard visualizations (stacked loss bars, cumulative funnel, doublet-vs-QC scatter, cohort donut chart). |
| 9 | Final Pipeline Footer | Closing log banner. |

**Default execution order** (`--stage all`): the **Doublet Detection stage runs first**, consuming the raw RDS directory; its singlet output then becomes the input to the **QC Filtering stage**. If `run_h5ad_conversion: TRUE` is set in the config, the **H5AD Conversion stage runs immediately afterward**, regardless of the `--stage` value, converting whichever directory represents the most recently completed stage's output (see Section 5.1.2 for the exact input-directory resolution rule).

> ** New Capability — Dynamic Mixed Chemistry Batching (Stage 1):** The doublet detection stage now robustly handles cohorts where a single dataset contains cells generated under **more than one 10x chemistry** (e.g., a merged object containing both v2 and v3 cells under one `sample_id`). When multiple `sc_protocol` values are detected within the same batch, the pipeline **automatically forces internal batch-mode** and resolves the mathematically correct multiplet-rate (DBR) table **per sub-batch, on the fly** — rather than requiring the operator to manually split the dataset or apply a single averaged DBR across a mixed-chemistry sample. This removes a previously manual, error-prone pre-processing step for cohorts built from merged runs.

> ** New Capability — Native Seurat v5 Compatibility (Stage 2):** The QC filtering stage now safely handles Seurat v5 objects whose counts are stored across **fragmented `@layers`** (e.g., multiple `counts.<sample>` layers from a merged, unjoined object) rather than a single monolithic matrix. The pipeline automatically detects this layout and **rebuilds a unified barcode × feature matrix on the fly** before applying thresholds, preventing the silent data loss or hard crashes that fragmented v5 layers caused in earlier pipeline versions and in other, less v5-aware QC pipelines.

---

## 3. Prerequisites & Requirements

### 3.1 Operating System

| OS | Support Status |
|---|---|
| Linux (Ubuntu 20.04+ / CentOS 7+ / Debian 11+) | **Recommended.** All Bioconductor and Python/reticulate binaries are most reliably built here. |
| macOS (Intel & Apple Silicon) | Supported, provided Xcode Command Line Tools are installed and a compatible Conda/Miniforge distribution is used for environment creation. |
| Windows (10/11) | Supported only via WSL2. The `reticulate`-based Python interop that the H5AD conversion stage now depends on is substantially more reliable under Linux/WSL2 than native Windows R. |

### 3.2 Runtime & Language Versions

The pipeline is now explicitly **dual-language** (R + Python), coordinated through `reticulate`. Both language runtimes are pinned in the provided `environment_qc.yaml`.

| Component | Minimum Version | Notes |
|---|---|---|
| R | 4.2.0 | Both the analytical stages and the new native H5AD conversion require modern R syntax (`tryCatch`, S4 dispatch on `Assay5`). |
| Python | 3.10 | Consumed transparently through `reticulate` by the R `anndata` package during the H5AD conversion stage. The operator does not invoke Python directly. |
| Bioconductor | 3.16+ | Required to install `SingleCellExperiment` and `scDblFinder`. |

### 3.3 Required Packages

All required packages are declared in `environment_qc.yaml` and are verified at runtime by the script's `load_packages_safely()` startup check.

#### R (CRAN)

| Package | Purpose |
|---|---|
| `Seurat` (≥4.3.0) | Core single-cell object model, normalization, PCA, UMAP. |
| `MASS` | 2-D kernel density estimation (`kde2d`) for scatterplot point-density coloring; also `mad()` for the `mad3` threshold method. |
| `ggplot2`, `patchwork`, `cowplot`, `scales`, `reshape2`, `ggridges`, `ggrepel` | Plotting, multi-panel layout composition, and axis/theming utilities across all audit PDFs and dashboards. |
| `yaml`, `dplyr`, `jsonlite`, `digest` | Configuration parsing, data manipulation, and parameter-hash checkpointing (`.meta.json` files). |
| **`anndata`** | **New.** Provides the R-native `AnnData()` constructor and `$write_h5ad()` method used by the integrated H5AD conversion stage (Section 6.5). |
| **`reticulate`** | **New.** Bridges R and the underlying Python `anndata` package; required by the R `anndata` package to function. |

#### R (Bioconductor)

| Package | Purpose |
|---|---|
| `SingleCellExperiment` | The object class consumed by `scDblFinder`. |
| `scDblFinder` | The core doublet detection algorithm. |

#### Python (consumed via `reticulate`, not invoked directly by the operator)

| Package | Purpose |
|---|---|
| `anndata` | The underlying Python AnnData implementation that the R `anndata` package wraps. |
| `umap-learn` | Provides UMAP support consumed transparently through the R/Python interop layer. |

> **Optional but recommended:** `ggrastr` — if installed, scatterplots with more than 50,000 points are automatically rasterized to keep PDF file sizes and rendering times manageable. The script degrades gracefully (vector rendering) if `ggrastr` is absent. This package is not part of the pinned `environment_qc.yaml` and may be added separately if desired.

### 3.4 Retired Dependencies

The following packages, previously required by the standalone `rds_to_h5ad.r` script, are **no longer needed** anywhere in this pipeline:

| Retired Package | Reason |
|---|---|
| `SeuratDisk` | Replaced by the native `anndata`/`reticulate` conversion path in Section 6.5. |
| `hdf5r` | Was used only for post-conversion metadata-parity verification against the intermediate `h5Seurat` file; this verification step does not exist in the new conversion path (see Section 6, Best Practices, for the implication of this). |

### 3.5 System-Level Dependencies

* A working **Conda / Miniforge / Mamba** installation to build the environment from `environment_qc.yaml`.
* Sufficient **RAM**: Seurat objects for cohorts with hundreds of thousands of cells can occupy many gigabytes in memory. The pipeline invokes `gc()` aggressively between stages, but very large cohorts should be run on a machine with at least 32–64 GB of RAM. Because H5AD conversion now happens **in the same process**, immediately after QC/doublet processing, peak memory usage per sample may be marginally higher than under the old two-process design — plan headroom accordingly.

---

## 4. Installation & Setup

### 4.1 Clone the Repository

```bash
git clone https://github.com/hassansaei/HKOCA-tools/QC/single-cell-qc-pipeline.git
cd single-cell-qc-pipeline
```

### 4.2 Build the Environment from `environment_qc.yaml`

The pipeline ships a single Conda environment specification covering **both** the R analytical stack and the Python packages required for native H5AD export. This is now the canonical, supported installation method — it replaces the previously separate CRAN/Bioconductor/GitHub install steps for `SeuratDisk`/`hdf5r`.

```bash
conda env create -f environment_qc.yaml -n sc_qc_pipeline
conda activate sc_qc_pipeline
```

<details>
<summary><b>Click to expand the environment_qc.yaml file</b></summary>

```yaml
name: sc_qc_pipeline
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - r-base>=4.2.0
  - anndata
  - umap-learn
  - r-seurat>=4.3.0
  - r-mass
  - r-ggplot2
  - r-patchwork
  - r-cowplot
  - r-scales
  - r-reshape2
  - r-ggridges
  - r-ggrepel
  - r-yaml
  - r-dplyr
  - r-jsonlite
  - r-digest
  - r-reticulate
  - r-anndata
  - bioconductor-singlecellexperiment
  - bioconductor-scdblfinder
```

</details>



### 4.3 Verify the `reticulate` ↔ Python Binding

Because the H5AD conversion stage depends on R correctly locating the Conda-provisioned Python interpreter through `reticulate`, it is strongly recommended to verify this binding once before running the full pipeline:

```r
library(reticulate)
reticulate::py_config()   # Confirm it resolves to the sc_qc_pipeline Conda env's Python
library(anndata)
```

If `py_config()` resolves to an unexpected Python installation (e.g., a system Python rather than the Conda environment's), set the `RETICULATE_PYTHON` environment variable to the correct interpreter path before invoking `Rscript`, or activate the Conda environment in the same shell session used to launch the pipeline.

### 4.4 Configure the Pipeline

`QC_scdbl_Combined.R` is driven by a single DCF configuration file, `qc_config.dcf`, expected by default in the same directory as the script. All stages — doublet detection, QC filtering, and now H5AD conversion — read from this one file.

```dcf
rds_dir: /path/to/raw_rds_input
output_dir: /path/to/pipeline_output
summary_subdir: Summary
filtered_subdir: qc_filtered_rds
doublet_subdir: doublet_filtered_rds
log_file: logs/qc_run.log
rds_pattern: \.rds$
recursive_discovery: FALSE
dbl_default_platform: 10x
dbl_default_chemistry: v3
dbl_batch_col: sample_id
dbl_ambiguous_min_overlap: 0.05
reverse_mode: TRUE

# Skip specific problem datasets from QC entirely
skip_noisy_datasets: TRUE
noisy_datasets: sample_lowqual_1,sample_lowqual_2

# Non-10x platform DBR estimates
platform_dbr.dropseq: 0.049
platform_dbr.pipseq: 0.07

# Per-sample chemistry overrides
sample.sample_A.chemistry: v4

# Option A: point to an external filter CSV
# filters: /path/to/my_qc_filters.csv

# Option B: embed the filter table inline (omit if using --filters or accepting MAD3 defaults)
qc_decisions_table:
 Dataset_Name,Lower_Feature_Method,Upper_Feature_Method,Lower_Count_Method,Upper_Count_Method,Upper_Mito_Method
 sample_A,lower_tri,upper_tri,lower_tri,upper_tri,renyi
 sample_B,manual_200,knee,manual_500,upper_tri,mad3

# ── H5AD conversion (new — replaces the standalone rds_to_h5ad.r script) ────
run_h5ad_conversion: TRUE
h5ad_output_subdir: h5ad_converted
```

> **No More Manual Path Synchronization:** In the previous architecture, the operator had to manually edit `INPUT_DIR`/`OUTPUT_DIR` constants inside `rds_to_h5ad.r` to match the QC pipeline's `filtered_subdir` output. This is no longer necessary — `H5AD_DIR` is automatically derived from `output_dir` and `h5ad_output_subdir`, and the conversion stage automatically targets the correct upstream output directory based on which stage(s) just ran (see Section 5.1.2).

---

## 5. Usage Guide & Configuration Options

### 5.1 Running the Pipeline

```bash
Rscript QC_scdbl_Combined.R [--config PATH] [--stage all|qc|doublet] [OPTIONS]
```

There is no separate script or command to invoke for H5AD conversion; it runs automatically within this same invocation whenever `run_h5ad_conversion: TRUE` is set in the active configuration.

#### 5.1.1 Command-Line Flags

| Option/Flag | Data Type | Required (Yes/No) | Default Value | Description |
|---|---|---|---|---|
| `--config` | Path (string) | No | `<script_dir>/qc_config.dcf` | Path to the DCF configuration file. |
| `--filters` | Path (string) | No | None | Path to an **external CSV file** containing per-dataset QC cutoff methods. Takes precedence over the `qc_decisions_table` inline block in the config. |
| `--stage` | Enum string (`all`, `qc`, `doublet`) | No | `all` | Restricts execution to a single stage, or runs the full pipeline. **Does not gate H5AD conversion** — see the note below. |
| `--rds_dir` | Path (string) | No | Value of `rds_dir` in config | Overrides the directory containing raw input `.rds` files. |
| `--output_dir` | Path (string) | No | Value of `output_dir` in config | Overrides the base directory for all pipeline outputs. |
| `--summary_subdir` | String | No | `Summary` | Overrides the QC summary output subdirectory name. |
| `--filtered_subdir` | String | No | `qc_filtered_rds` | Overrides the QC-filtered RDS output subdirectory name. |
| `--doublet_subdir` | String | No | `doublet_filtered_rds` | Overrides the doublet-detection output subdirectory name. |
| `--log_file` | Path (string) | No | `<output_dir>/logs/qc_pipeline.log` | Overrides the log file location. |
| `--recursive_discovery` | Boolean flag | No | `FALSE` | If present, enables recursive search for `.rds` files under `rds_dir`. |
| `--rds_pattern` | Regex (string) | No | `\.rds$` | Overrides the filename regex used to discover input RDS files. |
| `--force_overwrite` | Boolean flag | No | `FALSE` | Forces reprocessing of samples even if a matching `.meta.json` checkpoint indicates no parameter changes, **and** forces re-conversion of existing `.h5ad` files in the H5AD stage. |
| `--help` | Flag | No | — | Prints the usage banner and exits without running the pipeline. Note: the printed banner does not currently list `run_h5ad_conversion`/`h5ad_output_subdir`; these remain config-file (or ad hoc CLI override) keys as described below. |

> **CLI Overrides Extend to Any Config Key:** The CLI parser merges every recognized `--key value` pair directly into the configuration object, not just the flags enumerated in the `--help` banner. This means `--run_h5ad_conversion TRUE` and `--h5ad_output_subdir <name>` are also valid on the command line even though they are not printed by `--help`.

#### 5.1.2 Key DCF Configuration File Keys

| Configuration Key | Data Type | Required (Yes/No) | Default Value | Description |
|---|---|---|---|---|
| `rds_dir` | Path | Yes (or via `--rds_dir`) | None | Source directory of raw input `.rds` Seurat objects. |
| `output_dir` | Path | Yes (or via `--output_dir`) | None | Root directory for all pipeline outputs. |
| `filters` | Path (string) | No | None | Path to an external CSV file for QC cutoff methods. Equivalent to the `--filters` CLI flag; takes priority over `qc_decisions_table` when set. |
| `qc_decisions_table` | Inline CSV block | No | Auto-generated MAD3 | Per-dataset threshold method assignments embedded directly in the DCF file. If neither `filters` nor `qc_decisions_table` is provided, the pipeline automatically generates a `mad3`-for-all-datasets fallback table from the discovered RDS files. |
| `skip_noisy_datasets` | Boolean | No | `TRUE` | When `TRUE`, any datasets listed in `noisy_datasets` are excluded from the QC filtering stage entirely. |
| `noisy_datasets` | Comma-separated string | No | `""` (empty) | Comma-separated list of dataset base names (without `.rds` extension) to exclude from QC processing when `skip_noisy_datasets` is `TRUE`. |
| `summary_subdir` | String | No | `Summary` | QC plots/CSV subdirectory name. |
| `filtered_subdir` | String | No | `qc_filtered_rds` | Final QC-filtered RDS subdirectory name. |
| `doublet_subdir` | String | No | `doublet_filtered_rds` | Doublet-detection output subdirectory name. |
| **`run_h5ad_conversion`** | Boolean | **No (New)** | `FALSE` | Enables the integrated H5AD Conversion stage (Section 6.5). When `TRUE`, the stage runs automatically after the main stage dispatch, independent of the `--stage` value. |
| **`h5ad_output_subdir`** | String | **No (New)** | `h5ad_converted` | Subdirectory (relative to `output_dir`) where converted `.h5ad` files are written. Resolved path is `H5AD_DIR = file.path(output_dir, h5ad_output_subdir)`. |
| `integrated_summary_subdir` | String | No | `""` (i.e., `output_dir` itself) | Subdirectory for the Section 8 integrated summary; if empty, written directly into `output_dir`. |
| `log_file` | Path | No | `<output_dir>/logs/qc_pipeline.log` | Log file destination. |
| `dbl_batch_col` | String | No | `sample_id` | Metadata column used to define `scDblFinder` batches within a multi-sample object. |
| `dbl_min_count` | Integer | No | `100` | Minimum `nCount_RNA` for the ghost-cell pre-filter prior to doublet detection. |
| `dbl_min_feature` | Integer | No | `50` | Minimum `nFeature_RNA` for the ghost-cell pre-filter. |
| `dbl_umap_dims` | Integer | No | `20` | Number of PCA dimensions computed for the optional UMAP visualization (`1:N`). |
| `dbl_default_platform` | String | No | `10x` | Default sequencing platform assumed when `sc_protocol` metadata is absent. |
| `dbl_default_chemistry` | String (`v2`,`v3`,`v4`,`5p`) | No | `v3` | Default 10x chemistry assumed when not auto-detected. |
| `dbl_min_cells_run` | Integer | No | `25` | Minimum post-ghost-filter cell count required to attempt `scDblFinder`; below this, all cells are marked singlets. |
| `dbl_ambiguous_min_overlap` | Float (0–1) | No | `0.05` | Minimum normalized KDE density fraction defining the singlet/doublet score "ambiguous zone." Set to `1.0` to disable ambiguous classification entirely. |
| `dbl_ambiguous_kde_n` | Integer | No | `512` | Number of grid points used for the KDE overlap computation. |
| `sample.<name>.platform` | String | No | None | Per-sample override of sequencing platform. |
| `sample.<name>.chemistry` | String | No | None | Per-sample override of 10x chemistry. |
| `sample.<name>.dbr` | Float (0–0.3] | No | None | Per-sample manual override of the expected doublet rate. |
| `platform_dbr.<platform>` | Float | No | None | Fixed expected doublet rate for a non-10x platform. The reference `qc_config.dcf` ships with `platform_dbr.dropseq: 0.049` and `platform_dbr.pipseq: 0.07`. |
| `skip_scdblfinder_samples` | Comma-separated string | No | `""` | Sample names (canonicalized) to bypass `scDblFinder` entirely; all their cells are marked singlets. |
| `recursive_discovery` | Boolean | No | `FALSE` | Enables recursive `.rds` file discovery. |
| `rds_pattern` | Regex | No | `\.rds$` | Filename pattern for input file discovery. |
| `force_overwrite` | Boolean | No | `FALSE` | Forces reprocessing regardless of checkpoint state; also forces H5AD re-conversion. |
| `reverse_mode` | Boolean | No | `TRUE` | Controls the Section 8 integrated summary's assumed processing order (`raw → ghost → doublet → QC`) for the funnel chart and donut chart. |

**Input-directory resolution for the H5AD stage:** When `run_h5ad_conversion: TRUE`, the conversion stage automatically selects its input directory based on `RUN_STAGE`:

* `--stage all` or `--stage qc` → converts from `FILTERED_DIR` (the final QC-filtered output).
* `--stage doublet` → converts from `DOUBLET_DIR` (the post-doublet-detection, pre-QC output).

This means H5AD conversion always targets the output of whichever stage(s) were actually executed in that invocation — it does not require `--stage all`.

#### 5.1.3 QC Filter Criteria: Three-Tier Resolution

The pipeline resolves per-dataset QC cutoff methods using the following **priority order**, stopping at the first source that provides valid data:

1. **External CSV file** — supplied via `--filters PATH` on the CLI or the `filters` key in the DCF config. This is the recommended approach for large cohorts where the filter table is managed separately from pipeline configuration.
2. **Inline DCF block** — the `qc_decisions_table` key embedded directly in `qc_config.dcf`. Convenient for smaller cohorts where all configuration lives in one file.
3. **Auto-generated MAD3 fallback** — when neither of the above is provided, the pipeline discovers all input RDS files and automatically constructs a decisions table applying `mad3` to all five threshold dimensions for every discovered dataset. A `WARN`-level log message is emitted to alert the operator that defaults are in use.

Regardless of which source is used, the decisions table must conform to the following schema — a header row and one row per dataset with these **required** columns:

| Column | Description |
|---|---|
| `Dataset_Name` | Must match the canonicalized base filename of the corresponding `.rds` file. |
| `Lower_Feature_Method` | Threshold algorithm for the lower `nFeature_RNA` bound. |
| `Upper_Feature_Method` | Threshold algorithm for the upper `nFeature_RNA` bound. |
| `Lower_Count_Method` | Threshold algorithm for the lower `nCount_RNA` bound. |
| `Upper_Count_Method` | Threshold algorithm for the upper `nCount_RNA` bound. |
| `Upper_Mito_Method` | Threshold algorithm for the upper `percent.mito` bound. |

**Available threshold method strings** (used in any of the five method columns):

| Method String | Algorithm |
|---|---|
| `lower_tri` | Triangle-threshold method, searching from the histogram peak toward the lower endpoint. |
| `upper_tri` | Triangle-threshold method, searching from the histogram peak toward the upper endpoint. |
| `renyi` | Rényi entropy maximization split point across a 256-bin histogram. |
| `knee` | Kneedle algorithm applied to the log-log rank-ordered distribution. |
| `mad3` | Median ± 3 × Median Absolute Deviation (MAD). |
| `none` | No bound applied (`-Inf`/`Inf` depending on direction). |
| `manual_<value>` | A fixed, analyst-specified numeric cutoff (e.g., `manual_500`). |

> **⚠️ STRICT FORMAT REQUIREMENTS — READ BEFORE SUBMITTING A FILTER CSV**
>
> Malformed filter tables (whether supplied via `--filters` or `qc_decisions_table`) are a common cause of pipeline crashes. The parser enforces the following rules with **no tolerance for deviation**:
>
> 1. **The first column header must be exactly `Dataset_Name`** — case-sensitive, with no leading/trailing spaces, no prefixes, and no alternate spellings (e.g., `dataset_name`, `Dataset Name`, or `Sample_Name` will all fail).
> 2. **Only the following six literal strings are accepted** in the `Lower_Feature_Method`, `Upper_Feature_Method`, `Lower_Count_Method`, `Upper_Count_Method`, and `Upper_Mito_Method` columns: **`lower_tri`, `upper_tri`, `renyi`, `knee`, `mad3`, `none`**. Any other bare word (typos, alternate casing, or unsupported algorithm names) will not be recognized.
> 3. **Manual cutoffs require the exact `manual_<number>` syntax** — e.g., `manual_5000` or `manual_15`. A **raw number alone** (e.g., just `5000` in the cell) is **not valid** and will cause the pipeline to fail; the `manual_` prefix is mandatory and must directly precede the numeric value with no space (`manual_ 5000` is also invalid).
> 4. **You can find a CSV template in the directory**

#### 5.1.4 Stage 4 (Optional): H5AD Conversion

In Stage 4, the integrated conversion stage (Section 6.5) is intentionally lightweight compared to the retired `SeuratDisk`-based approach. For each `.rds` file in the resolved input directory, it:

1. Loads the Seurat object and, if necessary, casts a Seurat v5 `Assay5` to a legacy `Assay`.
2. Extracts the **raw counts matrix only** (`GetAssayData(..., slot = "counts")`), transposed to Cells × Genes as AnnData expects.
3. Extracts `meta.data` as-is and assigns it to `obs`.
4. Constructs the `AnnData` object (`X = counts_t, obs = meta_data`) and writes it via `$write_h5ad()`.

**Not currently carried over:** normalized/scaled data layers, dimensionality reductions (PCA/UMAP embeddings computed internally during doublet detection), gene-level (`var`) metadata beyond default AnnData indexing, and any `uns`-level unstructured metadata. There is also **no automated post-conversion parity check** (the old script's `hdf5r`-based column-count comparison has no equivalent here). See Section 6 for guidance on when this lightweight approach is sufficient versus when manual verification is warranted.

> **Naming Note:** This document refers to H5AD Conversion as **Stage 4**, reflecting its position as the fourth and fully optional stage in the end-to-end pipeline flow (Doublet → QC → Integrated Summary → H5AD Conversion, gated by `run_h5ad_conversion`). The pipeline's own internal log output for this stage still prints the legacy label `STAGE 3: H5AD CONVERSION` (a holdover from the code's internal section numbering, Section 6.5). Both labels refer to the same stage — if you are correlating this document against raw log output, expect to see `STAGE 3` there.

### 5.2 Concrete Usage Examples

**Example 1 — Full pipeline (doublet → QC), default configuration, no H5AD conversion:**

```bash
Rscript QC_scdbl_Combined.R --config ./qc_config.dcf
```

**Example 2 — Run only the doublet detection stage, forcing a clean reprocess:**

```bash
Rscript QC_scdbl_Combined.R --config ./qc_config.dcf --stage doublet --force_overwrite
```

**Example 3 — Run only the QC filtering stage on data already passed through doublet detection, using a non-default RDS directory and recursive search:**

```bash
Rscript QC_scdbl_Combined.R \
  --config ./qc_config.dcf \
  --stage qc \
  --rds_dir /data/raw_rds/batch_2024 \
  --recursive_discovery \
  --rds_pattern "\\.rds$"
```

**Example 4 — Supply QC cutoffs from an external CSV file, bypassing the inline config table:**

```bash
Rscript QC_scdbl_Combined.R --config ./qc_config.dcf --filters /data/cohort_qc_filters.csv
```

**Example 5 — Run QC with no filter table supplied, relying on the automatic MAD3 fallback for all datasets:**

```bash
Rscript QC_scdbl_Combined.R --config ./qc_config.dcf --stage qc --output_dir /scratch/experiment_42
```

> The pipeline will emit a `WARN`-level log message confirming that auto-generated `mad3` defaults are in use. Review the `qc_summary_detailed.csv` output and `median_inflation_ratio` column to verify the defaults are appropriate before treating results as final.

**Example 6 — Full pipeline with integrated H5AD export enabled via the config file (`run_h5ad_conversion: TRUE`):**

```bash
Rscript QC_scdbl_Combined.R --config ./qc_config.dcf
```

> No additional flag is required if `run_h5ad_conversion: TRUE` is already set in `qc_config.dcf`. The `.h5ad` files will be written to `<output_dir>/<h5ad_output_subdir>/` (default: `h5ad_converted`) automatically after the QC stage completes.

**Example 7 — Enable H5AD export ad hoc from the command line without editing the config file:**

```bash
Rscript QC_scdbl_Combined.R --config ./qc_config.dcf --run_h5ad_conversion TRUE --h5ad_output_subdir h5ad_export_v2
```

**Example 8 — Convert only the doublet-detection output to h5ad (skip QC filtering for this run) by combining `--stage doublet` with H5AD conversion enabled in the config:**

```bash
Rscript QC_scdbl_Combined.R --config ./qc_config.dcf --stage doublet --run_h5ad_conversion TRUE
```

**Example 9 — Override the output base directory for a single experimental run without modifying the config file:**

```bash
Rscript QC_scdbl_Combined.R --config ./qc_config.dcf --output_dir /scratch/experiment_42
```

**Example 10 — Display the full usage banner:**

```bash
Rscript QC_scdbl_Combined.R --help
```

---

## 6. Best Practices & Preferred Usage

1. **Treat `run_h5ad_conversion` as an explicit, opt-in switch.** It defaults to `FALSE`. Enable it only once you are satisfied with the QC and doublet-detection results for a cohort, since the conversion stage runs against whatever RDS files currently exist in the resolved input directory — it will happily convert partially tuned or exploratory results if invoked prematurely.

2. **Understand that the new conversion path is intentionally minimal.** Unlike the retired `SeuratDisk`-based script, the integrated H5AD stage transfers **only the raw counts matrix and `meta.data`** — it does not carry over normalized data, scaled data, or any PCA/UMAP embeddings, and it performs no automated metadata-parity verification after writing. If your downstream Python workflow requires normalized layers, embeddings, or gene-level annotations, perform those steps in Python post-conversion (e.g., via `scanpy.pp.normalize_total`) rather than expecting them to be present in the exported `.h5ad`.

3. **Verify the `reticulate` Python binding before a production run**, especially on shared or multi-user systems where multiple Python installations may be discoverable on the `PATH`. A silent binding to the wrong Python environment will surface as cryptic `anndata`-related errors rather than an obvious configuration error.

4. **Choose the appropriate filter input method for your workflow.** Use an external `--filters` CSV file when QC thresholds are managed by a separate curation process or shared across pipeline runs; use the inline `qc_decisions_table` block for self-contained, single-project configurations; and treat the automatic MAD3 fallback as an **exploratory convenience only** — never as a production-grade default for publication-quality data.

5. **Use the parameter-hash checkpointing system (`.meta.json` files) to your advantage.** Re-running the script with an unchanged configuration will automatically skip already-processed samples. The H5AD conversion stage has its own, simpler skip rule: it skips any `.rds` file whose corresponding `.h5ad` already exists and is non-empty, unless `force_overwrite` is set — it does **not** use parameter-hash checkpointing, so changing an unrelated QC parameter will not by itself trigger reconversion of an existing `.h5ad`.

6. **Use `noisy_datasets` and `skip_noisy_datasets` to gate problematic samples explicitly.** Rather than deleting entries from the filter table, add dataset names to the `noisy_datasets` comma-separated list and set `skip_noisy_datasets: TRUE`. Because H5AD conversion reads whatever `.rds` files exist in the resolved output directory, excluding a dataset from QC this way also naturally excludes it from the QC-stage H5AD output.

7. **Set `dbl_ambiguous_min_overlap` deliberately.** The default (`0.05`) introduces a third classification category, `"ambiguous"`, for cells whose doublet score falls in the KDE overlap zone between the singlet and doublet score distributions. Set to `1.0` to force purely binary singlet/doublet calls.

8. **Always supply `sc_protocol` metadata** (e.g., `"10x_3_v3.1"`, `"10x_5"`) in your raw Seurat objects' metadata when possible, and use per-sample `sample.<name>.chemistry`/`sample.<name>.dbr` DCF overrides when metadata is absent or known to be incorrect.

9. **Inspect the `median_inflation_ratio` diagnostic** in the QC summary CSV before trusting automatically derived thresholds, including the MAD3 auto-fallback, over manually tuned ones.

> **Performance Note:** Because H5AD conversion now runs in the **same R process** as QC/doublet detection rather than as a separate script invocation, memory is not fully released between the analytical stages and the conversion stage the way it was when they ran as independent processes. For very large cohorts, consider running with `--stage qc` (or `doublet`) and H5AD conversion **disabled**, then invoking the pipeline a second time with `--stage qc --run_h5ad_conversion TRUE` and `force_overwrite: FALSE` to let conversion run as its own, freshly started process against already-computed RDS outputs.

> **Data Fidelity Warning:** Because the integrated conversion stage performs no automated metadata-parity verification (unlike the retired script's `hdf5r`-based column comparison), it is good practice to spot-check at least one converted `.h5ad` file per cohort in Python (`anndata.read_h5ad(...).obs.columns`) against the source Seurat object's `meta.data` column names, particularly after changing which metadata columns are attached upstream.

> **Data Integrity Warning:** The doublet detection stage applies a **ghost-cell filter** (`nCount_RNA > dbl_min_count & nFeature_RNA > dbl_min_feature`) immediately before running `scDblFinder`. This filtering is **separate from and precedes** the main QC thresholds. Ensure `dbl_min_count`/`dbl_min_feature` are set conservatively (low) if your QC thresholds are intended to be the sole arbiter of cell inclusion.

> **Security/Provenance Note:** The `.meta.json` checkpoint files embed a hash of the run parameters but do not validate the *content* of the input `.rds` file itself, and — as noted above — the H5AD stage's own skip logic does not consult these checkpoints at all. Use `--force_overwrite` whenever source data provenance is uncertain for either the analytical stages or the conversion stage.

---

## 7. Troubleshooting & Common Errors

| # | Problem | Cause | Solution |
|---|---|---|---|
| 1 | `Missing required packages: anndata, reticulate` (or other packages) at startup | One or more R packages from the `load_packages_safely()` dependency check are not installed in the active R library path — commonly `anndata`/`reticulate` if the environment was built before this revision's `environment_qc.yaml` was adopted. | Rebuild or update the Conda environment from the provided `environment_qc_yml.yaml` (`conda env create -f environment_qc_yml.yaml` or `conda env update -f environment_qc_yml.yaml --prune`), and confirm the correct environment is activated before invoking `Rscript`. |
| 2 | `Error: Unable to find conda binary` / `Error: Python not found` or similarly cryptic errors originating from `anndata`/`reticulate` calls during the H5AD stage | `reticulate` could not locate a usable Python interpreter with the `anndata` package installed — often because the Conda environment was not activated in the shell used to launch `Rscript`, or a competing Python installation is earlier on the `PATH`. | Run `reticulate::py_config()` interactively in the same environment to confirm which Python is being used. If incorrect, either activate the `sc_qc_pipeline` Conda environment before running the pipeline, or explicitly set `Sys.setenv(RETICULATE_PYTHON = "/path/to/conda/envs/sc_qc_pipeline/bin/python")` before the script's package-loading step. |
| 3 | `Filter CSV file not found: <path>` | The path supplied to `--filters` (or the `filters` key in the DCF config) does not exist at the specified location. | Verify the path is absolute or correctly relative to the working directory from which `Rscript` is invoked. |
| 4 | Pipeline emits `WARN: No filter criteria provided. Auto-generating global 'mad3' default criteria` unexpectedly | Neither `--filters`, the `filters` config key, nor `qc_decisions_table` was set, triggering the automatic MAD3 fallback. | Confirm the DCF file contains either a `filters` path or a populated `qc_decisions_table` block with properly indented continuation lines (each continuation line must begin with at least one leading space). |
| 5 | `[ZERO CELLS] <sample>: all N cells removed by QC thresholds` (logged as an error; sample is skipped) | The computed lower/upper thresholds for nFeature, nCount, or percent.mito are too aggressive for that specific dataset's distribution. | Review the logged threshold values against the logged data ranges. Override the problematic row in the filter CSV or `qc_decisions_table` with a `manual_<value>` method, or switch to a different automatic method. If the dataset is known to be low quality, add it to `noisy_datasets` instead of removing it from the filter table. |
| 6 | `[FAILED H5AD] <sample>: <error message>` (logged during Stage 3; that sample's `.h5ad` is not produced, pipeline continues) | The `AnnData()` constructor raised an error — commonly because `meta.data` contains a column type that cannot be coerced cleanly into a Python `pandas.DataFrame` (e.g., nested list-columns, unusual factor encodings), or because the counts matrix is unexpectedly empty after the v5→v3 assay cast. | Inspect the sample's `meta.data` for non-scalar or exotic column types and coerce them to simple character/numeric vectors upstream. Confirm the assay actually contains a non-empty `counts` slot by checking the preceding pipeline stage's cell/feature counts in the log. |
| 7 | `.h5ad` output directory (`h5ad_output_subdir`) is empty after a run despite `run_h5ad_conversion: TRUE` being set | The H5AD stage's input-directory resolution rule selected an upstream directory that is itself empty — e.g., `--stage doublet` was used but the doublet stage produced zero singlet files for every sample, or `--stage qc` was used before any QC-filtered RDS files existed. | Confirm the resolved input directory (`FILTERED_DIR` for `--stage all`/`qc`, `DOUBLET_DIR` for `--stage doublet`) actually contains `.rds` files before re-running with H5AD conversion enabled; check the "STAGE 3: H5AD CONVERSION" log block for the exact input path it used. |
| 8 | `[H5AD SKIPPED] <sample>.h5ad already exists.` for every sample even after upstream data changed | The H5AD stage's skip rule only checks for the **existence and non-zero size** of the target `.h5ad` file — it does not use the parameter-hash checkpointing (`.meta.json`) system that the analytical stages use, so it cannot detect that upstream RDS content changed. | Pass `--force_overwrite` (or set `force_overwrite: TRUE` in the config) to force re-conversion, or manually delete the stale `.h5ad` files in `h5ad_output_subdir` before re-running. |
| 9 | `scDblFinder failed (...) — marking all N cells as singlets.` (warning, pipeline continues) | `scDblFinder` raised an internal error, commonly due to a near-singular size-factor estimation on a sample with very low cell counts or extreme batch imbalance under `dbl_batch_col`. | This is a non-fatal, intentional fallback. If the sample is structurally too small for reliable doublet calling, add it to `skip_scdblfinder_samples` to suppress the warning and document the decision explicitly. |
