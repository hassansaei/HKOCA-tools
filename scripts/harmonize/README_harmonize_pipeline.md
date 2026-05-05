# harmonize_pipeline.py

A standalone Python script that processes scRNA-seq datasets from a metadata CSV, harmonizes all samples to a common 24,100-gene reference (GRCh38.104 protein-coding + lncRNA), and outputs ready-to-use `.h5ad` files. It also includes an optional, fully automated bridge to convert harmonized data directly into Seurat `.rds` objects.



## Pipeline Overview


![Pipeline Overview](./harmonize_pipeline_overview.svg)



## What it does

For each study defined in your metadata CSV, the pipeline:

1. **Reads the GTF reference file** and builds the 24,100-gene target set
2. **Loads each sample** — auto-detecting H5, H5AD, or MTX triplet format
3. **Attaches standardized metadata** columns (Age, protocol, species, tissue, etc.) per cell
4. **Concatenates all samples** belonging to the same study
5. **Saves a raw `.h5ad`** (all original genes)
6. Harmonizes the gene space: Drops non-coding genes, zero-fills missing genes, and standardizes to a float32 sparse matrix for memory efficiency.
7. Preserves Transgenes: Integrates user-specified transgene/reporter names (EGFP, TTC21B, etc..) into the final feature set by merging them into the reference gene set before processing.
8. Automated GEO Recovery: If an MTX sample with a GSM ID fails to load, the script automatically attempts to re-download the data from the GEO FTP servers.
9. **Saves a `_harmonized.h5ad`** (24,100 reference genes + user transgenes)
10. **(Optional)** Safely bridges to R via `rpy2` to build and save a native Seurat `.rds` object
11. **(Optional)** Generates publication-ready summary plots

---

## Installation & Environment Setup (Recommended)

Because this pipeline bridges Python (Scanpy) and R (Seurat), managing dependencies can be tricky. We **highly recommend** creating a dedicated Conda environment to keep everything isolated and stable.

### Step 1: Create and activate a new Conda environment

```bash
conda create -n harmonize_env python=3.11 -y
conda activate harmonize_env
```

### Step 2: Install R and Seurat dependencies (For the --to-rds bridge)

> **Note:** Installing R via Conda ensures rpy2 can find the correct system libraries.

```bash
conda install -c conda-forge -c bioconda r-base r-seurat bioconductor-singlecellexperiment r-seuratobject -y
```

### Step 3: Install all Python dependencies

```bash
pip install scanpy anndata pandas numpy scipy h5py seaborn matplotlib 'rpy2==3.6.3' 'anndata2ri==2.0'
```
### Step4: Install all required R Bioconductor packages for anndata2ri

```bash
R -e 'if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("S4Vectors","SingleCellExperiment","SummarizedExperiment","DelayedArray","Matrix"), update=FALSE, ask=FALSE); if (!require("Seurat", quietly=TRUE)) install.packages("Seurat", repos="https://cloud.r-project.org")'
```

> **Note:** The script uses "lazy loading." If you absolutely cannot install R, you can skip Step 2 and just install the Python packages. The script will safely skip the .rds conversion step without crashing.

---

## Quick Start

### Create a datasets_metadata.csv manifest

See `config/CSV_metadata_template.csv` for the template.

### Configure your GTF path

Set it in `harmonize.config`.

### Run the pipeline

```bash
# Basic run (Generates raw and harmonized .h5ad)
python scripts/harmonize/harmonize.py --csv config/datasets_metadata.csv --output ./results/

# Full run with Seurat RDS conversion
python scripts/harmonize/harmonize.py --csv config/datasets_metadata.csv --output ./results/ --to-rds

# Run pipeline and generate summary plots
python scripts/harmonize/harmonize.py --csv config/datasets_metadata.csv --output ./results/ --summary

# Run with custom transgenes (e.g., GFP and MCHERRY)
python scripts/harmonize/harmonize.py --csv config/datasets_metadata.csv --transgenes EGFP,MCHERRY
```

---

## Configuration & Paths

The pipeline uses relative paths for maximum portability. You can define your environment using a combination of the CLI, the config file, and the working directory flag.

### 1. Command-line arguments (Highest Priority)

```bash
python harmonize.py \
  --csv    ./config/datasets_metadata.csv \
  --gtf    ./data/meta/gtf/Homo_sapiens.GRCh38.104.gtf \
  --output ./results \
  --transgenes EGFP,MCHERRY \
  --to-rds
```

### 2. The Configuration File (harmonize.config)

Instead of typing the GTF path every time, you can set it in a configuration file in the same directory as the script:

```ini
[paths]
gtf_file = ./data/meta/gtf/Homo_sapiens.GRCh38.104.gtf

[transgenes]
# Comma-separated list of transgene/reporter names to preserve
names = EGFP,MCHERRY,TTC21B,G7

[summary]
figure_dpi = 300
figure_extensions = png, pdf
report_subdir = reports/atlas_summary
```

### 3. The Working Directory Flag (-w)

If you are running the script from a different folder, use `-w` to tell the script where to base its relative paths (`./`) from:

```bash
python harmonize.py -w /path/to/project/root --csv config/datasets_metadata.csv
```

---

## Output Structure

The pipeline automatically creates a clean, standardized folder structure for every study to keep your raw, harmonized, and R-objects strictly separated.

For each study in the CSV, files are written under `OUTPUT_ROOT/{study}/`:

```
results/
└── study_name/
    ├── raw/
    │   └── study_name.h5ad                   ← raw (all original genes)
    ├── harmonized/
    │   └── study_name_harmonized.h5ad        ← harmonized (24,100 genes, float64)
    └── rds/
        └── study_name_harmonized.rds         ← Seurat object (Only if --to-rds is passed)
```

### The harmonized file guarantees:

-**Consistent Metadata**: Mandatory biological columns (Age, protocol, etc.) are strictly enforced, with missing values filled as "Unknown".
-**Flexible Metadata**: Any additional columns in your CSV (such as `disease` or `condition`) are automatically detected and added to the `.obs` slot.
**Consistent Genes**: Exactly the 24,100 protein-coding + lncRNA genes from the GRCh38.104 GTF.
- **Sparse Matrix**: The final output is stored as a `float32` CSR matrix to minimize file size while maintaining precision for downstream analysis.

---

## Summary Plots (--summary)

When `--summary` or `--summary-only` is passed, the script scans all `*_harmonized.h5ad` files under the output directory and writes a report to `OUTPUT_ROOT/reports/atlas_summary/`:

| File | Description |
|------|-------------|
| `dataset_summary.csv` | One row per h5ad file with cell count, gene count, protocol, age, technology |
| `cells_per_dataset.png/pdf` | Horizontal bar chart — cells per study |

> **Note:** Plots are saved to disk only — safe to run on a remote server without a display.

---

## CSV Format

The metadata CSV must be comma-separated. The script requires only three strict columns to run, but accepts many optional metadata columns.

### Required Columns

| Column | Description |
|--------|-------------|
| `data_path` | Exact relative or absolute path to the `.h5/.h5ad` file, or the MTX folder |
| `sample_id` | Unique sample identifier (GEO accession preferred) |
| `study` | Grouping key — all rows with the same value are concatenated together! |

### Optional System Columns

| Column | Description |
|--------|-------------|
| `file_prefix` | MTX filename prefix (e.g., `GSM8194405_d15_`). Leave empty for H5/H5AD |
| `output_dir` | Where to write outputs for this study (Overridden by the `--output` CLI flag) |
| `skip` | Set to `True` on any row to temporarily exclude that sample |

### Standard Metadata Columns (Attached to cells)

The following columns are required in your metadata CSV:
*  `source` 
*  `diff_protocol` 
*  `sc_protocol` 
*  `sequencing` 
*  `genome_build` 
*  `Age` 
*  `type` 

**Automatically Detected Columns**:
 Additional columns found in your CSV are automatically propagated to the AnnData `.obs`. For example:
* `disease`
* `condition`

---

## Supported File Formats

| Format | `data_path` | `file_prefix` |
| :--- | :--- | :--- |
| **10x H5 (`.h5`)** | Full path to the `.h5` file |  Leave empty  |
| **H5AD (`.h5ad`)** | Full path to the `.h5ad` file |  Leave empty  |
| **MTX Folder** | Path to the directory containing `.mtx.gz` |  The prefix including trailing underscores (e.g., `GSM5112197_KR01_`)  |

### Format Auto-Detection

The script auto-detects the format from `data_path`. For MTX samples with a GEO sample_id (starting with `GSM`), a failed read automatically triggers:
1. A fallback manual MTX loader
2. If that fails, a re-download from the GEO FTP

---

## Logging & Error Handling

This script features dual-level logging:

- **Console:** Prints clean `[INFO]` level updates to your terminal so you can track progress
- **File:** Writes a highly detailed `harmonize_pipeline.log` file in the working directory containing line-specific `[DEBUG]` traces and full crash reports

If an individual sample fails to load, the script will gracefully log the error, skip the sample, and continue concatenating the rest of the study.

---

## Common Errors

| Error | Cause | Fix |
| :--- | :--- | :--- |
| `CSV uses semicolons` | The script detected semicolon delimiters |  The script will attempt to read it, but re-saving as a standard comma-separated CSV is recommended. |
| `allowed_genes is empty` | GTF biotypes don't match script filters | Ensure your GTF uses `protein_coding` and `lncRNA` as the `gene_biotype`. |
| `unsafe for R interpolation` | RDS path contains quotes or newlines |  Ensure the `output_root` or study names do not contain special characters that break R string literals. |

---

## License & Questions

For questions or issues, please refer to the project documentation or submit an issue to the repository.
