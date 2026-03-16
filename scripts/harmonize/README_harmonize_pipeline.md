# harmonize_pipeline.py

A standalone Python script that processes scRNA-seq datasets from a metadata CSV, harmonizes all samples to a common 24,100-gene reference (GRCh38.104 protein-coding + lncRNA), and outputs ready-to-use `.h5ad` files.

---

## What it does

For each study defined in your metadata CSV, the pipeline:

1. Reads the GTF reference file and builds the 24,100-gene target set
2. Loads each sample — auto-detecting H5, H5AD, or MTX triplet format
3. Attaches all metadata columns (Age, protocol, condition, etc.) per cell
4. Concatenates all samples belonging to the same study
5. Saves a raw `.h5ad` (all original genes)
6. Harmonizes the gene space: drops non-coding genes, zero-fills missing genes, casts to float32 sparse
7. Saves a `_harmonized.h5ad` (exactly 24,100 genes, float32 CSR)
8. Optionally generates publication-ready summary plots

---

## Requirements

```bash
pip install scanpy anndata pandas numpy scipy h5py seaborn matplotlib
```

---

## Quick start

```bash
# 1. Edit the USER CONFIG section at the top of the script
#    (or set environment variables — see Configuration below)

# 2. Run
python harmonize_pipeline.py
```

---

## Configuration

There are three ways to set the four required paths, in order of priority:

### 1. Command-line arguments (highest priority)
```bash
python harmonize_pipeline.py \
  --csv    /path/to/datasets_metadata.csv \
  --gtf    /path/to/Homo_sapiens.GRCh38.104.gtf \
  --data   /path/to/raw_data_root \
  --output /path/to/results
```

### 2. Environment variables
```bash
export DATA_ROOT="/data-master/.../data"
export GTF_FILE="/data-master/.../Homo_sapiens.GRCh38.104.gtf"
export METADATA_CSV="/data-master/.../datasets_metadata.csv"
export OUTPUT_ROOT="/data-master/.../CSV_driven_results"

python harmonize_pipeline.py
```

Add those `export` lines to your `~/.bashrc` to make them permanent:
```bash
echo 'export GTF_FILE="/data-master/.../Homo_sapiens.GRCh38.104.gtf"' >> ~/.bashrc
source ~/.bashrc
```

### 3. USER CONFIG block (default fallback)

Edit the four variables at the top of the script directly:

```python
DATA_ROOT    = "/data-master/pure-workspace/labss/hmami/new_data/data"
GTF_FILE     = "/data-master/pure-workspace/labss/hmami/new_data/meta/gtf/Homo_sapiens.GRCh38.104.gtf"
METADATA_CSV = "/data-master/pure-workspace/labss/hmami/new_data/notebooks/data_loaders/datasets_metadata.csv"
OUTPUT_ROOT  = "/data-master/pure-workspace/labss/hmami/new_data/notebooks/data_loaders/CSV_driven_results"
```

---

## All run modes

```bash
# Run the full pipeline (default)
python harmonize_pipeline.py

# Run pipeline + generate summary plots afterwards
python harmonize_pipeline.py --summary

# Only generate summary plots from already-processed files (skip pipeline)
python harmonize_pipeline.py --summary-only

# Override any path on the fly
python harmonize_pipeline.py --csv my_new_datasets.csv

# Full override example
python harmonize_pipeline.py \
  --csv    datasets_metadata.csv \
  --gtf    Homo_sapiens.GRCh38.104.gtf \
  --data   /data/raw \
  --output /data/results \
  --summary
```

---

## Output structure

For each study in the CSV, two files are written under `OUTPUT_ROOT/{study}/`:

```
CSV_driven_results/
└── d10_1038_s42003_024_07069_6/
    ├── d10_1038_s42003_024_07069_6.h5ad                ← raw (all original genes)
    └── d10_1038_s42003_024_07069_6_harmonized.h5ad     ← harmonized (24,100 genes, float32)
```

The harmonized file contains:
- **Rows (obs)**: all cells with `sample_id`, `study`, `Age`, `condition`, and all other metadata columns from the CSV attached per cell
- **Columns (var)**: exactly the 24,100 protein-coding + lncRNA genes from `Homo_sapiens.GRCh38.104.gtf`
- **Matrix**: sparse float32 CSR — missing genes are zero-filled, non-coding genes are dropped

---

## Summary plots (--summary)

When `--summary` or `--summary-only` is passed, the script scans all `.h5ad` files under `OUTPUT_ROOT` and writes to `OUTPUT_ROOT/reports/atlas_summary/`:

| File | Description |
|---|---|
| `dataset_summary.csv` | One row per h5ad file with cell count, gene count, protocol, age, technology |
| `cells_per_dataset.png/pdf` | Horizontal bar chart — cells per study |
| `cells_by_protocol.png/pdf` | Donut chart — cell contribution by differentiation protocol |
| `cells_by_age.png/pdf` | Bar chart — cell contribution by organoid timepoint |

Plots are saved to disk only — safe to run on a remote server without a display.

---

## CSV format

The metadata CSV must be comma-separated with at least these five columns:

| Column | Description |
|---|---|
| `data_path` | Exact path to the `.h5`/`.h5ad` file, or the MTX folder |
| `file_prefix` | MTX filename prefix (e.g. `GSM8194405_d15_kidney_organoids_`). Leave empty for H5/H5AD or unprefixed MTX. |
| `sample_id` | Unique sample identifier (GEO accession preferred) |
| `output_dir` | Where to write outputs for this study (used only if `OUTPUT_ROOT` is not set) |
| `study` | Grouping key — all rows with the same value are concatenated |

Optional columns attached as per-cell metadata: `source`, `diff_protocol`, `sc_protocol`, `sequencing`, `genome_build`, `Age`, `type`, `disease`, `condition`.

Add a `skip` column and set it to `True` on any row to exclude that sample without deleting it from the CSV.

**One row = one sample. Each timepoint in a time-course must be its own row** with the appropriate `Age` value.

See `README_CSV_pipeline.md` for the full CSV guide with examples.

---

## Supported file formats

| Format | `data_path` | `file_prefix` |
|---|---|---|
| 10x H5 (`.h5`) | Full path to the `.h5` file | empty |
| H5AD (`.h5ad`) | Full path to the `.h5ad` file | empty |
| MTX triplet with prefix | Path to the folder containing the triplet | GSM prefix string |
| MTX triplet, no prefix | Path to the triplet folder | empty |

The script auto-detects the format from `data_path`. For MTX samples with a GEO `sample_id` (starting with `GSM`), a failed read automatically triggers a re-download from the GEO FTP before retrying.

---

## Exit codes

| Code | Meaning |
|---|---|
| `0` | All studies processed successfully |
| `1` | One or more studies or samples failed (failures listed in terminal output) |

---

## Common errors

| Error | Cause | Fix |
|---|---|---|
| `FileNotFoundError: GTF file not found` | `GTF_FILE` path is wrong or not set | Set `GTF_FILE` in USER CONFIG or as env var |
| `FileNotFoundError: CSV not found` | `METADATA_CSV` path is wrong | Set `METADATA_CSV` in USER CONFIG or pass `--csv` |
| `CSV missing required columns` | CSV is semicolon-separated or has whitespace in headers | Re-save as comma-separated; the script auto-detects `;` and warns |
| `data_path not found` | Path in CSV does not exist on the server | Check the path is absolute or relative to `DATA_ROOT` |
| `No features/genes file for prefix` | `file_prefix` is wrong or files are missing | Verify the prefix matches the actual filenames in the MTX folder |
| `[Errno 28] No space left on device` | Output partition is full | Point `OUTPUT_ROOT` to a partition with enough space |
