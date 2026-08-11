# Human Kidney Organoid Cell Atlas (NephAtlas)

Toolkit to process **new** single-cell datasets with the **same analysis pipeline** used to build the Human Kidney Organoid Cell Atlas, then place those datasets onto the atlas for comparison with the reference.

Typical path: CellRanger outputs -> ambient RNA removal (CellBender) -> gene-space harmonization and QC/doublet filtering -> annotation -> integration /
projection onto the atlas.

## Modules

| CLI | Package | Role |
|---|---|---|
| `hkoca pipeline` | `hkoca.pipeline` | End-to-end A-Z analysis (orchestrator) |
| `hkoca cellbender h5` | `hkoca.cellbender` | CellBender on CellRanger raw H5 |
| `hkoca cellbender mtx` | `hkoca.cellbender` | CellBender on CellRanger MTX dirs |
| `hkoca qc-filter` | `hkoca.qc_filter` | Harmonize (`--to-rds`) then doublet/QC |
| `hkoca annotation` | `hkoca.annotation` | Cell-type annotation |
| `hkoca integration` | `hkoca.integration` | Batch integration / atlas projection |

```
hkoca/
  pipeline/                 # A-Z orchestrator
  cellbender/               # CellBender only
  qc_filter/
    harmonize/              # internal: gene-space harmonization (step 1)
    r/                      # internal: doublet detection + QC (step 2)
  annotation/
  integration/
conda/                      # Stage-specific conda / Apptainer environments
```

`hkoca qc-filter` is a single command: harmonization always runs first, then
the QC R script. There are no separate harmonize/qc CLIs.

## Install

```bash
conda env create -f conda/environment_harmonize.yaml
conda env create -f conda/environment_qc.yaml
conda activate hkoca_harmonize
pip install -e .
```

`hkoca qc-filter` runs harmonization in the active env, then automatically
invokes the QC R script with `sc_qc_pipeline` (no mid-run `conda activate`).
Override with `HKOCA_QC_ENV` or `HKOCA_RSCRIPT` if needed.

## CLI

```bash
hkoca --version
hkoca pipeline --list-stages
hkoca pipeline --print-template
hkoca cellbender h5 --help
```

### Full pipeline

```bash
# CellBender -> harmonize/QC -> annotation/integration (stubs skipped with a log)
hkoca pipeline \
  --csv sample_info.csv \
  --gtf genes.gtf \
  --output /data/out

# Skip CellBender when inputs are already filtered or you use raw counts directly
hkoca pipeline \
  --csv sample_info.csv \
  --gtf genes.gtf \
  --output /data/out \
  --skip-cellbender

# Resume is the default everywhere (CellBender, harmonize, QC).
# Finished outputs under --output are skipped. Force a full re-run:
hkoca pipeline \
  --csv sample_info.csv \
  --gtf genes.gtf \
  --output /data/out \
  --force

# Run only through QC
hkoca pipeline \
  --csv sample_info.csv \
  --gtf genes.gtf \
  --output /data/out \
  --to-stage qc_filter
```

The **sample_info CSV** uses the same columns as harmonize metadata, plus optional
pipeline columns: `sample_dir`, `run_cellbender`, `file_prefix`, `skip`.
Print the packaged template path with `hkoca pipeline --print-template`.

Outputs are written under `--output`:

- CellBender: `{sample_dir}/{sample_id}_filtered.h5` (when enabled)
- Harmonize: `{output}/{study}/raw|harmonized|rds/`
- QC: `{output}/qc_filter/` (or `--qc-output`)
- Annotation / integration: reserved subdirs (stages not implemented yet)

### QC-filter

```bash
# Always: harmonize (--to-rds) then doublet detection + QC
hkoca qc-filter \
  --csv meta.csv --gtf genes.gtf --output results \
  --qc-output results/qc_filter --run all
```

### CellBender quick start

```bash
conda activate hkoca_cellbender
pip install -e .

# Edit samples + tunables in hkoca/cellbender/cellbender.config
# CLI flags override config values
hkoca cellbender h5 \
  --samples-dir /data/samples \
  --samples sampleA,sampleB \
  --epochs 100 \
  --total-droplets-included 35000 \
  --learning-rate 0.00005 \
  --cpu-threads 24

# Same parameters, MTX input (directory with matrix/features/barcodes)
hkoca cellbender mtx \
  --samples-dir /data/samples \
  --samples sampleA,sampleB \
  --dry-run
```
