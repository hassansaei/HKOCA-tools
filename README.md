# Human Kidney Organoid Cell Atlas (NephAtlas)

Modular toolkit for ambient RNA removal, QC/doublet filtering (with
harmonization), annotation, and integration / atlas projection.

## Modules

| CLI | Package | Role |
|---|---|---|
| `hkoca pipeline` | `hkoca.pipeline` | End-to-end A–Z analysis (orchestrator) |
| `hkoca cellbender h5` | `hkoca.cellbender` | CellBender on CellRanger raw H5 |
| `hkoca cellbender mtx` | `hkoca.cellbender` | CellBender on CellRanger MTX dirs |
| `hkoca qc-filter run` | `hkoca.qc_filter` | Harmonize (`--to-rds`) then doublet/QC |
| `hkoca qc-filter harmonize` | `hkoca.qc_filter.harmonize` | Harmonization only |
| `hkoca qc-filter qc` | `hkoca.qc_filter` | QC R script only (existing RDS) |
| `hkoca annotation` | `hkoca.annotation` | Cell-type annotation |
| `hkoca integration` | `hkoca.integration` | Batch integration / atlas projection |

```
hkoca/
  pipeline/                 # A–Z orchestrator
  cellbender/               # CellBender only
  qc_filter/
    harmonize/              # gene-space harmonization (before QC)
    r/                      # doublet detection + QC (R / scDblFinder)
  annotation/
  integration/
conda/                      # Stage-specific conda / Apptainer environments
```

**QC-filter always harmonizes before doublet/QC** — harmonization is not a
separate public module.

## Install

```bash
conda env create -f conda/environment_harmonize.yaml
conda activate hkoca_harmonize
pip install -e .
```

## CLI

```bash
hkoca --version
hkoca pipeline --list-stages
hkoca cellbender h5 --help
hkoca cellbender mtx --help
hkoca qc-filter
hkoca qc-filter run --help
hkoca qc-filter harmonize --help
hkoca qc-filter qc --print-script
hkoca annotation
hkoca integration
```

### QC-filter order

```bash
# 1 + 2 in one command (harmonize always first, forces --to-rds)
hkoca qc-filter run \
  --csv meta.csv --gtf genes.gtf --output results \
  --qc-output results/qc_filter --stage all

# Or stepwise:
hkoca qc-filter harmonize --csv meta.csv --gtf genes.gtf --output results --to-rds
hkoca qc-filter qc --rds-dir results --output-dir results/qc_filter --stage all
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
