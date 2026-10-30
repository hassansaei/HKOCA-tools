# Running HKOCA with Nextflow

Nextflow schedules the same `hkoca` stages and conda envs as the CLI pipeline. Defaults stay in the packaged `hkoca` configs.

## Prerequisites

1. Compute node (not login)
2. Java 17+ and Nextflow >= 23.04
3. Conda envs from `conda/*.yaml`
4. `pip install -e .` in the harmonize env; `export HKOCA_ROOT=$PWD`

```bash
conda env create -f conda/environment_harmonize.yaml
conda activate hkoca_harmonize
pip install -e .
```

Also create stage envs as needed: cellbender, qc, integration, projection.

## Quick start

```bash
export HKOCA_ROOT="$PWD"
export NXF_HOME="$HKOCA_ROOT/.nextflow_home"
export TMPDIR="$HKOCA_ROOT/tmp"
mkdir -p "$TMPDIR" "$NXF_HOME"

nextflow run . -profile conda \
  --sample_info path/to/sample_info.csv \
  --gtf path/to/genes.gtf \
  --outdir path/to/results \
  --skip_cellbender
```

Projection (GPU):

```bash
nextflow run . -profile conda,gpu \
  --sample_info path/to/sample_info.csv \
  --gtf path/to/genes.gtf \
  --outdir path/to/results \
  --run_projection \
  --atlas path/to/atlas.h5ad \
  --model_dir path/to/scPoli_Reference_Model
```

## Layout

```
main.nf
nextflow.config
workflows/hkoca.nf
modules/                 # wrappers + prepare_harmonize_csv.py
docs/
use_cases/
```

## Stage mapping

| Process | Command |
|---|---|
| `CELLBENDER` | `python -m hkoca.cellbender h5\|mtx` |
| `PREPARE_HARMONIZE_META` | path rewrite (same rules as `hkoca pipeline`) |
| `HARMONIZATION` | `python -m hkoca.qc_filter.harmonize --to-rds` |
| `QC` | `QC_scdbl_Combined.R` via qc runner |
| `ANNOTATION` | `python -m hkoca.annotation run` |
| `INTEGRATION` | `hkoca integration prep` then `run` |
| `PROJECTION` | `python -m hkoca.projection map` |

Projection takes `prep/sct_prepared.rds` plus `--atlas` and `--model_dir`, not the Harmony/RPCA/CCA objects.

## Existing conda prefixes

```bash
nextflow run . -profile conda \
  --conda_harmonize /path/to/envs/hkoca_harmonize \
  --conda_qc /path/to/envs/sc_qc_pipeline \
  ...
```

Conda cache defaults to `${projectDir}/.nextflow_conda`.

## Equivalence checks

After running both `hkoca pipeline` and `nextflow run .` on the same input, compare cell/gene counts, QC barcodes, Level labels, and projection summaries. Aim for scientific equivalence, not byte-identical files.

## Resume

```bash
nextflow run . -profile conda -resume --sample_info ... --gtf ... --outdir ...
```
