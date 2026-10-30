# HKOCA-tools Nextflow

Nextflow DSL2 wrapper around the existing HKOCA-tools (`hkoca`) pipeline. Analysis code, thresholds, and configs are unchanged; this layer schedules the same CLIs and writes standard Nextflow reports.

Place these files at the root of an [HKOCA-tools](https://github.com/haythem03/HKOCA-tools) checkout (next to `hkoca/` and `conda/`).

## Stages

| Process | Calls |
|---|---|
| `CELLBENDER` | `hkoca cellbender` (optional) |
| `PREPARE_HARMONIZE_META` | rewrite `data_path` after CellBender |
| `HARMONIZATION` | `hkoca.qc_filter.harmonize --to-rds` |
| `QC` | QC / scDblFinder |
| `ANNOTATION` | `hkoca annotation run` |
| `INTEGRATION` | `hkoca integration prep` + `run` |
| `PROJECTION` | `hkoca projection map` (optional, GPU) |

## Requirements

- Compute node (not a login node)
- Nextflow >= 23.04, Java 17+
- Conda envs from `conda/*.yaml`, with `hkoca` installed in the harmonize env (`pip install -e .`)
- GRCh38 GTF
- For projection: GPU, atlas `.h5ad`, and `scPoli_Reference_Model`

## Run

```bash
export HKOCA_ROOT="$PWD"
export NXF_HOME="${HKOCA_ROOT}/.nextflow_home"
export TMPDIR="${HKOCA_ROOT}/tmp"
mkdir -p "$NXF_HOME" "$TMPDIR"

nextflow run . -profile conda \
  --sample_info path/to/sample_info.csv \
  --gtf path/to/Homo_sapiens.GRCh38.104.gtf \
  --outdir path/to/results \
  --skip_cellbender

nextflow run . -profile conda,gpu \
  --sample_info path/to/sample_info.csv \
  --gtf path/to/Homo_sapiens.GRCh38.104.gtf \
  --outdir path/to/results \
  --skip_cellbender \
  --run_projection \
  --atlas path/to/atlas.h5ad \
  --model_dir path/to/scPoli_Reference_Model
```

Resume with `-resume`. Details: `docs/nextflow.md`, `docs/DEPENDENCY_MAP.md`.

## Layout

```
main.nf
nextflow.config
workflows/hkoca.nf
modules/          # process wrappers + prepare_harmonize_csv.py
conda/            # env YAMLs used by -profile conda
docs/
tests/
use_cases/GSM9712682/
```

## Example: GSM9712682

`use_cases/GSM9712682/` is a public GEO sample run through this workflow. It keeps plots and summary tables only (no `.h5ad` / `.rds`). See that folder’s README for dataset metadata and what each stage produced.
