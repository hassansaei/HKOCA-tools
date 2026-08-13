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
| `hkoca integration` | `hkoca.integration` | Batch integration (Harmony / RPCA / CCA) |
| `hkoca projection` | `hkoca.projection` | Map query onto the HKOCA atlas (scPoli) |

```
hkoca/
  pipeline/                 # A-Z orchestrator
  cellbender/               # CellBender only
  qc_filter/
    harmonize/              # internal: gene-space harmonization (step 1)
    r/                      # internal: doublet detection + QC (step 2)
  annotation/
  integration/
  projection/
conda/                      # Stage-specific conda / Apptainer environments
```

`hkoca qc-filter` is a single command: harmonization always runs first, then
the QC R script. There are no separate harmonize/qc CLIs.

## Install

```bash
conda env create -f conda/environment_harmonize.yaml
conda env create -f conda/environment_qc.yaml
conda env create -f conda/environment_integration.yaml
conda env create -f conda/environment_projection.yaml
conda activate hkoca_harmonize
pip install -e .
# Snapseed (annotation) is installed from https://github.com/hassansaei/snapseed
# via environment_harmonize.yaml. To refresh only Snapseed later:
#   pip install -U git+https://github.com/hassansaei/snapseed.git
# Atlas projection (GPU): conda activate hkoca_projection && pip install -e .
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
# CellBender -> harmonize/QC -> annotation -> integration (projection is separate)
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
  - `{output}/qc_filter/qc_filtered_rds/{study}*_filtered.rds`
  - `{output}/qc_filter/h5ad_converted/{study}*_filtered.h5ad`
- Annotation: `{output}/annotation/annotated_obj/{stem}_annotated.h5ad`
- Integration: `{output}/integration/prep/sct_prepared.rds` and
  `{output}/integration/objects/integrated_{method}.rds`
  Prep diagnostics: `{output}/integration/prep/*.png`
  Non-integrated UMAPs: `{output}/integration/nonintegrated/*.png`
  Method UMAPs: `{output}/integration/{harmony,rpca,cca}/*.png`

Projection onto the atlas is **not** part of `hkoca pipeline`; run
`hkoca projection map` separately after integration or annotation.

### Annotation (Snapseed)

Uses the same Python stack as harmonize (`hkoca_harmonize`), with
[hassansaei/snapseed](https://github.com/hassansaei/snapseed) installed from GitHub
(not PyPI).

```bash
conda activate hkoca_harmonize
pip install -e .

# Editable copy of the packaged marker hierarchy
hkoca annotation markers export ./snapseed_markers.yaml

# QC-filtered or post-integration h5ad; default Leiden resolutions 0.4, 0.6, 1.0
hkoca annotation run \
  --input sample_filtered.h5ad \
  --output-dir results/annotation \
  --markers ./snapseed_markers.yaml
```

Outputs under `--output-dir`: `annotated_obj/*_annotated.h5ad`, `clustered/*_clustered.h5ad`,
optional `figures/`. UMAP PNGs are written by default (flat under `figures/`, no
per-sample subfolder) for Leiden and Level 1-3 per resolution; pass `--no-save-plots`
to skip. Each annotated object stores labels for all three resolutions
(`leiden_res_0.4` / `Level_*_res0.4`, etc.).

### Integration (Seurat prep)

Uses the dedicated R conda env `hkoca_integration` (Seurat, Harmony, clustree,
glmGamPoi, etc.). Set `R_MAX_VSIZE=150Gb` during prep for large objects.

```bash
conda env create -f conda/environment_integration.yaml
conda activate hkoca_integration
pip install -e .

# Optional GitHub-only R packages (also installed in the Docker image):
# Rscript -e "remotes::install_github('mojaveazure/seurat-disk', upgrade='never')"
# Rscript -e "remotes::install_github('samuel-marsh/scCustomize', upgrade='never')"

hkoca integration prep \
  --input-rds results/qc_filter/qc_filtered_rds/sample_harmonized_singlets_filtered.rds \
  --annotated-h5ad results/annotation/annotated_obj/sample_annotated.h5ad \
  --output-dir results/integration
```

Stage 1 (`prep`) normalizes RNA, runs SCTransform (glmGamPoi, regress
`percent.mito`), PCA elbow, clustree, and silhouette-based resolution selection.
Cell-type colors are shared with annotation via `hkoca/config/celltype_colors.yaml`
(same hex codes in Python UMAPs and R ggplot/dittoSeq plots).
Outputs: `prep/sct_prepared.rds`, `prep/elbow_plot.png`, `prep/clustree_sct.png`,
`prep/silhouette_sct.png`, `nonintegrated/*.png`, `tables/silhouette_scores.csv`,
`tables/resolution_summary.csv`. Silhouette selection ignores resolutions below
`silhouette_min_resolution` (default 0.4) to avoid overly coarse clusterings.
Per-method figures live under `{output}/harmony/`, `{output}/rpca/`, and
`{output}/cca/`; integrated RDS files are in `{output}/objects/`. Logs are
written directly under the integration output root (`integration_prep.log`,
`integration_methods.log`).

### Projection (atlas mapping)

Uses the dedicated GPU conda env `hkoca_projection` (PyTorch + scArches/scPoli).
This is **not** scanpy ingest: query cells are mapped by scPoli surgery onto the
HKOCA reference model, then labeled by prototype matching to
`Level_3_Integrated` (rolled up to Level 2 / Level 1).

Provide:

- `--query` — integration prep `sct_prepared.rds` (RNA counts are exported; SCT
  residuals are not used) or a query `.h5ad` with raw UMIs in `layers['counts']`
  or `.X`, gene symbols in `var_names`, and `obs['sample_id']`
- `--atlas` — HKOCA atlas `.h5ad` (`Level_*_Integrated`, `X_scpoli` /
  `X_umap_scpoli`)
- `--model-dir` — scPoli weights (`model_params.pt`, `attr.pkl`, `var_names.csv`)

```bash
conda env create -f conda/environment_projection.yaml
conda activate hkoca_projection
pip install -e .

hkoca projection map \
  --query results/integration/prep/sct_prepared.rds \
  --atlas reference/Master_Atlas_scPoli_Integrated_Reannotated_fullgenes.h5ad \
  --model-dir reference/scPoli_Reference_Model \
  --output-dir results/projection
```

Run surgery on a GPU node. RDS conversion uses `Rscript` from
`hkoca_integration` (or `HKOCA_RSCRIPT`). Resume skips an existing projected
h5ad; pass `--force` to re-run. Optional `--joint-umap` adds an exploratory
joint latent UMAP (atlas subsample + query); default figures place the query on
the **fixed atlas UMAP** by kNN in scPoli latent space.

`Level_3_uncert` is prototype uncertainty scaled per query batch (0 = confident,
1 = uncertain; relative, not an absolute atlas-wide score).

Outputs:

- `projected_obj/<query>_projected.h5ad` — query with `X_scpoli`,
  `Level_3_pred` / `Level_2_pred` / `Level_1_pred`, `Level_3_uncert`
- `tables/query_predictions_*.tsv`, `tables/projection_summary.json`
- `figures/query_on_atlas_umap_*.png` — query on atlas UMAP (labels, similarity,
  uncertainty)
- `figures/composition_*.png`, optional `figures/confusion_*.png`
- `models/scpoli_query_surgery/` — fine-tuned query embeddings
- `logs/projection.log`

Packaged defaults: `hkoca projection --print-config`.

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

## License

This project is released under the [MIT License](LICENSE).

When you use, copy, modify, or redistribute this toolkit (or parts of it), you must:

- Keep the copyright notice and MIT license text with the source or in documentation
- Credit the HKOCA / NephAtlas project and point to the source repository:
  https://github.com/hassansaei/HKOCA-tools
- If you distribute a modified version, clearly state that changes were made and retain the original copyright and license notice

Example attribution in a paper, report, or derivative tool:

> Analysis performed with HKOCA-tools (Human Kidney Organoid Cell Atlas / NephAtlas).
> Source: https://github.com/hassansaei/HKOCA-tools
