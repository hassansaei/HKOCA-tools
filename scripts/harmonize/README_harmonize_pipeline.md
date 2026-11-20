# data_loader_harmonized_21datasets

Config-driven loaders and gene harmonisation for the **21 atlas datasets** (scRNA-seq only).

This is the cleaned, scripted form of `data_loader_harmonized_21datasets.ipynb` (derived from the freeze notebook `data_loader_harmonized.ipynb`, without `d10_1016_j_stem_2021_12_010`).

## Files

| File | Role |
|---|---|
| `data_loader_harmonized_21datasets.py` | Entry point: helpers + per-study loaders + CLI |
| `data_loader_harmonized_21datasets.config` | Paths, reference genes, metadata defaults, study list, LMRH/clones |
| `data_loader_harmonized_21datasets.ipynb` | Notebook twin (same logic) |

## What it does

For each enabled study:

1. Loads raw data from `data_root/<study>/...`
2. Assigns shared `obs` metadata columns
3. Writes a raw `.h5ad` (when missing; overwrite is refused by default)
4. Harmonises genes to the Ensembl GRCh38.104 panel (protein-coding + lncRNA from the GTF)
5. Writes `<study>_harmonized.h5ad` under `harmonized_subdir`

**LMRH special case:** `LMRH_H5AD` is LMRH-only. G7 and TTC21B DMSO controls are loaded from their MTX folders and concatenated into the freeze file `dno_doi_kidney_organoid_LMRH_harmonized.h5ad` (the object QC used).

## Configuration

Edit `data_loader_harmonized_21datasets.config`:

### `[paths]`

- `data_root`: raw study folders
- `output_root`: parent of harmonised outputs
- `harmonized_subdir` / `h5ad_subdir`: subfolders under `output_root`
- `gtf_file`: Ensembl GTF used for the reference gene set

### `[reference]`

- `gene_biotypes`: GTF biotypes to keep (default `protein_coding,lncRNA`)

### `[metadata]`

- `standard_cols`: required `obs` columns
- `fill_value`: filler for missing columns (default `unknown`)
- `default_type` / `default_genome_build` / `default_disease`: defaults (scRNA-seq atlas)

### `[io]`

- `refuse_overwrite=true`: skip or refuse writing when a target `.h5ad` already exists

### `[studies]`

- `list`: the 21 atlas study keys
- `enabled=all` or a comma-separated subset

### `[lmrh]`, `[clone_g7]`, `[clone_ttc21b]`

- LMRH input subdir, protocols, and clone MTX locations / sample IDs

Study-specific sample maps (GSM lists, age maps, archive layouts) stay in the Python loaders; only shared / path-like settings live in the config.

## Usage

Use a compute node and the project `atlas` (or equivalent) env. Do not run the full atlas on the login node.

Keep caches under the labss tree (not `/home`):

```bash
export NUMBA_CACHE_DIR=/data-master/pure-workspace/labss/hmami/new_data/.cache/numba
export MPLCONFIGDIR=/data-master/pure-workspace/labss/hmami/new_data/.cache/matplotlib
export XDG_CACHE_HOME=/data-master/pure-workspace/labss/hmami/new_data/.cache
```

```bash
cd /data-master/pure-workspace/labss/hmami/new_data/HKOCA-tools/scripts/harmonize

# List configured studies
python data_loader_harmonized_21datasets.py --list-studies

# Run all enabled studies (from config)
python data_loader_harmonized_21datasets.py

# Custom config
python data_loader_harmonized_21datasets.py --config /path/to/my.config

# Subset
python data_loader_harmonized_21datasets.py --studies d10_1016_j_cell_2020_04_004,dno_doi_kidney_organoid_LMRH
```

## Outputs

Under `output_root/harmonized_subdir/` (default `Updated_data/updated_harmonized_files/`):

- `{study}_harmonized.h5ad` for each study
- For LMRH: `{study}_only_harmonized.h5ad` (LMRH alone) then `{study}_harmonized.h5ad` (freeze: LMRH + G7 + TTC21B)
- Clone controls also write `{clone_controle_*}_harmonized.h5ad` when `include_clone_controls=true`

## Relation to QC / annotation

1. Harmonised `.h5ad` → Seurat `.rds` (separate conversion step)
2. QC: `HKOCA-tools/scripts/qc/scripts/QC_scdbl_Combined.R`
3. Annotation: `snapseed_annotation_after_harmonization.ipynb` on QC-filtered h5ads

## Notes

- scRNA-seq only (`type=scRNAseq`); no snRNA-seq claims
- Default IO policy refuses overwriting existing `.h5ad` outputs
- Sibling generic pipeline: `harmonize_pipeline.py` + `harmonize.config` (CSV-driven). This 21-dataset script keeps the custom per-study loaders from the freeze notebook.
