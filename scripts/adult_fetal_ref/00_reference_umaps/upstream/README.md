# Upstream — building the adult and fetal references

These scripts document (and optionally rebuild) the primary atlases used in section 00. They sit upstream of the UMAP figure scripts.

## Adult (diseased tissue removed)

Source: HCKA `Single_cell_human_kidney_cell_atlas.h5ad`.

We retain cells with normal disease status and healthy or adjacent sampling sites, and exclude diseased samples. `process_adult.py` audits that filter, writes cell-count tables, and a summary figure of kept vs excluded cells. The full healthy+adjacent scANVI object and organoid projection were built earlier in `human_kidney_reference_projection/`; pass `--rebuild-projection` only if you need to regenerate that heavy run.

## Fetal (collection → annotation → Harmony)

Three inputs:

- **Stewart** and **Hochane** — prebuilt `.h5ad` files already in the fetal project folder  
- **Tran** (GSE124472) — raw 10x MTX → QC (`tran_qc`) → DevKidCC in R (`run_devkidcc_tran.R`) → assemble `tran.h5ad`

Author annotations are remapped to `harmonized_celltype`, weeks outside 9–17 are dropped, and the three datasets are Harmony-integrated. `process_fetal.py` is the single Python driver (`--stage sync|tran_qc|assemble|integrate|all`). By default it only syncs the already-computed integration UMAPs (dataset, week, cell type, Leiden) into `00_reference_umaps/figures/fetal_integration/`.

```bash
python scripts/00_reference_umaps/upstream/process_adult.py
python scripts/00_reference_umaps/upstream/process_fetal.py
# Rscript scripts/00_reference_umaps/upstream/run_devkidcc_tran.R   # only when rebuilding Tran
```
