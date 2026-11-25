# Upstream — building the adult and fetal references

These scripts build (or audit) the primary atlases used in section 00. Outputs stay under `organoid_vs_primary_atlas/00_reference_umaps/`.

## Adult (diseased tissue removed)

Source: HCKA `Single_cell_human_kidney_cell_atlas.h5ad`.

We retain cells with normal disease status and healthy or adjacent sampling sites, and exclude diseased samples. `process_adult.py` audits that filter, writes cell-count tables, and a summary figure of kept vs excluded cells.

```bash
python scripts/00_reference_umaps/upstream/process_adult.py
```

## Fetal (local inputs → weeks 9–17 Harmony)

All fetal work is local to this project:

| Path | Role |
|---|---|
| `00_reference_umaps/inputs/fetal/{stewart,hochane,tran}.h5ad` | source AnnData (place here) |
| `00_reference_umaps/objects/fetal_weeks_9_17_harmony_integrated.h5ad` | integrated reference |
| `00_reference_umaps/figures/fetal_integration/` | UMAPs |
| `00_reference_umaps/tables/fetal_integration/` | counts + summary |

Author annotations are remapped to `harmonized_celltype`. Only cells with a **numeric** post-fertilization week in **9–17** are kept; embryonic stage, Carnegie stage, and other non-week labels are dropped. Then Harmony (batch = dataset) → UMAP + Leiden.

```bash
# 1) put stewart.h5ad, hochane.h5ad, tran.h5ad in 00_reference_umaps/inputs/fetal/
# 2) run
python scripts/00_reference_umaps/upstream/process_fetal.py
```

`process_fetal.py` calls `integrate_fetal_harmony.py`.


```bash
Rscript scripts/00_reference_umaps/upstream/run_devkidcc_tran.R
```

(Update that R script’s `ROOT`/`RAW`/`OUT` to the local `inputs/fetal/` paths if rebuilding Tran from scratch.)
