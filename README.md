# HKOCA — Analysis Archive

This branch stores the **analysis scripts and notebooks** used to build the Human Kidney Organoid Cell Atlas. It is the working archive for manuscript Materials and Methods: each stage has a dedicated folder with code, config, and a local README.

The reusable toolkit lives on `main`. Use this branch to recover how the atlas was generated and to lock parameters for the paper.

## Analysis Pipeline (Order Used for the Atlas)

```
1. harmonize          Gene-space harmonization to GRCh38.104 (~24,100 genes)
2. QC                 Doublet detection (scDblFinder) + QC filters + h5ad export
3. Annotation         Hierarchical Snapseed cell-type labels (Level 1–3)
4. Integration        Method benchmark (scIB) → full-atlas scPoli integration
5. Reannotation       High-resolution Leiden + majority-vote Level_*_Integrated
6. References         Adult/fetal marker extraction and organoid↔reference mapping
```

## Repository Layout

| Path | Role |
|------|------|
| [`scripts/`](scripts/) | Primary analysis code for the atlas (canonical) |
| [`notebooks/`](notebooks/) | Exploratory / legacy notebooks (see note below) |
| [`data/meta/`](data/meta/) | Small metadata assets (gene table, older marker YAML) |
| [`plots/`](plots/) | Example QC and integration comparison figures |
| [`docker/`](docker/) | Minimal Docker stub from early setup |
| [`src/Snapseed/`](src/Snapseed/) | Older Snapseed marker YAML copy |
| `containers/`, `test/` | Placeholders |

## Scripts

### `scripts/harmonize/`

Loads study samples from a metadata CSV, attaches standardized metadata, concatenates by study, and harmonizes each object to the Ensembl GRCh38.104 protein-coding + lncRNA gene set (optional Seurat RDS export).

| File | Purpose |
|------|---------|
| `harmonize_pipeline.py` | Main pipeline |
| `harmonize.config` | Paths and summary plot settings |
| `meta/datasets_metadata.csv` | Sample / study manifest |
| `README_harmonize_pipeline.md` | Full documentation |

### `scripts/QC/`

R pipeline for ghost-cell prefilter, scDblFinder doublet calling, per-sample QC thresholds, audit plots, and native RDS→h5ad conversion via `anndata`/`reticulate`.

| File | Purpose |
|------|---------|
| `QC_scdbl_Combined.R` | Doublet + QC + optional h5ad export |
| `qc_config.dcf` | Stage paths and doublet/QC settings |
| `qc_filters_template.csv` | Per-dataset threshold method template |
| `environment_qc.yaml` | Conda environment |
| `QC_README.md` | Full documentation |

### `scripts/Annotation/`

Hierarchical Snapseed annotation of QC-filtered, harmonized `.h5ad` objects (per-dataset Leiden, marker hierarchy, optional manual overrides and cross-dataset summaries).

| File | Purpose |
|------|---------|
| `snapseed_annotation.ipynb` | End-to-end annotation notebook (**canonical**) |
| `snapseed_annotation_config.yaml` | Paths and Scanpy parameters |
| `snapseed_markers_v4.yaml` | Level 1–3 marker hierarchy |
| `README.md` | Full documentation |

### `scripts/Integration/`

| Subfolder | Purpose |
|-----------|---------|
| [`Benchmark/`](scripts/Integration/Benchmark/) | Twelve-method integration benchmark with scIB metrics (`workflow_scheduler`, Python/R methods, visualisation) |
| [`Integration/`](scripts/Integration/Integration/) | Production full-atlas scPoli training (`run_integration_scpoli_parameters.py`) |
| [`Reannotation/`](scripts/Integration/Reannotation/) | Leiden resolution 10 majority-vote Level_*_Integrated labels |
| [`conda/`](scripts/Integration/conda/) | Shared `environment.yml` for integration stages |

See `scripts/Integration/Integration/integration_README.md`, `Benchmark/README.md`, and `Reannotation/reannotation_README.md`.

### `scripts/adult_fetal_ref/`

Derives canonical adult and fetal kidney marker sets from CELLxGENE reference atlases (Wilcoxon + stringent score + elbow cutoff).

| File | Purpose |
|------|---------|
| `adult_fetal_atlas_preprocess.py` | Marker extraction pipeline |
| `adult_fetal_atlas_preprocess_config.yaml` | Paths, mappings, thresholds |
| `README_atlas_preprocess.md` | Full documentation |

### `scripts/organoid_mappping/`

Cluster-level marker extraction for organoid datasets and Pearson correlation mapping to adult/fetal reference marker tables. Folder name retains the historical spelling `organoid_mappping`.

| File | Purpose |
|------|---------|
| `organoid_mapping.py` | Phase 1 markers + Phase 2 reference mapping |
| `organoid_mapping_config.yaml` | Resolutions and scoring parameters |
| `environment_organoid_mapping.yml` | Conda environment |
| `organoid_mapping_README.md` | Full documentation |

### `scripts/old/`

Archived earlier integration, QC, and annotation prototypes. **Not** the production path for the manuscript; kept for provenance only.

## Notebooks

| Notebook | Notes |
|----------|-------|
| [`notebooks/snapseed_annotation.ipynb`](notebooks/snapseed_annotation.ipynb) | Duplicate of the Annotation notebook; prefer `scripts/Annotation/` |
| [`notebooks/data_loaders.ipynb`](notebooks/data_loaders.ipynb) | Early sfaira / GTF harmonization experiments; contains machine-local paths |

## Materials and Methods Checklist

Open items for manuscript Methods (each is a GitHub issue; collaborator @hythem03 is tagged):

| Issue | Topic |
|------:|-------|
| [#18](https://github.com/hassansaei/HKOCA-tools/issues/18) | Atlas cohort, inclusion criteria, sample inventory |
| [#19](https://github.com/hassansaei/HKOCA-tools/issues/19) | Gene harmonization (GRCh38.104) parameters |
| [#20](https://github.com/hassansaei/HKOCA-tools/issues/20) | QC and doublet settings (incl. CellBender provenance) |
| [#21](https://github.com/hassansaei/HKOCA-tools/issues/21) | Snapseed markers, resolutions, manual overrides |
| [#22](https://github.com/hassansaei/HKOCA-tools/issues/22) | Integration benchmark design and method choice |
| [#23](https://github.com/hassansaei/HKOCA-tools/issues/23) | Production scPoli hyperparameters and Level-3 reannotation |
| [#24](https://github.com/hassansaei/HKOCA-tools/issues/24) | Adult/fetal references and organoid↔reference mapping |
| [#25](https://github.com/hassansaei/HKOCA-tools/issues/25) | Environments, seeds, compute, data deposition |

## How to Use This Archive

1. Start from the stage README in each `scripts/` subfolder.
2. Point configs at your local data roots (committed configs may still contain cluster-absolute paths).
3. Prefer stage conda YAMLs over ad hoc notebook installs when reproducing results.
4. Treat `scripts/old/` and `notebooks/data_loaders.ipynb` as historical only.

## License

See the repository `LICENSE` on `main` (BSD 3-Clause) unless this branch carries a local copy.
