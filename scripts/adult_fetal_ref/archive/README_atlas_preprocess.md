# Adult/Fetal Kidney Atlas Marker Extraction Pipeline

A configuration-driven Scanpy pipeline for deriving **canonical, cross-atlas
cell-type marker gene sets** from paired adult and fetal human kidney
single-cell reference atlases. The pipeline standardizes cell-type
nomenclature across developmental stages, performs differential expression
testing, applies a reproducible statistical cutoff to marker lists, and
generates a full deliverables bundle for downstream annotation, reference
mapping, or atlas-integration work.

`adult_fetal_atlas_preprocess.py` is entirely parameterized by
`adult_fetal_atlas_preprocess_config.yaml` : no thresholds, paths, or
mappings are hard-coded, making the pipeline reusable across atlas
versions or other tissues without touching the code.

---

## Table of Contents

1. [Reference Data](#1-reference-data)
2. [Pipeline Architecture](#2-pipeline-architecture)
3. [Key Features](#3-key-features)
4. [The Stringent Score & Elbow Cutoff](#4-the-stringent-score--elbow-cutoff)
5. [Installation](#5-installation)
6. [Usage](#6-usage)
7. [Configuration Reference](#7-configuration-reference)
8. [Output Structure](#8-output-structure)

---

## 1. Reference Data

This pipeline expects two **preprocessed** `.h5ad` files as input: one
adult human kidney reference and one fetal (developing) human kidney
reference, each with a categorical cell-type annotation column in `.obs`.

The cell-type vocabulary encoded in the config's `mappings` section
(`epithelial_cell_of_proximal_tubule`, `kidney_loop_of_Henle`,
`endothelial_cell_of_lymphatic_vessel`, `glomerular_capillary_endothelial_cell`,
etc.) is aligned with the Cell Ontology labels used by the following
CZ CELLxGENE Discover collections:

| Atlas | Source | CELLxGENE Collection |
|---|---|---|
| **Adult** | Xu *et al.* (2023) *Cell* : "Automatic cell-type harmonization and integration across Human Cell Atlas datasets" | https://cellxgene.cziscience.com/e/8dc29642-be06-46dc-a289-8075915745d4.cxg/ |
| **Fetal** | Hochane *et al.* (2019) *PLoS Biol* : "Single-cell transcriptomics reveals gene expression dynamics of human fetal kidney development" | https://cellxgene.cziscience.com/e/f95d8919-1f2a-405f-8776-bfecc0ab0f3f.cxg/ |

**To obtain the input files:**

1. Open the link above for each atlas.
2. Download the `.h5ad` export directly from the CELLxGENE Discover
   dataset page, or via the `cellxgene-census` / `cellxgene` Python API for
   programmatic access.
3. Place the downloaded files under your project's data directory and
   point `paths.adult_h5ad` / `paths.fetal_h5ad` in the config at their
   locations (relative to `paths.base_dir`).

> Always confirm that the cell-type annotation column and label spelling
> in the file you download match `columns.adult_celltype_col` /
> `columns.fetal_celltype_col` and the raw labels referenced in
> `mappings.adult_canonical_raw` / `mappings.fetal_canonical_raw` :
> CELLxGENE dataset versions are periodically revised and annotation
> granularity can change between releases.

---

## 2. Pipeline Architecture

The pipeline executes as a single linear workflow orchestrated by `main()`:

```
Config & I/O setup
        │
        ▼
Load adult.h5ad + fetal.h5ad
        │
        ▼
Standardize cell-type column → obs["cell_type"]
        │
        ▼
Checkpoint overview (cell/gene counts, type frequencies)
        │
        ▼
Endothelial subtype consolidation
        │
        ▼
Low-abundance cell-type QC filter
        │
        ▼
Wilcoxon rank-sum marker extraction (one-vs-rest, per atlas)
        │
        ├──▶ Optional elbow-threshold diagnostic (single test cell type per atlas)
        │
        ▼
Canonical cell-type label mapping (adult ∪ fetal)
        │
        ▼
Per-cell-type marker filtering → stringent score → elbow cutoff → TSV export
        │
        ▼
Pan-endothelial consensus marker merge
        │
        ▼
Deliverables generation:
   • label-mapping table
   • top-N marker tables (adult / fetal)
   • cell-type presence comparison table
   • adult ↔ fetal marker overlap table
   • pseudobulk correlation heatmap + QC checks
   • run-parameter record (YAML)
```

Every stage logs structured progress and QC information to stdout, so the
pipeline is safe to run unattended (e.g. in a batch job or workflow
manager) with full traceability from the log alone.

---

## 3. Key Features

- **Fully config-driven** : every path, threshold, and label mapping lives
  in YAML; the script itself never needs to be edited between runs or
  atlas versions.
- **Cross-atlas label harmonization** : raw, atlas-specific cell-type
  labels (which differ in naming convention between the adult and fetal
  references) are mapped onto a single shared canonical vocabulary,
  enabling direct adult-vs-fetal comparison.
- **Endothelial subtype consolidation** : fine-grained vascular subtypes
  (lymphatic, glomerular capillary, peritubular capillary, arterial) are
  optionally collapsed into a unified `endothelial_cell` group prior to
  differential expression testing, and a separate **pan-endothelial
  consensus marker list** is also generated by combining evidence across
  the individual subtypes.
- **Statistically principled marker selection** : markers are not simply
  top-*N* by fold change; each candidate list is filtered by
  significance/effect-size thresholds, ranked by a composite score, and
  truncated at a **data-driven elbow point** rather than an arbitrary
  fixed count.
- **Built-in QC diagnostics** : automatic checks flag missing expected
  cell lineages and low adult/fetal pseudobulk correlation for supposedly
  matching cell types, surfacing likely annotation or mapping problems
  before they propagate downstream.
- **Reproducibility record** : every run writes out the exact filtering
  thresholds, tolerances, and scoring formula used, alongside the marker
  tables themselves.
- **Cross-atlas visual summary** : a pseudobulk transcriptomic correlation
  heatmap, annotated with shared marker-gene counts, gives an at-a-glance
  view of adult/fetal cell-type correspondence.

---

## 4. The Stringent Score & Elbow Cutoff

Marker gene selection is the statistical core of this pipeline. Rather
than exporting a fixed top-*N* gene list per cell type, every candidate
gene list goes through a four-stage funnel:

### Step 1 : Baseline significance/effect-size filter
Genes must pass, simultaneously:

| Criterion | Config key | Purpose |
|---|---|---|
| `pct_in ≥ min_pct_in` | `deliverables.min_pct_in` | The gene must be expressed in a substantial fraction of cells *within* the target cell type. |
| `logFC ≥ min_logfc` | `deliverables.min_logfc` | The gene must show a minimum magnitude of up-regulation. |
| `pval_adj < max_pval_adj` | `deliverables.max_pval_adj` | The result must be statistically significant after Benjamini–Hochberg correction. |
| gene does not start with `exclude_gene_prefix` | `deliverables.exclude_gene_prefix` | Excludes uninformative or technically confounded gene families (e.g. mitochondrial genes) from marker consideration. |

### Step 2 : Stringent score

Surviving genes are ranked by a composite **stringent score**:

```
stringent_score = logFC × (pct_in ^ 2) / (pct_out + tolerance)
```

| Term | Meaning |
|---|---|
| `logFC` | Log fold-change of expression in the target group vs. all other groups : captures effect *magnitude*. |
| `pct_in` | Fraction of cells within the target cell type expressing the gene : captures within-group *penetrance*. Squaring this term disproportionately rewards genes that are broadly expressed across the target population, rather than only in a subset of cells. |
| `pct_out` | Fraction of cells in *all other* cell types expressing the gene : captures background/off-target expression (specificity). |
| `tolerance` | A stabilizing constant added to the denominator (see below). |

The result is a single score that jointly rewards genes that are strongly
up-regulated, expressed consistently throughout the target population, and
rarely expressed elsewhere : the combination of properties that makes a
gene a good discriminative marker, rather than optimizing for any one
property in isolation.

**Why `tolerance` is needed:** The **`tolerance`** parameter serves a critical dual purpose in the stringent marker scoring formula: `logFC * (pct_in²) / (pct_out + tolerance)`. 

Mathematically, it prevents division-by-zero errors in cases where a gene has absolutely zero background expression (`pct_out = 0`). Biologically, it acts as a stabilizing buffer against ambient RNA "soup":a ubiquitous technical artifact in single-cell sequencing where free-floating transcripts cause negligible, false-positive reads across unintended clusters. 

Without this buffer, the algorithm becomes ruthlessly hypersensitive to background noise, over-penalizing highly specific markers for microscopic, biologically irrelevant differences (exp, 1% vs. 2% background expression). 

Crucially, this parameter acts as a "biological forgiveness" dial and should be tuned depending on the reference atlas you are mapping against:

* **For Adult Atlases (Low Tolerance: exp,  `0.01` – `0.05`)**
    Adult tissues are fully mature and terminally differentiated, meaning their cells exist in highly specialized, discrete states. Setting a low tolerance acts as a strict, ruthless filter. It aggressively penalizes any gene that "bleeds" into other clusters, ensuring that only the most exclusive, canonical mature markers survive the elbow cutoff.

* **For Fetal & Organoid Atlases (High Tolerance: exp,  `0.15` – `0.25`)**
    Organoids model active fetal development, which is characterized by continuous differentiation trajectories rather than rigid, discrete bins. A transitioning cell (exp,  a nephron progenitor committing to a podocyte lineage) will naturally share a gradient of overlapping transcripts with its neighbors. A higher tolerance provides the necessary mathematical forgiveness, ensuring that vital developmental markers are not wrongly discarded just because they are naturally co-expressed across neighboring lineages on the UMAP.
### Step 3 : Elbow-point detection

Genes are sorted by `stringent_score` in descending order and capped to
the top `deliverables.max_n` candidates. Rank and score are each
normalized to `[0, 1]`, and the gene at maximum distance from the
diagonal line (`x_norm + y_norm = 1`) of the rank-vs-score curve is taken
as the natural inflection point ("elbow") : the point beyond which
additional genes add diminishing discriminative value.

### Step 4 : Final cutoff

The exported marker list keeps:

```
n_kept = max(elbow_rank, deliverables.min_kept)
```

i.e. the pipeline never drops below a configured floor of markers per cell
type (`min_kept`), even for cell types whose elbow point falls
unexpectedly early, but otherwise defers entirely to the data-driven
elbow detection.

This same four-step procedure is applied consistently in three places:
the optional diagnostic tool (`check_threshold_elbow`), the per-cell-type
TSV export, and the deliverables' top-N marker tables : so diagnostics
run before a full export are guaranteed to reflect the exact cutoff logic
used in the final output.

---

## 5. Installation

```bash
pip install scanpy pandas numpy matplotlib seaborn pyyaml
```

Python ≥3.9 and Scanpy ≥1.9 are recommended. No GPU is required.

---

## 6. Usage

```bash
python adult_fetal_atlas_preprocess.py --config adult_fetal_atlas_preprocess_config.yaml
```

### Command-line options

| Flag | Required | Description |
|---|---|---|
| `--config` | Yes | Path to the pipeline's YAML configuration file. This is the pipeline's only entry point for parameters : there are no other CLI flags. |

All run behavior : input/output locations, QC thresholds, differential
expression settings, scoring parameters, and label mappings : is
controlled exclusively through the YAML file passed to `--config`,
allowing multiple pipeline configurations (e.g. different tissues, atlas
versions, or threshold sensitivities) to be maintained side-by-side
without modifying the script.

Logging is written to stdout in real time, including per-atlas dataset
overviews, QC filtering decisions, diagnostic elbow reports (if enabled),
and a final completion message : suitable for redirecting to a log file
in automated environments (`... 2>&1 | tee run.log`).

---

## 7. Configuration Reference

### `paths`
| Key | Description |
|---|---|
| `base_dir` | Root directory that all other paths are resolved against. |
| `adult_h5ad` | Path (relative to `base_dir`) to the preprocessed adult `.h5ad` reference. |
| `fetal_h5ad` | Path (relative to `base_dir`) to the preprocessed fetal `.h5ad` reference. |
| `out_ref` | Output directory reserved for reference-level artifacts. |
| `out_adult` | Output directory for adult per-cell-type marker TSVs. |
| `out_fetal` | Output directory for fetal per-cell-type marker TSVs. |
| `out_req` | Output directory for the final deliverables bundle. |

### `columns`
| Key | Description |
|---|---|
| `adult_celltype_col` | `obs` column in the adult `.h5ad` holding cell-type labels. |
| `fetal_celltype_col` | `obs` column in the fetal `.h5ad` holding cell-type labels. |

### `scanpy_settings`
| Key | Description |
|---|---|
| `verbosity` | Scanpy logging verbosity (0 = errors only, 3 = debug). |
| `dpi` | Figure resolution for Scanpy and Matplotlib output. |
| `facecolor` | Background color for generated figures. |
| `frameon` | Whether plots draw an axes frame/border. |

### `qc`
| Key | Description |
|---|---|
| `min_cells` | Minimum number of cells a cell type must have to be retained prior to marker calling. |
| `min_corr` | Minimum acceptable pseudobulk Pearson correlation between adult and fetal profiles of the same canonical cell type; flags likely mapping issues below this threshold. |
| `expected_lineages` | Canonical lineage names expected to be present in both atlases; missing entries are flagged as QC warnings. |

### `wilcoxon`
| Key | Description |
|---|---|
| `n_genes` | Number of genes per group requested from the Wilcoxon rank-sum test. |

### `run_diagnostics` / `diagnostics`
| Key | Description |
|---|---|
| `run_diagnostics` | Boolean toggle for running the elbow-threshold diagnostic report. |
| `diagnostics.adult_test_cell_type` | Cell type used for the adult diagnostic pass. |
| `diagnostics.fetal_test_cell_type` | Cell type used for the fetal diagnostic pass. |

### `elbow_params`
| Key | Description |
|---|---|
| `adult_tolerance` | Stringent-score denominator tolerance for the adult atlas. |
| `fetal_tolerance` | Stringent-score denominator tolerance for the fetal atlas. |

### `pan_endothelial`
| Key | Description |
|---|---|
| `files_to_merge` | Adult per-cell-type marker TSV filenames to combine into a consensus endothelial marker list. |
| `top_n` | Number of top-ranked consensus genes to retain in the merged output. |
| `merged_group_name` | Label assigned to the merged marker group. |
| `merged_file_name` | Output filename for the merged pan-endothelial marker TSV. |

### `deliverables`
| Key | Description |
|---|---|
| `min_pct_in` | Minimum within-group expression fraction for marker inclusion. |
| `min_logfc` | Minimum log fold-change for marker inclusion. |
| `max_pval_adj` | Maximum adjusted p-value for marker inclusion. |
| `exclude_gene_prefix` | Gene-name prefix excluded from marker consideration (e.g. mitochondrial genes). |
| `export_n` | Number of top-ranked genes per cell type in the deliverable marker tables. |
| `overlap_n` | Number of top-ranked genes per cell type used for adult/fetal overlap and heatmap computation. |
| `max_n` | Maximum candidate gene pool size per cell type before elbow detection. |
| `min_kept` | Minimum number of genes always retained per cell type. |

### `plotting`
| Key | Description |
|---|---|
| `save_figs` | Whether to generate and save the pseudobulk correlation heatmap and its associated QC checks. |

### `mappings`
| Key | Description |
|---|---|
| `endothelial_consolidation` | Raw-label → consolidated-label map applied to both atlases before marker extraction, unifying endothelial subtypes. |
| `adult_canonical_raw` | Override map from normalized adult raw labels to final canonical cell-type names. |
| `fetal_canonical_raw` | Override map from normalized fetal raw labels to final canonical cell-type names, including renaming of fetal-specific ontology terms to align with adult nomenclature. |

---

## 8. Output Structure

```
out_adult/
├── adult_<canonical_cell_type>.tsv        # one file per canonical adult cell type
└── <merged_file_name>                     # pan-endothelial consensus markers

out_fetal/
└── fetal_<canonical_cell_type>.tsv        # one file per canonical fetal cell type

out_req/
├── cell_type_label_mapping.tsv            # raw → canonical label map, both atlases
├── adult_kidney_markers_top20.tsv         # top-N ranked markers per adult cell type
├── fetal_kidney_markers_top20.tsv         # top-N ranked markers per fetal cell type
├── cell_type_comparison.tsv               # canonical cell-type presence, adult vs. fetal
├── marker_overlap_adult_vs_fetal.tsv      # shared marker genes for matching cell types
├── marker_overlap_figure.png              # pseudobulk correlation heatmap + shared-gene counts
└── marker_method_params.yaml              # exact thresholds/formula used for this run

out_ref/                                   # reserved for reference-level artifacts
```

Each marker TSV includes, at minimum: `group`, `gene`, `stringent_score`,
`pct_in`, `pct_out`, `logFC`, and `pval_adj`, giving downstream users
everything needed to independently re-rank, re-filter, or audit the
marker calls.
