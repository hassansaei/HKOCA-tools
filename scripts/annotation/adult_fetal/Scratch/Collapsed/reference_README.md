# Unified Adult-Fetal Kidney Reference Pipeline

This folder contains a production pipeline script designed with explicit, testable, and reproducible logic.

- Script: `adult_fetal_ref_overlap.py`
- Config: `config.yaml`

This README is intentionally explanatory. It focuses on:

1. What the pipeline actually does in each phase.
2. Why specific logic was chosen over simpler alternatives.
3. What is robust/sophisticated about the implementation.
4. How to run it safely and interpret outputs.

## 1) What the script introduces


- Explicit phase boundaries and contracts.
- Fail-fast checks for missing/invalid inputs.
- Resume logic with deterministic fingerprints.
- Structured diagnostics and reproducibility artifacts.
- Config-driven behavior with cross-section validation.


## 2) What You Get (Practically)


- Repeatable reruns with automatic invalidation when relevant settings/data/code change.
- Traceability (manifest + logs + parameter snapshots).
- Safer downstream-only runs (upstream compatibility guards).
- More robust label harmonization and overlap interpretation.
- Cleaner output separation between deliverables and diagnostics.


## 3) Phase 1 Deep Dive: Preprocess, Marker Extraction, Fetal Rescue

### 3.1 What it does

Phase 1:

- Loads adult and fetal references (`.h5ad`).
- Resolves gene symbol source robustly.
- Filters cells/genes, normalizes/log-transforms, computes HVGs/PCA/neighbors/UMAP.
- Extracts markers using Wilcoxon + BH correction.
- Applies optional stress mitigation logic (adult).
- Applies fetal rescue logic for generic fetal labels.
- Writes:
  - `adult_kidney_markers_all.tsv` (aux)
  - `fetal_kidney_markers_all.tsv` (aux)
  - `fetal_kidney_markers_top{K}.tsv` (main)
  - preprocessed h5ad checkpoints
  - 2-panel UMAP figure

Important: the fetal top table uses `hpa_validation.top_n_validated` as `K`, not `marker_extraction.top_n_markers`.

Why: it aligns top-N size across adult validated and fetal marker lists for downstream comparison consistency.


### 3.2 Sophisticated logic in Phase 1

#### Gene symbol source decision (`var_names` vs `feature_name`)

The script does not blindly overwrite with `feature_name`.
It compares symbol quality heuristically:

- Ensembl-like ratio in each source.
- Numeric-like ratio in each source.
- Invalid/blank rate in candidate source.

Switch only when there is clear quality gain; otherwise keep current names.
Then it logs and resolves duplicates using `var_names_make_unique()`.

Why this logic: naive overwrite can corrupt marker identity and create silent gene name collisions.


#### Adult stress handling in two layers

1. Pre-DE mitigation (`adult_stress_mitigation_mode`):
   - `none`
   - `exclude_high_stress` (quantile-based)
   - `regress_out`
2. Post-DE stress-gene exclusion from adult marker tables.

Why: stress confounding can enter at both cell-level and marker-level; dual control is more robust than either alone.


#### Marker quality gating beyond p-value/logFC

Additional quality gates:

- `min_pct_cells`
- `max_pct_rest`

Why: avoids markers that are statistically significant but biologically weak or non-specific.


#### Fetal rescue with confidence model

Generic fetal labels are scored using configured distal and LoH marker sets.
Rescue uses:

- Score margin threshold (`fetal_rescue_min_margin`)
- Absolute confidence threshold (`fetal_rescue_min_score`)
- Ambiguous strategy (`keep_original` or `assign_label`)

Then residual generic labels are dropped before final fetal marker extraction.

Why: argmax-only assignment overstates certainty in close-score or low-signal cells.


## 4) Phase 2 Deep Dive: Adult HPA Validation And Joint Ranking

### 4.1 What it does

Phase 2:

- Loads adult all-marker table from Phase 1.
- Loads HPA TSV into a gene-indexed lookup.
- Validates each marker against mapped HPA cell type expression (`nCPM`).
- Applies threshold + manual override logic.
- Builds a weighted joint score and ranks within each adult cell type.
- Writes:
  - validation audit TSV
  - ranked candidate diagnostics TSV
  - top validated adult markers in both aux and output root


### 4.2 Core validation logic

For each marker row:

1. Resolve adult cell type mapping to HPA label.
2. Find gene in HPA.
3. Find mapped HPA cell-type entries.
4. Compute:
   - `hpa_nCPM` (max nCPM for mapped type)
   - `hpa_specificity` (mapped max / sum of per-type maxima)
5. Mark as validated if `hpa_nCPM >= validation_ncpm_threshold`.

Manual override behavior:

- Cell type must still have explicit HPA mapping.
- Listed manual markers can be forced validated.

Why this logic: keeps manual curation explicit but prevents hidden fallback behavior.


### 4.3 Joint score rationale

Final ranking combines:

- extraction strength (from DE rank)
- HPA expression strength
- HPA specificity

via weighted sum of scaled components:

`weight_extraction_strength`, `weight_hpa_strength`, `weight_hpa_specificity`.

Why: single-metric ranking (only DE or only HPA) often over-favors one axis and reduces practical interpretability.


## 5) Phase 3 Deep Dive: Mapping-aware Overlap And YAML Export

### 5.1 What it does

Phase 3:

- Loads adult validated markers (Phase 2) and fetal all markers (Phase 1).
- Harmonizes labels via `overlap.harmonization_map`.
- Builds mapping-aware cell type comparison tables.
- Computes exhaustive adult x fetal overlap (`n_shared`, Jaccard, shared genes).
- Applies optional acceptance criteria and fails run if criteria are violated.
- Generates overlap figure and two YAML exports.


### 5.2 Why mapping-aware overlap (not exact-label-only)

Human label granularity differs between references. Exact-label intersection alone is too strict and often biologically misleading.

The script combines:

- Exact label matches.
- Manual mapping seed (config or sensible default).

Then computes shared/adult-only/fetal-only summaries under this mapping-aware perspective.


### 5.3 Acceptance criteria as quality gates

Optional `overlap.acceptance_criteria` can enforce expected overlap behavior:

- minimum best mapped Jaccard
- minimum best mapped shared genes
- expected ranges for shared/adult-only/fetal-only counts

If configured criteria fail, Phase 3 raises an error.

Why: it turns overlap expectations into explicit, machine-checkable QA.


### 5.4 YAML exports and ordering policy

The script writes:

- Detailed YAML: includes per-cell-type adult/fetal/shared/adult-only/fetal-only marker lists.
- Official YAML: includes only non-empty overlap genes per cell type.

Marker ordering is rank-preserving (not alphabetical).
Adult markers are capped at 25 per harmonized cell type in YAML; fetal markers are uncapped.

Why: rank-preserving export keeps biological signal priority visible to downstream tools/users.


## 6) Resume, Rerun, And Upstream-Compatibility Logic

The resume mechanism uses `<aux_dir>/.meta.json`.

Additional guard:

- If running Phase 2/3 without rerunning upstream phases, fingerprints are checked for upstream compatibility.
- If incompatible, execution fails with an explicit rerun message.

Why: prevents accidentally combining stale upstream artifacts with new downstream config/code.


## 7) Placeholders And Directory Semantics

Use these placeholders in examples:

- `<pipeline_dir>`: directory containing script and config
- `<config_path>`: path to config file
- `<output_dir>`: value of `paths.output_dir`
- `<aux_dir>`: `<output_dir>/<aux_subdir>`
- `<figures_dir>`: `<aux_dir>/<figures_subdir>`
- `<diagnostics_dir>`: `<aux_dir>/<diagnostics_subdir>`


## 8) Output Model (Main vs Aux)

### Main deliverables in `<output_dir>`

- `fetal_kidney_markers_top{K}.tsv`
- `adult_kidney_markers_top{K}_validated.tsv`
- `cell_type_comparison.tsv`
- `marker_overlap_adult_vs_fetal.tsv`
- `marker_overlap_figure.png`

`K = hpa_validation.top_n_validated`


### Provenance and diagnostics in `<aux_dir>`

- `dataset_manifest.tsv`
- `marker_method_params.yaml`
- `adult_kidney_markers_all.tsv`
- `fetal_kidney_markers_all.tsv`
- `fetal_generic_kidney_rescue_assignments.tsv`
- `adult_preprocessed.h5ad`
- `fetal_preprocessed.h5ad`
- `cell_type_label_mapping.tsv`
- `snapseed_markers_adult_fetal_overlap_detailed.yaml`
- `snapseed_markers_adult_fetal_overlap_official.yaml`
- `logs/pipeline_<RUN_ID>.log`
- `run_manifests/run_manifest_<RUN_ID>.json`
- diagnostics subfolders for HPA and overlap tables


## 9) How To Run (Concise)

From `<pipeline_dir>`:

```bash
python adult_fetal_ref_overlap.py --config <config_path> --dry-run
python adult_fetal_ref_overlap.py --config <config_path>
python adult_fetal_ref_overlap.py --config <config_path> --verify
```

Selected phases:

```bash
python adult_fetal_ref_overlap.py --config <config_path> --only-phases 1
python adult_fetal_ref_overlap.py --config <config_path> --only-phases 2 3
python adult_fetal_ref_overlap.py --config <config_path> --skip-phases 3
```

Force rerun (ignore resume state):

```bash
python adult_fetal_ref_overlap.py --config <config_path> --no-resume
```


## 10) Requirements

Python:

- 3.8+ recommended

Minimal for startup/dry-run:

- PyYAML

Runtime by phase:

- Phase 1: matplotlib, numpy, pandas, anndata, scanpy
- Phase 2: numpy, pandas
- Phase 3: matplotlib, numpy, pandas, seaborn, anndata

Install common stack:

```bash
pip install pyyaml numpy pandas matplotlib seaborn anndata scanpy
```


## 11) Configuration Guide By Intent

### `paths`

Defines input locations and output root/subfolders.
Use explicit paths when possible; fallback logic is robust but intentionally fails on ambiguous matches.


### `preprocessing`

Controls QC and dimensionality setup before marker calling.
`n_jobs` supports `-1` (all cores) or positive int.


### `marker_extraction`

Controls DE filtering, ranking behavior, stress handling, and fetal rescue.
Most impactful quality knobs:

- `min_pct_cells`
- `max_pct_rest`
- `adult_stress_mitigation_mode`
- `adult_stress_signature_genes`
- `fetal_rescue_min_margin`
- `fetal_rescue_min_score`
- `fetal_rescue_ambiguous_strategy`


### `hpa_validation`

Controls external validation thresholds and weighted final ranking.
Manual overrides are supported but must remain explicitly mapped.


### `overlap`

Controls harmonization, mapping seeds, acceptance criteria, and YAML schema version.
If you introduce new manual override labels upstream, ensure harmonization keys include them.


### `dataset_manifest` (optional)

Adds provenance metadata for adult/fetal inputs into `dataset_manifest.tsv`.


## 12) Tuning Playbook

If top markers look too broad/non-specific:

- Increase `min_pct_cells`.
- Decrease `max_pct_rest`.

If adult stress signatures dominate:

- Set `adult_stress_mitigation_mode` to `exclude_high_stress` or `regress_out`.
- Curate `adult_stress_signature_genes`.

If fetal rescue looks unstable:

- Increase `fetal_rescue_min_margin`.
- Increase `fetal_rescue_min_score`.
- Use `assign_label` for ambiguous handling if you want explicit ambiguity labeling.

If overlap seems unexpectedly sparse:

- Revisit `overlap.harmonization_map`.
- Add/adjust `overlap.manual_mapping_seed`.
- Inspect `cell_type_label_mapping.tsv` and overlap diagnostics.


## 13) Verification And Troubleshooting

Run:

```bash
python adult_fetal_ref_overlap.py --config <config_path> --verify
```

`--verify` checks:

- presence/non-empty status of expected outputs
- YAML parse structure
- absence of legacy notebook artifacts

Common issues:

- Missing package errors: install phase dependencies in active environment.
- Missing upstream artifacts for selected phases: run required upstream phase(s).
- Ambiguous input path resolution: use explicit absolute or config-relative path.
- No phases selected: occurs when all selected phases are skipped.


## 14) Practical Usage Patterns

First run:

```bash
python adult_fetal_ref_overlap.py --config <config_path> --dry-run
python adult_fetal_ref_overlap.py --config <config_path>
python adult_fetal_ref_overlap.py --config <config_path> --verify
```

Iterate overlap only after overlap-config edits:

```bash
python adult_fetal_ref_overlap.py --config <config_path> --only-phases 3
```

Regenerate everything in place:

```bash
python adult_fetal_ref_overlap.py --config <config_path> --no-resume
```


## 15) Summary

This pipeline is not just a scripted notebook export. It is a reproducible execution framework with explicit contracts, robust harmonization logic, confidence-aware rescue, external validation integration, and deterministic resume behavior.

Use it when you want biologically interpretable outputs that remain stable and auditable across reruns, team members, and evolving configs.
