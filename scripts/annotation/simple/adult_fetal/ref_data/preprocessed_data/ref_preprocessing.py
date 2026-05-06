"""
ref_preprocessing.py
--------------------
Preprocessing pipeline for adult and fetal kidney single-cell RNA-seq reference datasets.

What this does, in order:
  1. Load adult + fetal .h5ad references and standardize gene names.
  2. Audit raw cell-type labels and save the distribution to TSV.
  3. Remove known immune / non-kidney contaminant cell types.
  4. QC filter → normalize → log1p → HVG selection → scale → PCA → neighbors → UMAP.
  5. Strip high-ubiquity technical/housekeeping genes that dominate marker tests.
  6. Optionally drop cell types that have fewer than N cells (too small for robust DE).
  7. Save an intermediate preprocessed .h5ad for each stage.
  8. Rescue "generic label" fetal cells (e.g. "kidney cell") by scoring them against
     known marker panels and re-labelling the best-matching type above a score threshold.
  9. Harmonize nomenclature across the two stages (endothelial consolidation, fibroblast
     unification, etc.) and save a full harmonization mapping table.
 10. Write final harmonized .h5ad files ready for downstream marker extraction / annotation.

All tunable parameters live in config.yaml.  Paths are also there — edit them before
running rather than touching this file.
"""

import gc
import logging
import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

# ---------------------------------------------------------------------------
# Silence AnnData / scanpy chatter that doesn't add value in script mode.
# Actual pipeline progress is handled by our own logger below.
# ---------------------------------------------------------------------------
warnings.filterwarnings("ignore")
ad.settings.allow_write_nullable_strings = True


# ===========================================================================
#  Logging setup
# ===========================================================================

def setup_logging(log_level: str = "INFO", log_file: str | None = None) -> logging.Logger:
    """
    Build a logger that writes to stdout and, optionally, a file.

    We use a named logger instead of the root logger so downstream libraries
    (scanpy, anndata) that touch the root logger don't pollute our output.
    """
    logger = logging.getLogger("ref_preprocessing")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))

    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Always log to stdout so you can follow progress in real time.
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(fmt)
    logger.addHandler(stdout_handler)

    if log_file:
        file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)
        logger.info("Logging to file: %s", log_file)

    return logger


# ===========================================================================
#  Config loader
# ===========================================================================

def load_config(config_path: str = "config.yaml") -> dict:
    """
    Read the YAML config file.  Raises a clear error if the file is missing
    so users know immediately what to fix rather than hunting through tracebacks.
    """
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(
            f"Config file not found: {p.resolve()}\n"
            "Make sure config.yaml is in the same directory as this script, "
            "or pass a path via the CONFIG_PATH environment variable."
        )
    with open(p, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    return cfg


# ===========================================================================
#  Gene / cell-type helpers
# ===========================================================================

def ensure_gene_symbols(adata: ad.AnnData) -> ad.AnnData:
    """
    Swap var_names to the human-readable 'feature_name' column when it exists,
    then deduplicate.  Many 10x outputs store Ensembl IDs as var_names by default;
    this makes sure we work in gene-symbol space throughout.
    """
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).values
    adata.var_names_make_unique()
    return adata


def summarize_cell_types(adata: ad.AnnData, stage_label: str) -> pd.DataFrame:
    """
    Return a tidy DataFrame with cell counts and fractions per cell type.
    Used both for QC audits and for the final composition summaries.
    """
    counts = adata.obs["cell_type"].astype(str).value_counts().rename("n_cells")
    summary = counts.to_frame()
    summary["fraction"] = summary["n_cells"] / summary["n_cells"].sum()
    summary["stage"] = stage_label
    return summary.reset_index().rename(columns={"index": "cell_type"})


# ===========================================================================
#  Step 1 – Load references
# ===========================================================================

def load_references(cfg: dict, logger: logging.Logger) -> tuple[ad.AnnData, ad.AnnData]:
    """
    Load the adult and fetal .h5ad files, standardize gene names, and run a
    quick sanity check that the cell-type column is present in both objects.
    """
    adult_path = Path(cfg["paths"]["adult_h5ad"])
    fetal_path = Path(cfg["paths"]["fetal_h5ad"])

    for p in (adult_path, fetal_path):
        if not p.exists():
            raise FileNotFoundError(
                f"Reference file not found: {p}\n"
                "Update 'paths' in config.yaml to point to your .h5ad files."
            )

    logger.info("Loading adult reference: %s", adult_path)
    adata_adult = ensure_gene_symbols(sc.read_h5ad(str(adult_path)))
    logger.info("Adult raw: %s cells x %s genes", f"{adata_adult.n_obs:,}", f"{adata_adult.n_vars:,}")

    logger.info("Loading fetal reference: %s", fetal_path)
    adata_fetal = ensure_gene_symbols(sc.read_h5ad(str(fetal_path)))
    logger.info("Fetal raw: %s cells x %s genes", f"{adata_fetal.n_obs:,}", f"{adata_fetal.n_vars:,}")

    # Both datasets need a cell_type column for everything downstream to work.
    adult_ct_col = cfg["obs_columns"]["adult_celltype"]
    fetal_ct_col = cfg["obs_columns"]["fetal_celltype"]
    for obj, col, label in (
        (adata_adult, adult_ct_col, "adult"),
        (adata_fetal, fetal_ct_col, "fetal"),
    ):
        if col not in obj.obs.columns:
            raise KeyError(
                f"Column '{col}' not found in {label} .obs.\n"
                f"Available columns: {sorted(obj.obs.columns.tolist())}\n"
                "Fix 'obs_columns' in config.yaml."
            )
        # Normalise the column name to 'cell_type' so the rest of the pipeline
        # doesn't need to know the original name.
        if col != "cell_type":
            obj.obs["cell_type"] = obj.obs[col]

    logger.info("Adult obs columns: %s", sorted(adata_adult.obs.columns.tolist()))
    logger.info("Fetal obs columns: %s", sorted(adata_fetal.obs.columns.tolist()))

    return adata_adult, adata_fetal


# ===========================================================================
#  Step 2 – Label audit
# ===========================================================================

def audit_labels(
    adata_adult: ad.AnnData,
    adata_fetal: ad.AnnData,
    out_dir: Path,
    logger: logging.Logger,
) -> None:
    """
    Write a TSV showing raw cell-type distributions for both stages.
    Useful to inspect before committing to the removal list in config.yaml.
    """
    adult_audit = summarize_cell_types(adata_adult, "adult_raw")
    fetal_audit = summarize_cell_types(adata_fetal, "fetal_raw")
    audit = pd.concat([adult_audit, fetal_audit], ignore_index=True)

    audit_path = out_dir / "label_audit_raw.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)
    logger.info("Raw label audit saved: %s", audit_path)
    logger.debug("Adult label counts:\n%s", adult_audit.to_string(index=False))
    logger.debug("Fetal label counts:\n%s", fetal_audit.to_string(index=False))


# ===========================================================================
#  Step 3 – Remove contaminant / unwanted cell types
# ===========================================================================

def remove_unwanted_cell_types(
    adata: ad.AnnData,
    labels_to_remove: set,
    stage_label: str,
    logger: logging.Logger,
) -> ad.AnnData:
    """
    Drop cells whose (lowercased, stripped) cell_type label is in labels_to_remove.
    Logs which labels were actually found vs. listed but absent — helpful for
    catching typos in config.yaml.
    """
    ct_norm = adata.obs["cell_type"].astype(str).str.strip().str.lower()
    labels_norm = {lbl.lower().strip() for lbl in labels_to_remove}

    found = sorted(set(ct_norm.unique()) & labels_norm)
    missing = sorted(labels_norm - set(ct_norm.unique()))
    n_removed = int(ct_norm.isin(labels_norm).sum())

    logger.info("[%s] Labels to remove that were found: %s", stage_label, found)
    if missing:
        logger.warning(
            "[%s] Labels listed for removal but NOT found in data: %s — check config.yaml for typos.",
            stage_label, missing,
        )
    logger.info("[%s] Removing %s cells across %s label(s).", stage_label, f"{n_removed:,}", len(found))

    adata = adata[~ct_norm.isin(labels_norm)].copy()
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    return adata


# ===========================================================================
#  Step 4 – Preprocessing pipeline (QC → norm → HVG → PCA → UMAP)
# ===========================================================================

def preprocess_for_markers(
    adata: ad.AnnData,
    label: str,
    cfg: dict,
    logger: logging.Logger,
) -> ad.AnnData:
    """
    Standard scRNA-seq preprocessing:
      - Cell/gene QC filters
      - Library-size normalization + log1p
      - Highly variable gene selection (batch-aware if dataset_id column exists)
      - Scaling, PCA, k-NN graph, UMAP

    We stash the normalized-but-unscaled data in adata.raw so that marker
    tests later can work with values closer to the original expression
    distribution (scaling distorts fold-change estimates).
    """
    qc = cfg["qc"]
    hvg_cfg = cfg["hvg"]
    pca_cfg = cfg["pca"]

    adata = adata.copy()
    logger.info("[%s] Starting preprocessing: %s cells x %s genes", label, f"{adata.n_obs:,}", f"{adata.n_vars:,}")

    # --- QC filters -----------------------------------------------------------
    sc.pp.filter_cells(adata, min_genes=qc["min_genes_per_cell"])
    sc.pp.filter_genes(adata, min_cells=qc["min_cells_per_gene"])
    logger.info("[%s] After QC: %s cells x %s genes", label, f"{adata.n_obs:,}", f"{adata.n_vars:,}")

    if adata.n_obs == 0:
        raise ValueError(f"[{label}] No cells remain after QC filtering — check your thresholds in config.yaml.")

    # --- Normalization + log transform ----------------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Snapshot before scaling — marker tests should use this, not scaled values.
    adata.raw = adata

    # --- Highly variable genes ------------------------------------------------
    batch_key = "dataset_id" if "dataset_id" in adata.obs.columns else None
    if batch_key:
        logger.info("[%s] HVG selection using batch_key='%s'", label, batch_key)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg_cfg["n_top_genes"], batch_key=batch_key)
    adata = adata[:, adata.var["highly_variable"]].copy()
    logger.info("[%s] After HVG selection: %s genes retained", label, f"{adata.n_vars:,}")

    # --- Scale → PCA → neighbors → UMAP --------------------------------------
    sc.pp.scale(adata, max_value=pca_cfg["scale_max_value"])
    sc.tl.pca(
        adata,
        svd_solver=pca_cfg["svd_solver"],
        n_comps=pca_cfg["n_comps"],
        random_state=pca_cfg["random_state"],
    )

    # Use fewer PCs if the data doesn't support the full requested number.
    n_pcs = min(pca_cfg["n_pcs_neighbors"], adata.obsm["X_pca"].shape[1])
    if n_pcs < 2:
        raise ValueError(f"[{label}] Only {n_pcs} PCs available — too few to build a meaningful neighbors graph.")

    sc.pp.neighbors(adata, n_pcs=n_pcs)
    sc.tl.umap(adata, random_state=pca_cfg["random_state"])

    logger.info(
        "[%s] Preprocessing done: %s cells x %s HVGs, UMAP embedded.",
        label, f"{adata.n_obs:,}", f"{adata.n_vars:,}",
    )
    return adata


# ===========================================================================
#  Step 5 – Remove technical / housekeeping genes
# ===========================================================================

def filter_noise_genes(
    adata: ad.AnnData,
    label: str,
    cfg: dict,
    logger: logging.Logger,
) -> ad.AnnData:
    """
    Strip genes that dominate one-vs-rest marker tests for the wrong reasons:
    mitochondrial genes, ribosomal proteins, HLA loci, housekeeping genes, etc.

    We filter both .X (the HVG-subset used for dimensionality reduction) and
    .raw (the full normalized matrix used for DE), so neither layer leaks
    technical signal into marker extraction.
    """
    noise_cfg = cfg["noise_genes"]
    noise_prefixes = tuple(noise_cfg["prefixes"])
    exact_noise = set(noise_cfg["exact"])

    def _noise_mask(var_names):
        return ~(
            var_names.str.startswith(noise_prefixes)
            | var_names.isin(exact_noise)
        )

    before_x = adata.n_vars
    mask_x = _noise_mask(adata.var_names)
    adata = adata[:, mask_x].copy()

    if adata.raw is not None:
        raw_ad = adata.raw.to_adata()
        before_raw = raw_ad.n_vars
        mask_raw = _noise_mask(raw_ad.var_names)
        adata.raw = raw_ad[:, mask_raw]
        logger.info(
            "[%s] Noise gene filter — .X: %s → %s, .raw: %s → %s",
            label, before_x, adata.n_vars, before_raw, adata.raw.n_vars,
        )
    else:
        logger.info("[%s] Noise gene filter — .X: %s → %s (no .raw found)", label, before_x, adata.n_vars)

    return adata


# ===========================================================================
#  Step 6 – Minimum cells-per-type filter
# ===========================================================================

def enforce_min_cells_per_type(
    adata: ad.AnnData,
    min_cells: int,
    logger: logging.Logger,
) -> ad.AnnData:
    """
    Drop cell types with fewer than min_cells observations.  One-vs-rest
    Wilcoxon tests become unreliable below ~20 cells, so we remove these
    rather than letting them produce noisy marker lists.
    """
    adata = adata.copy()
    counts = adata.obs["cell_type"].value_counts()
    keep = counts[counts >= min_cells].index
    removed = counts[counts < min_cells]

    if len(removed) > 0:
        logger.warning(
            "Removed %s low-count cell type(s) (< %s cells): %s",
            len(removed), min_cells, removed.index.tolist(),
        )
    else:
        logger.info("All cell types pass the minimum-cells threshold (%s).", min_cells)

    adata = adata[adata.obs["cell_type"].isin(keep)].copy()
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    return adata


# ===========================================================================
#  Step 7 – Fetal generic-label rescue
# ===========================================================================

def rescue_generic_fetal_cells(
    fetal: ad.AnnData,
    cfg: dict,
    logger: logging.Logger,
) -> ad.AnnData:
    """
    Some fetal cells carry uninformative labels like "kidney cell" because
    the original study couldn't resolve them at annotation time.  We score
    each such cell against curated marker gene sets and re-label it as the
    best-matching specific type — but only if the score clears a minimum
    threshold (to avoid forcing a label onto truly ambiguous cells).

    Cells that score below the threshold are simply dropped rather than
    passing garbage labels downstream.
    """
    rescue_cfg = cfg["rescue"]
    rescue_markers: dict = rescue_cfg["markers"]
    generic_labels: list = rescue_cfg["generic_labels"]
    score_threshold: float = rescue_cfg["score_threshold"]

    generic_mask = fetal.obs["cell_type"].isin(generic_labels)
    fetal_to_rescue = fetal[generic_mask].copy()
    fetal_specific = fetal[~generic_mask].copy()

    logger.info(
        "Rescue: targeting %s cells with generic labels: %s",
        f"{fetal_to_rescue.n_obs:,}", generic_labels,
    )

    if fetal_to_rescue.n_obs == 0:
        logger.info("No cells with generic labels found — skipping rescue step.")
        return fetal

    # Score each cell against each rescue marker panel.
    rescue_scores = pd.DataFrame(index=fetal_to_rescue.obs_names)

    for target_type, markers in rescue_markers.items():
        available = [m for m in markers if m in fetal_to_rescue.var_names]
        if not available:
            logger.warning(
                "Rescue: no markers from %s found in fetal var_names — assigning score 0.",
                target_type,
            )
            rescue_scores[target_type] = 0.0
            continue
        missing_markers = set(markers) - set(available)
        if missing_markers:
            logger.debug("Rescue [%s]: markers not in data (skipped): %s", target_type, missing_markers)

        # Mean expression across available markers as a simple, interpretable score.
        expr = fetal_to_rescue[:, available].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        rescue_scores[target_type] = np.asarray(expr).mean(axis=1)

    fetal_to_rescue.obs["rescue_type"] = rescue_scores.idxmax(axis=1)
    fetal_to_rescue.obs["rescue_score"] = rescue_scores.max(axis=1)

    rescued_mask = fetal_to_rescue.obs["rescue_score"] > score_threshold
    fetal_rescued = fetal_to_rescue[rescued_mask].copy()
    fetal_rescued.obs["cell_type"] = fetal_rescued.obs["rescue_type"]

    n_dropped = fetal_to_rescue.n_obs - fetal_rescued.n_obs
    logger.info(
        "Rescue results: %s rescued, %s dropped (score ≤ %.3f).",
        f"{fetal_rescued.n_obs:,}", f"{n_dropped:,}", score_threshold,
    )
    for ct in sorted(fetal_rescued.obs["cell_type"].unique()):
        n = (fetal_rescued.obs["cell_type"] == ct).sum()
        logger.info("  └─ %s: %s cells", ct, f"{n:,}")

    # Merge rescued cells back with those that already had specific labels.
    fetal_final = ad.concat([fetal_specific, fetal_rescued], axis=0)
    fetal_final.obs["cell_type"] = fetal_final.obs["cell_type"].astype("category")

    logger.info(
        "Fetal after rescue: %s cells (%s specific + %s rescued).",
        f"{fetal_final.n_obs:,}", f"{fetal_specific.n_obs:,}", f"{fetal_rescued.n_obs:,}",
    )
    return fetal_final


# ===========================================================================
#  Step 8 – Cross-stage nomenclature harmonization
# ===========================================================================

def harmonize_cell_types(
    adata: ad.AnnData,
    harmonization_map: dict,
    stage_name: str,
    logger: logging.Logger,
) -> ad.AnnData:
    """
    Apply a label→label mapping to collapse redundant or stage-specific
    nomenclature into a shared vocabulary.

    The original label is kept in 'cell_type'; the harmonized label goes
    into 'cell_type_harmonized'.  Unmapped labels pass through unchanged so
    you don't silently lose anything.
    """
    adata = adata.copy()
    before_types = adata.obs["cell_type"].nunique()

    adata.obs["cell_type_harmonized"] = (
        adata.obs["cell_type"]
        .map(lambda x: harmonization_map.get(x, x))  # identity for unmapped
        .astype("category")
    )

    after_types = adata.obs["cell_type_harmonized"].nunique()
    logger.info(
        "[%s] Harmonization: %s → %s cell types (%s label(s) consolidated).",
        stage_name, before_types, after_types, before_types - after_types,
    )

    # Log endothelial-specific stats since that's the main merge here.
    endo_sources = [k for k, v in harmonization_map.items() if v == "endothelial cell"]
    if endo_sources:
        before_counts = adata.obs["cell_type"].value_counts().to_dict()
        total_endo = sum(before_counts.get(lbl, 0) for lbl in endo_sources)
        logger.info(
            "[%s] Endothelial consolidation: %s subtypes → 1, %s cells total.",
            stage_name, len(endo_sources), f"{total_endo:,}",
        )

    return adata


def build_harmonization_table(
    harmonization_map: dict,
    adata_adult: ad.AnnData,
    adata_fetal: ad.AnnData,
) -> pd.DataFrame:
    """Build a documentation table showing how many adult/fetal cells each mapping affects."""
    records = []
    for source, target in harmonization_map.items():
        consol_type = "endothelial_unification" if target == "endothelial cell" else "nomenclature_harmonization"
        records.append({
            "source_label": source,
            "target_label": target,
            "consolidation_type": consol_type,
            "adult_cells_affected": int((adata_adult.obs["cell_type"] == source).sum()),
            "fetal_cells_affected": int((adata_fetal.obs["cell_type"] == source).sum()),
        })
    return pd.DataFrame(records)


# ===========================================================================
#  Main pipeline
# ===========================================================================

def run_pipeline(config_path: str = "config.yaml") -> None:
    """
    Orchestrate all preprocessing steps from loading raw .h5ad files to
    writing final harmonized references.
    """
    cfg = load_config(config_path)

    # Resolve output directories; create them if they don't exist yet.
    out_dir = Path(cfg["paths"]["out_dir"])
    preprocessed_dir = Path(cfg["paths"]["preprocessed_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    preprocessed_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(
        log_level=cfg.get("logging", {}).get("level", "INFO"),
        log_file=cfg.get("logging", {}).get("file", None),
    )

    # Let scanpy use all available cores for neighbor computation etc.
    sc.settings.verbosity = 0  # We handle our own progress output.
    sc.settings.n_jobs = cfg.get("scanpy", {}).get("n_jobs", -1)

    logger.info("=" * 70)
    logger.info("Reference preprocessing pipeline — starting")
    logger.info("Config: %s", Path(config_path).resolve())
    logger.info("Output dir: %s", out_dir.resolve())
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # 1. Load
    # ------------------------------------------------------------------
    logger.info("--- Step 1: Loading reference datasets ---")
    try:
        adata_adult, adata_fetal = load_references(cfg, logger)
    except (FileNotFoundError, KeyError) as exc:
        logger.error("Failed to load reference data: %s", exc)
        sys.exit(1)

    # ------------------------------------------------------------------
    # 2. Audit raw labels
    # ------------------------------------------------------------------
    logger.info("--- Step 2: Auditing raw cell-type labels ---")
    audit_labels(adata_adult, adata_fetal, out_dir, logger)

    # ------------------------------------------------------------------
    # 3. Remove contaminant cell types
    # ------------------------------------------------------------------
    logger.info("--- Step 3: Removing unwanted / contaminant cell types ---")
    removal_cfg = cfg["cell_type_removal"]
    adata_adult = remove_unwanted_cell_types(
        adata_adult, set(removal_cfg["adult"]), "adult", logger
    )
    adata_fetal = remove_unwanted_cell_types(
        adata_fetal, set(removal_cfg["fetal"]), "fetal", logger
    )

    # Save the post-removal audit so you can see what's left.
    audit_labels(adata_adult, adata_fetal, out_dir, logger)
    post_audit = pd.concat(
        [summarize_cell_types(adata_adult, "adult_raw"), summarize_cell_types(adata_fetal, "fetal_raw")],
        ignore_index=True,
    )
    post_audit.to_csv(out_dir / "label_audit_raw_after_removal.tsv", sep="\t", index=False)
    logger.info("Post-removal label audit saved.")

    # ------------------------------------------------------------------
    # 4. Preprocess (QC → norm → HVG → PCA → UMAP)
    # ------------------------------------------------------------------
    logger.info("--- Step 4: Preprocessing (QC, normalization, dimensionality reduction) ---")
    try:
        adata_adult = preprocess_for_markers(adata_adult, "adult", cfg, logger)
        adata_fetal = preprocess_for_markers(adata_fetal, "fetal", cfg, logger)
    except ValueError as exc:
        logger.error("Preprocessing failed: %s", exc)
        sys.exit(1)

    # ------------------------------------------------------------------
    # 5. Remove noise / housekeeping genes
    # ------------------------------------------------------------------
    logger.info("--- Step 5: Filtering technical / noise genes ---")
    adata_adult = filter_noise_genes(adata_adult, "adult", cfg, logger)
    adata_fetal = filter_noise_genes(adata_fetal, "fetal", cfg, logger)

    # ------------------------------------------------------------------
    # 6. Minimum cells-per-type filter
    # ------------------------------------------------------------------
    filter_cfg = cfg["cell_type_filter"]
    if filter_cfg["apply"]:
        logger.info("--- Step 6: Enforcing minimum cells per cell type (%s) ---", filter_cfg["min_cells"])
        adata_adult = enforce_min_cells_per_type(adata_adult, filter_cfg["min_cells"], logger)
        adata_fetal = enforce_min_cells_per_type(adata_fetal, filter_cfg["min_cells"], logger)
    else:
        logger.info("--- Step 6: Minimum cells-per-type filter SKIPPED (set apply: true in config to enable) ---")

    # Save intermediate preprocessed .h5ad (before rescue / harmonization).
    adult_pre_out = preprocessed_dir / "adult_reference_preprocessed.h5ad"
    fetal_pre_out = preprocessed_dir / "fetal_reference_preprocessed.h5ad"
    logger.info("Saving intermediate preprocessed references ...")
    adata_adult.write_h5ad(adult_pre_out)
    adata_fetal.write_h5ad(fetal_pre_out)
    logger.info("Adult preprocessed → %s", adult_pre_out)
    logger.info("Fetal preprocessed → %s", fetal_pre_out)

    # Save composition summary at this stage.
    composition_pre = pd.concat(
        [summarize_reference(adata_adult, "adult_preprocessed", logger),
         summarize_reference(adata_fetal, "fetal_preprocessed", logger)],
        ignore_index=True,
    )
    composition_pre.to_csv(
        out_dir / "reference_cell_type_summary_preprocessed.tsv", sep="\t", index=False
    )

    # ------------------------------------------------------------------
    # 7. Rescue generic fetal labels
    # ------------------------------------------------------------------
    logger.info("--- Step 7: Rescuing fetal cells with generic labels ---")
    adata_fetal = rescue_generic_fetal_cells(adata_fetal, cfg, logger)

    # Free memory we no longer need before the harmonization step.
    gc.collect()

    # ------------------------------------------------------------------
    # 8. Harmonize nomenclature across stages
    # ------------------------------------------------------------------
    logger.info("--- Step 8: Harmonizing cell-type nomenclature across adult and fetal ---")
    harm_map: dict = cfg["harmonization_map"]

    adata_adult_harm = harmonize_cell_types(adata_adult, harm_map, "adult", logger)
    adata_fetal_harm = harmonize_cell_types(adata_fetal, harm_map, "fetal", logger)

    # Save the mapping table — useful for methods sections and reproducibility.
    harm_table = build_harmonization_table(harm_map, adata_adult, adata_fetal)
    harm_table_path = out_dir / "cell_type_harmonization_mapping.tsv"
    harm_table.to_csv(harm_table_path, sep="\t", index=False)
    logger.info("Harmonization mapping table saved: %s", harm_table_path)

    # Final composition summary.
    composition_harm = pd.concat(
        [summarize_reference(adata_adult_harm, "adult_harmonized", logger),
         summarize_reference(adata_fetal_harm, "fetal_harmonized", logger)],
        ignore_index=True,
    )
    composition_harm.to_csv(
        out_dir / "reference_cell_type_summary_harmonized.tsv", sep="\t", index=False
    )

    # ------------------------------------------------------------------
    # 9. Write final harmonized references
    # ------------------------------------------------------------------
    adult_final_out = preprocessed_dir / "adult_preprocessed.h5ad"
    fetal_final_out = preprocessed_dir / "fetal_preprocessed.h5ad"

    logger.info("Saving final harmonized references ...")
    adata_adult_harm.write_h5ad(adult_final_out)
    adata_fetal_harm.write_h5ad(fetal_final_out)
    logger.info("Adult final → %s", adult_final_out)
    logger.info("Fetal final → %s", fetal_final_out)

    logger.info("=" * 70)
    logger.info("Pipeline complete.  Outputs in: %s", out_dir.resolve())
    logger.info("=" * 70)


def summarize_reference(adata: ad.AnnData, label: str, logger: logging.Logger) -> pd.DataFrame:
    """
    Print and return a per-cell-type count/fraction table.
    Re-used both for intermediate and final summaries.
    """
    ct_col = "cell_type_harmonized" if "cell_type_harmonized" in adata.obs.columns else "cell_type"
    counts = adata.obs[ct_col].astype(str).value_counts().rename("n_cells").to_frame()
    counts["fraction"] = counts["n_cells"] / counts["n_cells"].sum()
    counts["stage"] = label
    logger.info(
        "[%s] %s cells across %s cell types.", label, f"{adata.n_obs:,}", counts.shape[0]
    )
    return counts.reset_index().rename(columns={"index": "cell_type"})


# ===========================================================================
#  Entry point
# ===========================================================================

if __name__ == "__main__":
    import os

    # Allow overriding the config path via an environment variable so the same
    # script can be used with different configs in CI / batch jobs without
    # touching the code.
    config_path = os.environ.get("CONFIG_PATH", "config.yaml")
    run_pipeline(config_path=config_path)
