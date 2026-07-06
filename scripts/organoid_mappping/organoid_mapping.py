#!/usr/bin/env python3
"""
Organoid Marker Extraction Pipeline (scanpy / h5ad edition)
=============================================================

Python/scanpy port of the original Seurat/RDS notebook. Same logic, new
stack: load each dataset's .h5ad, cluster it, pull out Wilcoxon markers,
score and curate them, then (optionally) map the resulting clusters onto
adult/fetal reference cell types.

Clustering resolution can either be fixed per dataset (default, from
config.yaml) or found automatically via a silhouette + marker-coherence
sweep (opt-in, slower — see clustering.dynamic_resolution in the config).

Everything tunable lives in config.yaml — run with:

    python organoid_mapping.py --config organoid_mapping_config.yaml
"""

from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from sklearn.metrics import silhouette_score

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("organoid_pipeline")


# =============================================================================
# Config
# =============================================================================

def load_config(config_path: str | Path) -> dict:
    with open(config_path, "r") as fh:
        cfg = yaml.safe_load(fh)
    return cfg


@dataclass
class NoiseGeneFilter:
    prefix_pattern: re.Pattern
    exact: set[str]
    exclude_ensembl: bool

    @classmethod
    def from_config(cls, cfg: dict) -> "NoiseGeneFilter":
        ng = cfg["noise_genes"]
        prefixes = "|".join(re.escape(p) for p in ng["prefixes"])
        pattern = re.compile(f"^({prefixes})")
        return cls(
            prefix_pattern=pattern,
            exact=set(ng["exact"]),
            exclude_ensembl=ng.get("exclude_ensembl_prefixed", True),
        )

    def clean(self, genes: Iterable[str]) -> list[str]:
        out = []
        for g in genes:
            if self.exclude_ensembl and g.startswith("ENSG"):
                continue
            if g in self.exact:
                continue
            if self.prefix_pattern.match(g):
                continue
            out.append(g)
        return out


# =============================================================================
# 1 . Discovery & loading
# =============================================================================

def discover_h5ad_files(input_dir: str | Path) -> list[Path]:
    """Top-level .h5ad files only (subfolders ignored), sorted for determinism."""
    input_dir = Path(input_dir)
    files = sorted(p for p in input_dir.glob("*.h5ad") if p.is_file())
    log.info("Found %d .h5ad files in %s", len(files), input_dir)
    return files


def dataset_name_from_path(path: Path) -> str:
    return path.stem  # filename without .h5ad


def load_adata(path: Path, counts_layer: str = "counts") -> sc.AnnData:
    """
    Load an h5ad file and make sure raw counts are preserved.

    Assumes adata.X on disk holds raw (unnormalised) counts, which is the
    standard convention for QC'd h5ad exports. If a `counts` layer already
    exists it's trusted as-is; otherwise it's created from adata.X.
    """
    adata = sc.read_h5ad(path)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if counts_layer not in adata.layers:
        adata.layers[counts_layer] = adata.X.copy()

    # Raw counts are preserved in `adata.layers[counts_layer]`; avoid duplicating them in `adata.raw`.

    adata.obs["dataset_name"] = dataset_name_from_path(path)
    return adata


# =============================================================================
# 2 . Preprocessing: Normalise -> HVG -> Scale -> PCA -> Neighbors
# =============================================================================

def preprocess_adata(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    pp = cfg["preprocessing"]

    # Normalise
    sc.pp.normalize_total(adata, target_sum=pp["normalize_target_sum"])
    sc.pp.log1p(adata)

    # Highly variable genes
    sc.pp.highly_variable_genes(
        adata, n_top_genes=pp["n_hvg"], flavor="seurat", subset=False
    )

    # Scale HVGs only (mirrors Seurat's ScaleData on VariableFeatures)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=pp["scale_max_value"])

    # PCA
    sc.tl.pca(
        adata_hvg,
        n_comps=pp["n_pcs"],
        svd_solver="arpack",
        random_state=pp["random_state"],
    )
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    # kNN graph
    sc.pp.neighbors(
        adata,
        n_neighbors=pp["n_neighbors"],
        n_pcs=pp["n_pcs"],
        random_state=pp["random_state"],
    )
    return adata


# =============================================================================
# 3 . Clustering at a FIXED resolution (no dynamic sweep)
# =============================================================================

def cluster_adata(adata: sc.AnnData, resolution: float, cfg: dict) -> sc.AnnData:
    cl = cfg["clustering"]
    key = cl["cluster_key"]

    if cl["algorithm"] == "leiden":
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=key,
            flavor=cl.get("flavor", "igraph"),
            n_iterations=cl.get("n_iterations", 2),
            random_state=cfg["preprocessing"]["random_state"],
        )
    elif cl["algorithm"] == "louvain":
        sc.tl.louvain(
            adata,
            resolution=resolution,
            key_added=key,
            random_state=cfg["preprocessing"]["random_state"],
        )
    else:
        raise ValueError(f"Unknown clustering algorithm: {cl['algorithm']}")

    n_clusters = adata.obs[key].nunique()
    log.info("    clustered at resolution=%.2f -> %d clusters", resolution, n_clusters)
    return adata


# =============================================================================
# 3b . Dynamic resolution clustering (optional, off by default)
# =============================================================================
# Instead of a single pre-chosen resolution, this sweeps a range of
# resolutions, scores each with silhouette width, shortlists the best few,
# picks the one with the fewest "duplicate" clusters (near-identical marker
# profiles), then merges duplicate and micro clusters. It's the same
# strategy as the original R sweep, just driven off scanpy/Leiden objects.
#
# Turned on via clustering.dynamic_resolution.enabled: true in config.yaml.
# NB: this is a lot slower than fixed-resolution clustering — it clusters
# the dataset once per candidate resolution, plus a Wilcoxon pass per merge
# iteration. Only turn it on if you don't already know a good resolution.
# =============================================================================

def compute_silhouette_for_resolution(
    adata: sc.AnnData, resolution: float, cfg: dict
) -> tuple[float, pd.Series]:
    """
    Cluster at the given resolution and score the result with mean
    silhouette width (computed on a subsample of cells, for speed).
    Returns (score, cluster_labels). Score is -1 if the resolution
    collapses everything into a single cluster.
    """
    cl = cfg["clustering"]
    dr = cl["dynamic_resolution"]
    tmp_key = "_dyn_res_tmp"

    if cl["algorithm"] == "leiden":
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=tmp_key,
            flavor=cl.get("flavor", "igraph"),
            n_iterations=cl.get("n_iterations", 2),
            random_state=cfg["preprocessing"]["random_state"],
        )
    else:
        sc.tl.louvain(
            adata,
            resolution=resolution,
            key_added=tmp_key,
            random_state=cfg["preprocessing"]["random_state"],
        )

    labels = adata.obs[tmp_key].copy()
    if labels.nunique() < 2:
        return -1.0, labels

    n_cells = min(adata.n_obs, dr.get("silhouette_max_cells", 5000))
    rng = np.random.default_rng(cfg["preprocessing"]["random_state"])
    idx = rng.choice(adata.n_obs, size=n_cells, replace=False)

    pca = adata.obsm["X_pca"][idx]
    score = silhouette_score(pca, labels.values[idx])
    return float(score), labels


def find_duplicate_cluster_pairs(
    adata: sc.AnnData, cluster_key: str, cfg: dict
) -> list[tuple[str, str]]:
    """
    Quick-and-loose marker pass per cluster (relaxed thresholds, not the
    curated Phase 1 pipeline), then flag pairs of clusters whose top marker
    genes overlap by at least `merge_overlap` — a sign they're really the
    same population split in two.
    """
    dr = cfg["clustering"]["dynamic_resolution"]
    top_n = dr.get("merge_top_n", 15)
    overlap_thresh = dr.get("merge_overlap", 0.6)
    min_pct = dr.get("merge_min_pct", 0.1)
    min_logfc = dr.get("merge_min_logfc", 0.25)

    try:
        adata_tmp = adata.copy()
        sc.tl.rank_genes_groups(
            adata_tmp, groupby=cluster_key, method="wilcoxon", use_raw=False, pts=True
        )
        de = sc.get.rank_genes_groups_df(adata_tmp, group=None)
        de = de.rename(columns={"names": "gene", "group": "cluster", "logfoldchanges": "logFC"})
        de = de[(de["logFC"] >= min_logfc) & (de["pct_nz_group"] >= min_pct)]

        gene_sets: dict[str, set[str]] = {}
        for cluster, grp in de.groupby("cluster"):
            top = grp.sort_values("logFC", ascending=False).head(top_n)
            gene_sets[str(cluster)] = set(top["gene"])
    except Exception:
        return []

    duplicates = []
    for c1, c2 in combinations(sorted(gene_sets.keys()), 2):
        a, b = gene_sets.get(c1, set()), gene_sets.get(c2, set())
        if not a or not b:
            continue
        if len(a & b) / top_n >= overlap_thresh:
            duplicates.append((c1, c2))
    return duplicates


def merge_duplicate_clusters(adata: sc.AnnData, cluster_key: str, cfg: dict) -> pd.Series:
    """Iteratively fold together clusters that keep showing up as duplicates."""
    dr = cfg["clustering"]["dynamic_resolution"]
    max_iterations = dr.get("merge_max_iterations", 10)

    working_key = "_merge_tmp"
    adata.obs[working_key] = adata.obs[cluster_key].astype(str)

    for _ in range(max_iterations):
        dups = find_duplicate_cluster_pairs(adata, working_key, cfg)
        if not dups:
            break
        c1, c2 = dups[0]
        target, source = sorted((c1, c2), key=int)
        adata.obs[working_key] = adata.obs[working_key].replace(source, target)

    old_levels = sorted(adata.obs[working_key].unique(), key=int)
    remap = {old: str(i) for i, old in enumerate(old_levels)}
    merged = adata.obs[working_key].map(remap)
    adata.obs.drop(columns=[working_key], inplace=True)
    return merged


def merge_micro_clusters(adata: sc.AnnData, cluster_key: str, cfg: dict) -> pd.Series:
    """Absorb clusters smaller than `min_cluster_frac` into their dominant kNN neighbour."""
    dr = cfg["clustering"]["dynamic_resolution"]
    min_frac = dr.get("min_cluster_frac", 0.005)

    labels = adata.obs[cluster_key].astype(str).to_numpy().copy()
    threshold = int(np.floor(min_frac * len(labels)))
    counts = pd.Series(labels).value_counts()
    tiny = counts[counts < threshold].index.tolist()
    if not tiny:
        return pd.Series(labels, index=adata.obs_names)

    knn = adata.obsp["distances"].tocsr()  # binary-ish kNN graph from sc.pp.neighbors

    for cl in tiny:
        cl_idx = np.where(labels == cl)[0]
        neighbour_labels = []
        for i in cl_idx:
            nb_idx = knn.getrow(i).indices
            neighbour_labels.extend(labels[j] for j in nb_idx if labels[j] != cl)
        if neighbour_labels:
            dominant = pd.Series(neighbour_labels).value_counts().idxmax()
            labels[cl_idx] = dominant

    return pd.Series(labels, index=adata.obs_names)


def cluster_adata_dynamic(adata: sc.AnnData, cfg: dict, dataset_name: str) -> sc.AnnData:
    """
    Sweep clustering.dynamic_resolution.resolution_min..max, shortlist the
    top-scoring resolutions by silhouette, pick the one with the fewest
    duplicate-marker cluster pairs, then merge duplicate and micro clusters.
    """
    cl = cfg["clustering"]
    dr = cl["dynamic_resolution"]
    key = cl["cluster_key"]

    resolutions = np.round(
        np.arange(dr["resolution_min"], dr["resolution_max"] + 1e-9, dr["resolution_step"]), 2
    )
    log.info("  [%s] sweeping %d candidate resolutions", dataset_name, len(resolutions))

    scored = []
    for res in resolutions:
        score, labels = compute_silhouette_for_resolution(adata, float(res), cfg)
        log.info("    resolution=%.2f -> silhouette=%.4f", res, score)
        scored.append((float(res), score, labels))

    scored.sort(key=lambda item: item[1], reverse=True)
    shortlist = scored[: dr.get("top_k_by_silhouette", 3)]

    best_res, best_labels, best_n_dups = None, None, None
    for res, score, labels in shortlist:
        adata.obs[key] = labels
        n_dups = len(find_duplicate_cluster_pairs(adata, key, cfg))
        log.info(
            "    candidate resolution=%.2f silhouette=%.4f duplicate_pairs=%d",
            res, score, n_dups,
        )
        if best_n_dups is None or n_dups < best_n_dups:
            best_res, best_labels, best_n_dups = res, labels, n_dups

    log.info(
        "  [%s] chosen resolution=%.2f (fewest duplicate marker pairs)",
        dataset_name, best_res,
    )
    adata.obs[key] = best_labels
    adata.obs.drop(columns=["_dyn_res_tmp"], inplace=True, errors="ignore")

    adata.obs[key] = merge_duplicate_clusters(adata, key, cfg)
    adata.obs[key] = merge_micro_clusters(adata, key, cfg).astype("category")

    log.info(
        "  [%s] dynamic clustering final: %d clusters", dataset_name, adata.obs[key].nunique()
    )
    return adata


# =============================================================================
# 4 . Wilcoxon marker extraction (pct_in / pct_out from RAW counts)
# =============================================================================

def compute_pct_in_out(
    adata: sc.AnnData, cluster_key: str, genes: list[str], counts_layer: str
) -> pd.DataFrame:
    """
    For every gene x cluster pair, compute:
        pct_in  = fraction of cells IN the cluster with count > 0
        pct_out = fraction of cells OUTSIDE the cluster with count > 0

    Computed from raw counts, mirroring Seurat's slot='counts' percentages.
    """
    counts = adata[:, genes].layers[counts_layer]
    if not isinstance(counts, np.ndarray):
        counts = counts.toarray()
    expressed = counts > 0  # cells x genes boolean

    clusters = adata.obs[cluster_key].astype(str).values
    unique_clusters = pd.unique(clusters)

    rows = []
    for cl in unique_clusters:
        in_mask = clusters == cl
        pct_in = expressed[in_mask].mean(axis=0)
        pct_out = expressed[~in_mask].mean(axis=0)
        rows.append(
            pd.DataFrame(
                {
                    "gene": genes,
                    "cluster": cl,
                    "pct_in": pct_in,
                    "pct_out": pct_out,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def extract_wilcoxon_markers(
    adata: sc.AnnData,
    dataset_name: str,
    cluster_key: str,
    noise_filter: NoiseGeneFilter,
    cfg: dict,
) -> pd.DataFrame:
    """
    Run Wilcoxon rank_genes_groups on log-normalised data (adata.X) for
    logFC/significance, and separately compute pct_in/pct_out from raw
    counts, matching the original pipeline's rationale for using counts to
    calculate percentages while ranking on normalised expression.
    """
    me_cfg = cfg["marker_extraction"]
    clean_genes = noise_filter.clean(adata.var_names)
    if len(clean_genes) == 0:
        log.warning("    [%s] no genes left after noise filtering, skipping", dataset_name)
        return pd.DataFrame()

    adata_sub = adata[:, clean_genes].copy()

    sc.tl.rank_genes_groups(
        adata_sub,
        groupby=cluster_key,
        method=me_cfg["method"],
        corr_method=me_cfg["corr_method"],
        use_raw=False,
    )

    de = sc.get.rank_genes_groups_df(adata_sub, group=None)
    de = de.rename(
        columns={
            "names": "gene",
            "group": "cluster",
            "logfoldchanges": "logFC",
            "pvals": "pval",
            "pvals_adj": "pval_adj",
        }
    )[["gene", "cluster", "logFC", "pval", "pval_adj", "scores"]]

    if me_cfg.get("only_positive", True):
        de = de[de["logFC"] > 0]

    pct_df = compute_pct_in_out(
        adata, cluster_key, clean_genes, me_cfg["counts_layer"]
    )

    merged = de.merge(pct_df, on=["gene", "cluster"], how="left")
    merged.insert(0, "dataset", dataset_name)
    return merged


# =============================================================================
# 5 . Baseline filters + stringent score + elbow detection
# =============================================================================

def apply_baseline_filters(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    mf = cfg["marker_filtering"]
    return df[
        (df["pct_in"] >= mf["min_pct_in"]) & (df["logFC"] >= mf["min_logfc"])
    ].copy()


def compute_stringent_score(df: pd.DataFrame, tolerance: float) -> pd.DataFrame:
    df = df.copy()
    df["stringent_score"] = df["logFC"] * (df["pct_in"] ** 2) / (df["pct_out"] + tolerance)
    return df


def detect_elbow(scores: np.ndarray, max_show: int, min_kept: int) -> int:
    """
    Given scores sorted descending, detect the elbow using normalised
    distance (same 0-indexed geometry as the original R implementation).
    Returns the number of genes to keep (>= min_kept).
    """
    y = np.asarray(scores[:max_show], dtype=float)
    n = len(y)
    if n <= 2:
        return n

    x = np.arange(n, dtype=float)
    x_norm = x / (x.max() + 1e-9)
    y_norm = (y - y.min()) / (y.max() - y.min() + 1e-9)
    distances = 1.0 - (x_norm + y_norm)
    elbow_idx = int(np.argmax(distances))  # 0-based

    return max(elbow_idx + 1, min_kept)


def curate_markers(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Per dataset x cluster: sort by stringent_score desc, apply elbow
    detection, keep the resulting top genes.
    """
    mf = cfg["marker_filtering"]
    curated_parts = []

    for (dataset, cluster), grp in df.groupby(["dataset", "cluster"], sort=False):
        grp_sorted = grp.sort_values("stringent_score", ascending=False).reset_index(drop=True)
        n_keep = detect_elbow(
            grp_sorted["stringent_score"].values,
            max_show=mf["max_n_shown"],
            min_kept=mf["min_kept"],
        )
        kept = grp_sorted.head(n_keep).copy()
        kept["rank"] = np.arange(1, len(kept) + 1)
        curated_parts.append(kept)

    if not curated_parts:
        return pd.DataFrame()

    curated = pd.concat(curated_parts, ignore_index=True)
    curated["stringent_score"] = curated["stringent_score"].round(2)
    curated["pct_in"] = curated["pct_in"].round(3)
    curated["pct_out"] = curated["pct_out"].round(3)
    curated["logFC"] = curated["logFC"].round(2)
    return curated


# =============================================================================
# 6 . Export
# =============================================================================

def export_per_dataset_tsv(curated: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for dataset, grp in curated.groupby("dataset"):
        out_path = out_dir / f"{dataset}_markers.tsv"
        grp.sort_values(["cluster", "rank"]).to_csv(out_path, sep="\t", index=False)
        log.info("    wrote %s (%d rows)", out_path.name, len(grp))


# =============================================================================
# Orchestration
# =============================================================================

def process_dataset(path: Path, cfg: dict, noise_filter: NoiseGeneFilter) -> pd.DataFrame | None:
    dataset_name = dataset_name_from_path(path)
    dynamic_enabled = cfg["clustering"].get("dynamic_resolution", {}).get("enabled", False)

    resolution = None
    if not dynamic_enabled:
        resolutions = cfg["clustering"]["resolutions"]
        if dataset_name not in resolutions:
            log.warning("[%s] no resolution configured — skipping", dataset_name)
            return None
        resolution = resolutions[dataset_name]

    log.info("=" * 70)
    if dynamic_enabled:
        log.info("Dataset: %s (dynamic resolution search)", dataset_name)
    else:
        log.info("Dataset: %s (resolution=%.2f)", dataset_name, resolution)
    log.info("=" * 70)

    adata = load_adata(path, counts_layer=cfg["marker_extraction"]["counts_layer"])
    log.info("  loaded: %d cells x %d genes", adata.n_obs, adata.n_vars)

    adata = preprocess_adata(adata, cfg)

    if dynamic_enabled:
        adata = cluster_adata_dynamic(adata, cfg, dataset_name)
    else:
        adata = cluster_adata(adata, resolution, cfg)

    raw_markers = extract_wilcoxon_markers(
        adata,
        dataset_name,
        cfg["clustering"]["cluster_key"],
        noise_filter,
        cfg,
    )
    if raw_markers.empty:
        log.warning("  [%s] no markers extracted", dataset_name)
        return None

    filtered = apply_baseline_filters(raw_markers, cfg)
    scored = compute_stringent_score(filtered, cfg["marker_filtering"]["tolerance"])
    curated = curate_markers(scored, cfg)

    log.info(
        "  [%s] %d clusters, %d curated marker rows",
        dataset_name,
        curated["cluster"].nunique() if not curated.empty else 0,
        len(curated),
    )
    return curated


def run_phase1(config_path: str | Path) -> None:
    """Phase 1: per-dataset preprocessing, clustering, Wilcoxon marker extraction."""
    cfg = load_config(config_path)
    io_cfg = cfg["io"]

    input_dir = Path(io_cfg["input_dir"])
    output_dir = Path(io_cfg["output_dir"])
    per_dataset_dir = output_dir / io_cfg["markers_per_dataset_dir"]

    noise_filter = NoiseGeneFilter.from_config(cfg)
    h5ad_files = discover_h5ad_files(input_dir)

    if not h5ad_files:
        log.error("No .h5ad files found in %s — nothing to do.", input_dir)
        return

    dynamic_enabled = cfg["clustering"].get("dynamic_resolution", {}).get("enabled", False)
    if dynamic_enabled:
        log.warning(
            "Dynamic resolution clustering is ON — this searches a resolution "
            "per dataset instead of using a fixed one, so it will take noticeably "
            "longer than a fixed-resolution run."
        )
    else:
        configured = set(cfg["clustering"]["resolutions"].keys())
        discovered = {dataset_name_from_path(p) for p in h5ad_files}
        unmatched_config = configured - discovered
        if unmatched_config:
            log.warning(
                "%d resolution entries in config.yaml don't match any discovered "
                "dataset: %s",
                len(unmatched_config),
                sorted(unmatched_config),
            )

    all_curated = []
    for path in h5ad_files:
        try:
            curated = process_dataset(path, cfg, noise_filter)
        except Exception as exc:  # one bad dataset shouldn't kill the whole run
            log.exception("[%s] FAILED: %s", dataset_name_from_path(path), exc)
            continue
        if curated is not None and not curated.empty:
            all_curated.append(curated)

    if not all_curated:
        log.error("No curated markers produced for any dataset.")
        return

    combined = pd.concat(all_curated, ignore_index=True)

    export_per_dataset_tsv(combined, per_dataset_dir)

    log.info("=" * 70)
    log.info("PHASE 1 COMPLETE")
    log.info("  Datasets processed : %d", combined["dataset"].nunique())
    log.info("  Total marker rows  : %d", len(combined))
    log.info("  Per-dataset TSVs   -> %s", per_dataset_dir)
    log.info("=" * 70)


# =============================================================================
# PHASE 2 — Organoid <-> Adult/Fetal Reference Mapping (Pearson correlation)
# =============================================================================
# Each organoid cluster's curated marker list (Phase 1 output) is correlated
# against adult and fetal reference cell types. We use Pearson correlation
# over the genes shared between an organoid cluster and a reference cell
# type (organoid stringent_score vs reference score), rather than a plain
# Jaccard overlap — this rewards references whose relative marker strength
# tracks the organoid's, not just raw gene overlap.
#
# Each reference (adult / fetal) is a single TSV with one row per
# (cell_type, gene). Column names are configurable via phase2.reference_columns
# in config.yaml (defaults: canonical_cell_type, gene, score).
# =============================================================================


def load_reference_markers(path: str | Path, columns: dict) -> dict[str, pd.Series]:
    """
    Load a single reference marker TSV (adult or fetal) with one row per
    (cell_type, gene, score). Returns dict: cell_type -> pd.Series(gene -> score).
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t")

    ct_col, gene_col, score_col = columns["cell_type"], columns["gene"], columns["score"]
    missing = {ct_col, gene_col, score_col} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing expected column(s) {missing}")

    reference = {}
    for cell_type, sub in df.groupby(ct_col):
        sub = sub.dropna(subset=[gene_col, score_col]).drop_duplicates(subset=[gene_col])
        reference[cell_type] = pd.Series(
            sub[score_col].values, index=sub[gene_col].values
        )
    return reference


def load_organoid_marker_tsvs(folder: str | Path) -> dict[str, dict[str, pd.DataFrame]]:
    """
    Load phase-1 output (<dataset>_markers.tsv). Returns nested dict:
    dataset -> cluster -> DataFrame(gene, stringent_score, rank), sorted by rank.
    """
    folder = Path(folder)
    tsv_files = sorted(folder.glob("*_markers.tsv"))
    if not tsv_files:
        raise FileNotFoundError(f"No *_markers.tsv files found in {folder}")

    result: dict[str, dict[str, pd.DataFrame]] = {}
    for tsv_path in tsv_files:
        dataset = tsv_path.stem.replace("_markers", "")
        df = pd.read_csv(tsv_path, sep="\t")

        required = {"gene", "cluster", "stringent_score", "rank"}
        missing = required - set(df.columns)
        if missing:
            log.warning("    %s: missing columns %s — skipping", tsv_path.name, missing)
            continue

        cluster_dict = {}
        for cluster, sub in df.groupby("cluster"):
            sub = sub.sort_values("rank")
            cluster_dict[str(cluster)] = sub[["gene", "stringent_score"]].reset_index(drop=True)
        result[dataset] = cluster_dict

    return result


def pearson_match(
    cluster_scores: pd.Series,
    reference: dict[str, pd.Series],
    min_overlap_genes: int,
) -> list[tuple[str, float, int, list[str]]]:
    """
    Correlate the organoid cluster's gene scores against every reference cell
    type's gene scores, restricted to the genes they share.

    Returns a list of (cell_type, pearson_r, n_overlap, overlap_genes_in_organoid_rank_order),
    sorted best-first (highest r, then largest overlap). Cell types with fewer
    than `min_overlap_genes` shared genes are excluded (Pearson r is not
    meaningful / not computable on too few points).
    """
    matches = []
    for cell_type, ref_scores in reference.items():
        shared = cluster_scores.index.intersection(ref_scores.index)
        n_shared = len(shared)
        if n_shared < min_overlap_genes:
            continue

        r = cluster_scores.loc[shared].corr(ref_scores.loc[shared], method="pearson")
        if pd.isna(r):
            continue

        overlap_ordered = [g for g in cluster_scores.index if g in shared]
        matches.append((cell_type, float(r), n_shared, overlap_ordered))

    matches.sort(key=lambda x: (x[1], x[2]), reverse=True)
    return matches


def build_mapping_table(cfg: dict) -> pd.DataFrame:
    p2 = cfg["phase2"]
    io_cfg = cfg["io"]

    ref_columns = p2["reference_columns"]
    adult_ref = load_reference_markers(p2["adult_markers_path"], ref_columns)
    fetal_ref = load_reference_markers(p2["fetal_markers_path"], ref_columns)
    log.info(
        "Loaded references: %d adult cell types, %d fetal cell types",
        len(adult_ref), len(fetal_ref),
    )

    organoid_dir = Path(io_cfg["output_dir"]) / io_cfg["markers_per_dataset_dir"]
    organoid_markers = load_organoid_marker_tsvs(organoid_dir)
    total_clusters = sum(len(c) for c in organoid_markers.values())
    log.info(
        "Loaded organoid markers: %d datasets, %d clusters",
        len(organoid_markers), total_clusters,
    )

    top_n = p2["top_n_genes"]
    min_overlap = p2["min_overlap_genes"]

    rows = []
    for dataset, clusters in organoid_markers.items():
        for cluster, cluster_df in clusters.items():
            top_df = cluster_df.head(top_n)
            top_genes_ordered = top_df["gene"].tolist()
            cluster_scores = pd.Series(
                top_df["stringent_score"].values, index=top_df["gene"].values
            )

            adult_matches = pearson_match(cluster_scores, adult_ref, min_overlap)
            fetal_matches = pearson_match(cluster_scores, fetal_ref, min_overlap)

            if adult_matches:
                a_ct, a_r, _, a_overlap = adult_matches[0]
            else:
                a_ct, a_r, a_overlap = None, np.nan, []

            if fetal_matches:
                f_ct, f_r, _, f_overlap = fetal_matches[0]
            else:
                f_ct, f_r, f_overlap = None, np.nan, []

            trinity = [g for g in a_overlap if g in set(f_overlap)]

            rows.append(
                {
                    "Organoid_dataset": dataset,
                    "Cluster": cluster,
                    "Top_Marker_Genes": ", ".join(top_genes_ordered),
                    "Adult_Match": a_ct,
                    "Adult_Score": round(a_r, 4) if pd.notna(a_r) else np.nan,
                    "Adult_Overlap_Genes": ", ".join(a_overlap),
                    "Fetal_Match": f_ct,
                    "Fetal_Score": round(f_r, 4) if pd.notna(f_r) else np.nan,
                    "Fetal_Overlap_Genes": ", ".join(f_overlap),
                    "Trinity_Overlap_Genes": ", ".join(trinity),
                }
            )

    return pd.DataFrame(
        rows,
        columns=[
            "Organoid_dataset",
            "Cluster",
            "Top_Marker_Genes",
            "Adult_Match",
            "Adult_Score",
            "Adult_Overlap_Genes",
            "Fetal_Match",
            "Fetal_Score",
            "Fetal_Overlap_Genes",
            "Trinity_Overlap_Genes",
        ],
    )


def run_phase2(config_path: str | Path) -> None:
    """Phase 2: correlate organoid clusters against adult/fetal references."""
    cfg = load_config(config_path)
    io_cfg = cfg["io"]

    mapping_table = build_mapping_table(cfg)
    if mapping_table.empty:
        log.error("Phase 2 produced no rows — check phase-1 outputs and reference paths.")
        return

    out_dir = Path(io_cfg["output_dir"]) / cfg["phase2"]["mapping_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "mapping_table.tsv"
    mapping_table.to_csv(out_path, sep="\t", index=False)

    log.info("=" * 70)
    log.info("PHASE 2 COMPLETE")
    log.info("  Clusters mapped : %d", len(mapping_table))
    log.info("  Adult matches   : %d", mapping_table["Adult_Match"].notna().sum())
    log.info("  Fetal matches   : %d", mapping_table["Fetal_Match"].notna().sum())
    log.info("  Mapping table   -> %s", out_path)
    log.info("=" * 70)


# =============================================================================
# Orchestration — run phase 1 and/or phase 2
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(description="Organoid marker extraction & reference mapping pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument(
        "--phase",
        choices=["1", "2", "both"],
        default="both",
        help="Which phase to run: 1 (marker extraction), 2 (reference mapping), or both.",
    )
    args = parser.parse_args()

    if args.phase in ("1", "both"):
        run_phase1(args.config)
    if args.phase in ("2", "both"):
        run_phase2(args.config)


if __name__ == "__main__":
    main()
