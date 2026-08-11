"""Snapseed hierarchical annotation for QC-filtered or integrated h5ad objects."""

from __future__ import annotations

import gc
import logging
import random
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from hkoca.annotation.config import load_marker_hierarchy

logger = logging.getLogger("hkoca.annotation")


def resolution_tag(resolution: float) -> str:
    """Stable tag for filenames / obs keys, e.g. 0.4 -> '0.4', 1.0 -> '1.0'."""
    return f"{float(resolution):.1f}"


def cluster_key_for(resolution: float) -> str:
    return f"leiden_res_{resolution_tag(resolution)}"


def level_key_for(level: int, resolution: float) -> str:
    return f"Level_{level}_res{resolution_tag(resolution)}"


def latest_key_for(resolution: float) -> str:
    return f"Level_latest_res{resolution_tag(resolution)}"


def discover_h5ad_inputs(input_path: Path | str) -> list[Path]:
    path = Path(input_path).expanduser().resolve()
    if path.is_file():
        if path.suffix.lower() != ".h5ad":
            raise ValueError(f"Expected an .h5ad file, got: {path}")
        return [path]
    if path.is_dir():
        files = sorted(path.glob("*.h5ad"))
        if not files:
            raise FileNotFoundError(f"No .h5ad files found in: {path}")
        return files
    # Treat as glob relative to CWD if parent exists.
    matches = sorted(Path().glob(str(input_path)))
    matches = [m.resolve() for m in matches if m.suffix.lower() == ".h5ad" and m.is_file()]
    if not matches:
        raise FileNotFoundError(f"No .h5ad inputs matched: {input_path}")
    return matches


def _matrix_extremes(x) -> tuple[float, float]:
    if hasattr(x, "data") and getattr(x, "nnz", 1) >= 0 and not isinstance(x, np.ndarray):
        data = x.data
        if data.size == 0:
            return 0.0, 0.0
        return float(np.min(data)), float(np.max(data))
    arr = np.asarray(x)
    if arr.size == 0:
        return 0.0, 0.0
    return float(np.min(arr)), float(np.max(arr))


def _ensure_gene_names(adata) -> None:
    if "features" in adata.var.columns:
        adata.var_names = adata.var["features"].astype(str)
    adata.var_names_make_unique()


def _recover_unscaled_expression(adata):
    """If .X looks scaled (negatives), restore from .raw or a counts-like layer."""
    min_val, _ = _matrix_extremes(adata.X)
    if min_val >= 0:
        return adata

    logger.warning("Negative values in adata.X; attempting to restore unscaled expression.")
    if adata.raw is not None:
        adata = adata.raw.to_adata()
        adata.obsm.clear()
        adata.varm.clear()
        return adata

    for layer in ("counts", "raw", "normalized"):
        if layer in adata.layers:
            logger.info("Resetting adata.X from layer '%s'.", layer)
            adata.X = adata.layers[layer].copy()
            return adata

    raise ValueError(
        "adata.X appears scaled (negative values) and no adata.raw / counts layer is available."
    )


def _normalize_if_needed(adata, *, target_sum: float) -> None:
    import scanpy as sc

    _, max_val = _matrix_extremes(adata.X)
    if max_val > 100:
        logger.info("Raw-count footprint detected (max=%.1f); normalize_total + log1p.", max_val)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    else:
        logger.info("Matrix looks pre-normalized (max=%.1f); skipping log-norm.", max_val)


def _prepare_embedding(
    adata,
    *,
    seed: int,
    hvg_n_top_genes: int,
    pca_n_comps: int,
    neighbors_n_neighbors: int,
    neighbors_n_pcs: int,
    scale_max_value: float,
):
    import scanpy as sc

    logger.info("Selecting HVGs (n_top_genes=%d).", hvg_n_top_genes)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg_n_top_genes, subset=False)
    n_hvg = int(adata.var["highly_variable"].sum())
    logger.info("HVGs retained: %d", n_hvg)

    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=scale_max_value)
    x = adata_hvg.X
    if isinstance(x, np.ndarray):
        if np.isnan(x).any():
            adata_hvg.X = np.nan_to_num(x)
    elif hasattr(x, "data") and np.isnan(x.data).any():
        x.data = np.nan_to_num(x.data)

    sc.pp.pca(adata_hvg, n_comps=pca_n_comps, svd_solver="arpack", random_state=seed)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()
    del adata_hvg
    gc.collect()

    sc.pp.neighbors(adata, n_neighbors=neighbors_n_neighbors, n_pcs=neighbors_n_pcs)
    sc.tl.umap(adata, random_state=seed)


def _resolve_assignment_col(assignments: pd.DataFrame, level_key: str) -> str | None:
    if level_key in assignments.columns:
        return level_key
    lookup = {str(c).strip().lower(): c for c in assignments.columns}
    key = str(level_key).strip().lower()
    if key in lookup:
        return lookup[key]
    if key.startswith("level_"):
        alt = "level" + key[len("level_") :]
        if alt in lookup:
            return lookup[alt]
    if key.startswith("level") and not key.startswith("level_"):
        alt = key.replace("level", "level_", 1)
        if alt in lookup:
            return lookup[alt]
    return None


def _apply_manual_overrides(
    assignments: pd.DataFrame,
    manual_annotations: dict[str, dict[str, Any]] | None,
) -> pd.DataFrame:
    if not manual_annotations:
        return assignments
    if not isinstance(manual_annotations, dict):
        logger.warning("manual_annotations is not a dict; skipping overrides.")
        return assignments

    for cluster_id, levels in manual_annotations.items():
        cid = str(cluster_id)
        if cid not in set(map(str, assignments.index)):
            logger.warning("Cluster %s not in Snapseed assignments; skip override.", cid)
            continue
        for level_name, value in dict(levels).items():
            col = _resolve_assignment_col(assignments, str(level_name))
            if col is None:
                logger.warning("Level '%s' not in assignments columns; skip.", level_name)
                continue
            assignments.loc[cid, col] = value
    return assignments


def _map_levels_to_obs(
    adata,
    *,
    cluster_key: str,
    assignments: pd.DataFrame,
    resolution: float,
) -> list[str]:
    annotation_cols: list[str] = []
    for i, col in enumerate(assignments.columns):
        obs_col = level_key_for(i + 1, resolution)
        mapped = adata.obs[cluster_key].astype(str).map(assignments[col].to_dict())
        adata.obs[obs_col] = pd.Categorical(mapped.astype(object))
        annotation_cols.append(obs_col)

    latest = assignments.apply(
        lambda row: row.dropna().iloc[-1] if row.dropna().shape[0] > 0 else np.nan,
        axis=1,
    )
    fill = adata.obs[cluster_key].astype(str).map(latest.to_dict()).astype(object)
    for fallback in annotation_cols[::-1]:
        fill = fill.fillna(adata.obs[fallback].astype(object))
    fill = fill.fillna(adata.obs[cluster_key].astype(str))
    latest_col = latest_key_for(resolution)
    adata.obs[latest_col] = pd.Categorical(fill.astype(object))
    annotation_cols.append(latest_col)
    return annotation_cols


def _save_umap_png(adata, color_key: str, out_path: Path, *, title: str, dpi: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if "X_umap" not in adata.obsm:
        return
    umap = np.asarray(adata.obsm["X_umap"])
    labels = adata.obs[color_key].astype(str)
    cats = list(pd.Categorical(labels).categories)
    cmap = plt.get_cmap("tab20", max(len(cats), 1))
    colors = {c: cmap(i) for i, c in enumerate(cats)}

    fig, ax = plt.subplots(figsize=(7, 5.5))
    for c in cats:
        mask = labels == c
        if not mask.any():
            continue
        ax.scatter(umap[mask, 0], umap[mask, 1], s=4, c=[colors[c]], rasterized=True, alpha=0.85)
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def annotate_dataset(
    file_path: Path | str,
    marker_dict: dict | None = None,
    *,
    markers_path: Path | str | None = None,
    resolutions: Sequence[float] = (0.4, 0.6, 1.0),
    output_dir: Path | str,
    clustered_dir: Path | str | None = None,
    figures_dir: Path | str | None = None,
    manual_annotations: dict[str, dict[str, Any]] | None = None,
    seed: int = 42,
    hvg_n_top_genes: int = 3000,
    pca_n_comps: int = 50,
    neighbors_n_neighbors: int = 20,
    neighbors_n_pcs: int = 30,
    scale_max_value: float = 10,
    normalize_target_sum: float = 10000,
    save_plots: bool = False,
    dpi: int = 150,
    skip_existing: bool = True,
    force: bool = False,
) -> Path:
    """Cluster at each resolution, run Snapseed, write annotated (+ optional clustered) h5ad."""
    import anndata as ad
    import scanpy as sc

    try:
        import snapseed
    except ImportError as exc:
        raise ImportError(
            "Snapseed is required for annotation. Install with: "
            "pip install -U git+https://github.com/hassansaei/snapseed.git"
        ) from exc

    try:
        import igraph  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "python-igraph is required for Leiden clustering. Install with: "
            "conda install -c conda-forge python-igraph leidenalg "
            "or: pip install igraph leidenalg"
        ) from exc

    path = Path(file_path).expanduser().resolve()
    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    clustered_dir_p = (
        Path(clustered_dir).expanduser().resolve()
        if clustered_dir
        else out_dir.parent / "clustered"
    )
    figures_dir_p = (
        Path(figures_dir).expanduser().resolve()
        if figures_dir
        else out_dir.parent / "figures"
    )

    out_path = out_dir / f"{path.stem}_annotated.h5ad"
    if skip_existing and not force and out_path.is_file() and out_path.stat().st_size > 0:
        logger.info("[SKIP] Annotated output exists: %s", out_path)
        return out_path

    if marker_dict is None:
        marker_dict = load_marker_hierarchy(markers_path)

    random.seed(seed)
    np.random.seed(seed)
    if hasattr(sc.settings, "seed"):
        sc.settings.seed = seed
    try:
        ad.settings.allow_write_nullable_strings = True
    except Exception:
        pass

    logger.info("Processing %s", path.name)
    adata = sc.read_h5ad(path)
    _ensure_gene_names(adata)
    logger.info("Loaded %s cells x %s genes", f"{adata.n_obs:,}", f"{adata.n_vars:,}")

    adata = _recover_unscaled_expression(adata)
    _normalize_if_needed(adata, target_sum=normalize_target_sum)
    _prepare_embedding(
        adata,
        seed=seed,
        hvg_n_top_genes=hvg_n_top_genes,
        pca_n_comps=pca_n_comps,
        neighbors_n_neighbors=neighbors_n_neighbors,
        neighbors_n_pcs=neighbors_n_pcs,
        scale_max_value=scale_max_value,
    )

    if "_index" in adata.var.columns:
        adata.var.rename(columns={"_index": "var_index"}, inplace=True)
    if adata.raw is not None and "_index" in adata.raw.var.columns:
        adata.raw.var.rename(columns={"_index": "raw_index"}, inplace=True)

    res_list = [float(r) for r in resolutions]
    adata.uns["snapseed_resolutions"] = res_list
    adata.uns["snapseed_markers_path"] = str(markers_path or "")

    assignment_tables: dict[str, pd.DataFrame] = {}

    for resolution in res_list:
        ckey = cluster_key_for(resolution)
        logger.info("Leiden clustering resolution=%s -> %s", resolution_tag(resolution), ckey)
        sc.tl.leiden(adata, resolution=resolution, key_added=ckey, random_state=seed)
        adata.obs[ckey] = adata.obs[ckey].astype(str).astype("category")
        n_clusters = int(adata.obs[ckey].nunique())
        logger.info("  %d clusters", n_clusters)

        if save_plots:
            _save_umap_png(
                adata,
                ckey,
                figures_dir_p / path.stem / f"umap_{ckey}.png",
                title=f"{path.stem} | {ckey}",
                dpi=dpi,
            )

        logger.info("  Wilcoxon rank_genes_groups on %s", ckey)
        sc.tl.rank_genes_groups(adata, groupby=ckey, use_raw=False, method="wilcoxon", pts=True)

        logger.info("  Snapseed annotate_hierarchy on %s", ckey)
        results = snapseed.annotate_hierarchy(adata, marker_dict, group_name=ckey)
        assignments = results["assignments"].copy()
        metrics_df = results["metrics"]
        assignments.index = assignments.index.astype(str)
        assignments = _apply_manual_overrides(assignments, manual_annotations)
        assignments = assignments.fillna("Unassigned").astype(str)

        tag = resolution_tag(resolution)
        adata.uns[f"snapseed_assignments_res{tag}"] = assignments
        adata.uns[f"snapseed_metrics_res{tag}"] = metrics_df
        assignment_tables[tag] = assignments

        cols = _map_levels_to_obs(
            adata, cluster_key=ckey, assignments=assignments, resolution=resolution
        )
        logger.info("  Wrote obs columns: %s", ", ".join(cols))

        if save_plots:
            latest = latest_key_for(resolution)
            _save_umap_png(
                adata,
                latest,
                figures_dir_p / path.stem / f"umap_{latest}.png",
                title=f"{path.stem} | {latest}",
                dpi=dpi,
            )

    # Convenience aliases from the middle / last resolution for downstream tools.
    primary = res_list[1] if len(res_list) >= 2 else res_list[-1]
    for src, dst in (
        (cluster_key_for(primary), "leiden_clusters"),
        (latest_key_for(primary), "Level_3_latest"),
        (level_key_for(1, primary), "Level_1"),
        (level_key_for(2, primary), "Level_2"),
        (level_key_for(3, primary), "Level_3"),
    ):
        if src in adata.obs.columns:
            adata.obs[dst] = adata.obs[src].astype("category")

    clustered_dir_p.mkdir(parents=True, exist_ok=True)
    clustered_path = clustered_dir_p / f"{path.stem}_clustered.h5ad"
    adata.write_h5ad(clustered_path)
    logger.info("Clustered object: %s", clustered_path)

    adata.write_h5ad(out_path)
    logger.info("Annotated object: %s", out_path)

    # Compact assignment CSV for all resolutions.
    csv_rows = []
    for tag, assignments in assignment_tables.items():
        tmp = assignments.copy()
        tmp.index.name = "cluster"
        tmp = tmp.reset_index()
        tmp.insert(0, "resolution", tag)
        csv_rows.append(tmp)
    if csv_rows:
        csv_path = out_dir / f"{path.stem}_snapseed_assignments.csv"
        pd.concat(csv_rows, ignore_index=True).to_csv(csv_path, index=False)
        logger.info("Assignments CSV: %s", csv_path)

    return out_path


def run_annotation_batch(
    inputs: Iterable[Path | str],
    *,
    markers_path: Path | str | None = None,
    resolutions: Sequence[float] = (0.4, 0.6, 1.0),
    annotated_dir: Path | str,
    clustered_dir: Path | str | None = None,
    figures_dir: Path | str | None = None,
    parameters: dict[str, Any] | None = None,
    force: bool = False,
) -> list[Path]:
    params = dict(parameters or {})
    marker_dict = load_marker_hierarchy(markers_path)
    outputs: list[Path] = []
    for inp in inputs:
        out = annotate_dataset(
            inp,
            marker_dict=marker_dict,
            markers_path=markers_path,
            resolutions=params.get("resolutions", resolutions),
            output_dir=annotated_dir,
            clustered_dir=clustered_dir,
            figures_dir=figures_dir,
            seed=int(params.get("seed", 42)),
            hvg_n_top_genes=int(params.get("hvg_n_top_genes", 3000)),
            pca_n_comps=int(params.get("pca_n_comps", 50)),
            neighbors_n_neighbors=int(params.get("neighbors_n_neighbors", 20)),
            neighbors_n_pcs=int(params.get("neighbors_n_pcs", 30)),
            scale_max_value=float(params.get("scale_max_value", 10)),
            normalize_target_sum=float(params.get("normalize_target_sum", 10000)),
            save_plots=bool(params.get("save_plots", False)),
            dpi=int(params.get("dpi", 150)),
            skip_existing=bool(params.get("skip_existing", True)),
            force=force,
        )
        outputs.append(out)
    return outputs
