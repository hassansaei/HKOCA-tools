"""Project query datasets onto a reference atlas (scanpy ingest)."""

from __future__ import annotations

import json
import logging
import random
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from hkoca.projection.config import DEFAULT_LABEL_PRIORITY, load_projection_config

logger = logging.getLogger("hkoca.projection")


def _matrix_extremes(x) -> tuple[float, float]:
    if hasattr(x, "min"):
        return float(x.min()), float(x.max())
    arr = np.asarray(x)
    return float(arr.min()), float(arr.max())


def _ensure_gene_names(adata) -> None:
    if "features" in adata.var.columns:
        adata.var_names = adata.var["features"].astype(str)
    elif "gene_symbols" in adata.var.columns:
        adata.var_names = adata.var["gene_symbols"].astype(str)
    adata.var_names_make_unique()


def _normalize_if_needed(adata, *, target_sum: float) -> None:
    import scanpy as sc

    _, max_val = _matrix_extremes(adata.X)
    if max_val > 100:
        logger.info("Raw counts detected (max=%.1f); normalize_total + log1p.", max_val)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    else:
        logger.info("Matrix looks normalized (max=%.1f); skipping normalize_total.", max_val)


def _prepare_reference(adata, cfg: dict[str, Any]) -> None:
    import scanpy as sc

    _normalize_if_needed(adata, target_sum=cfg["normalize_target_sum"])
    if "X_pca" not in adata.obsm or "neighbors" not in adata.uns:
        logger.info("Reference missing PCA/neighbors; computing embedding on atlas.")
        sc.pp.highly_variable_genes(adata, n_top_genes=cfg["hvg_n_top_genes"], subset=True)
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=cfg["pca_n_comps"])
        sc.pp.neighbors(
            adata,
            n_neighbors=cfg["neighbors_n_neighbors"],
            n_pcs=cfg["neighbors_n_pcs"],
        )
    if "X_umap" not in adata.obsm:
        logger.info("Reference missing UMAP; computing UMAP on atlas.")
        sc.tl.umap(adata)


def _align_genes(reference, query):
    common = reference.var_names.intersection(query.var_names)
    if len(common) < 500:
        raise ValueError(
            f"Too few shared genes between atlas and query ({len(common)}). "
            "Check gene naming / harmonization."
        )
    logger.info("Shared genes for projection: %d", len(common))
    ref = reference[:, common].copy()
    qry = query[:, common].copy()
    qry = qry[:, ref.var_names].copy()
    return ref, qry


def _detect_label_columns(reference, requested: Sequence[str] | None) -> list[str]:
    if requested:
        cols = [c for c in requested if c in reference.obs.columns]
        missing = [c for c in requested if c not in reference.obs.columns]
        if missing:
            logger.warning("Reference label columns not found (skipped): %s", ", ".join(missing))
        if not cols:
            raise ValueError("None of the requested reference label columns exist in the atlas.")
        return cols
    for col in DEFAULT_LABEL_PRIORITY:
        if col in reference.obs.columns:
            logger.info("Auto-selected reference label column: %s", col)
            return [col]
    raise ValueError(
        "Could not detect a cell-type column on the atlas. "
        f"Expected one of: {', '.join(DEFAULT_LABEL_PRIORITY)}"
    )


def _setup_logging_file(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if any(isinstance(h, logging.FileHandler) and h.baseFilename == str(log_path) for h in logger.handlers):
        return
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"))
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)


def _save_umap_scatter(
    umap: np.ndarray,
    labels: pd.Series,
    out_path: Path,
    *,
    title: str,
    dpi: int,
    palette: dict[str, str] | None = None,
    alpha: float = 0.55,
    figsize: tuple[float, float] = (8, 6),
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    cats = [c for c in pd.Categorical(labels.astype(str)).categories if (labels.astype(str) == c).any()]
    if palette is None:
        cmap = plt.get_cmap("tab20", max(len(cats), 1))
        colors = {c: cmap(i) for i, c in enumerate(cats)}
    else:
        from hkoca.annotation.colors import fallback_color, lookup_color

        colors = {c: lookup_color(c, palette) for c in cats}

    fig, ax = plt.subplots(figsize=figsize)
    lbl = labels.astype(str)
    for cat in cats:
        mask = lbl == cat
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            s=4,
            c=[colors[cat]],
            rasterized=True,
            alpha=alpha,
        )
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    if palette is not None and cats:
        handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[c], markersize=6, label=c)
            for c in cats
        ]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)
        fig.tight_layout(rect=(0, 0, 0.72, 1))
    else:
        fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_overlay(
    ref_umap: np.ndarray,
    query_umap: np.ndarray,
    query_labels: pd.Series,
    out_path: Path,
    *,
    title: str,
    dpi: int,
    palette: dict[str, str] | None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(ref_umap[:, 0], ref_umap[:, 1], s=2, c="#DDDDDD", rasterized=True, alpha=0.35, label="Atlas")
    lbl = query_labels.astype(str)
    cats = [c for c in pd.Categorical(lbl).categories if (lbl == c).any()]
    if palette is None:
        cmap = plt.get_cmap("tab20", max(len(cats), 1))
        colors = {c: cmap(i) for i, c in enumerate(cats)}
    else:
        from hkoca.annotation.colors import lookup_color

        colors = {c: lookup_color(c, palette) for c in cats}
    for cat in cats:
        mask = lbl == cat
        ax.scatter(
            query_umap[mask, 0],
            query_umap[mask, 1],
            s=6,
            c=[colors[cat]],
            rasterized=True,
            alpha=0.75,
            label=cat,
        )
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0, 0.75, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _composition_table(labels: pd.Series, name: str) -> pd.DataFrame:
    counts = labels.astype(str).value_counts()
    total = counts.sum()
    return pd.DataFrame(
        {
            "dataset": name,
            "label": counts.index.astype(str),
            "n_cells": counts.values,
            "fraction": (counts.values / total).round(4),
        }
    )


def _plot_composition(ref_comp: pd.DataFrame, query_comp: pd.DataFrame, out_path: Path, *, dpi: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    merged = pd.concat([ref_comp, query_comp], ignore_index=True)
    labels = sorted(merged["label"].unique(), key=str)
    datasets = merged["dataset"].unique().tolist()
    x = np.arange(len(labels))
    width = 0.35 if len(datasets) == 2 else 0.8 / max(len(datasets), 1)

    fig, ax = plt.subplots(figsize=(max(10, len(labels) * 0.45), 5))
    for i, ds in enumerate(datasets):
        sub = merged[merged["dataset"] == ds].set_index("label").reindex(labels).fillna(0)
        offset = (i - (len(datasets) - 1) / 2) * width
        ax.bar(x + offset, sub["fraction"].values, width=width, label=ds)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Fraction of cells")
    ax.set_title("Cell-type composition: atlas vs projected query")
    ax.legend()
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _plot_confusion(
    truth: pd.Series,
    predicted: pd.Series,
    out_path: Path,
    *,
    dpi: int,
) -> pd.DataFrame:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.crosstab(truth.astype(str), predicted.astype(str))
    fig, ax = plt.subplots(figsize=(max(6, df.shape[1] * 0.5), max(5, df.shape[0] * 0.5)))
    im = ax.imshow(df.values, aspect="auto", cmap="Blues")
    ax.set_xticks(range(df.shape[1]))
    ax.set_yticks(range(df.shape[0]))
    ax.set_xticklabels(df.columns, rotation=45, ha="right")
    ax.set_yticklabels(df.index)
    ax.set_xlabel("Projected label")
    ax.set_ylabel("Query label (ground truth)")
    ax.set_title("Projection confusion matrix")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return df


def _palette_for_label(label_col: str) -> dict[str, str] | None:
    from hkoca.annotation.colors import (
        combined_color_map,
        level1_color_map,
        level2_color_map,
        level3_color_map,
        palette_for_obs_key,
    )

    mapped = palette_for_obs_key(label_col)
    if mapped is not None:
        return mapped
    if label_col in {"celltype_final", "Level_3", "Level_3_latest"}:
        return level3_color_map()
    if label_col == "Level_2":
        return level2_color_map()
    if label_col == "Level_1":
        return level1_color_map()
    return combined_color_map()


def project_query(
    *,
    atlas_h5ad: str | Path,
    query_h5ad: str | Path,
    output_dir: str | Path,
    label_columns: Sequence[str] | None = None,
    query_label_column: str | None = None,
    query_batch_key: str | None = "sample_id",
    cfg: dict[str, Any] | None = None,
    force: bool = False,
) -> Path:
    import scanpy as sc

    atlas_path = Path(atlas_h5ad).expanduser().resolve()
    query_path = Path(query_h5ad).expanduser().resolve()
    out_root = Path(output_dir).expanduser().resolve()
    cfg = cfg or load_projection_config()

    projected_dir = out_root / cfg["projected_subdir"]
    figures_dir = out_root / cfg["figures_subdir"]
    tables_dir = out_root / cfg["tables_subdir"]
    logs_dir = out_root / cfg["logs_subdir"]
    for d in (projected_dir, figures_dir, tables_dir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)

    out_h5ad = projected_dir / f"{query_path.stem}_projected.h5ad"
    meta_json = projected_dir / f".{query_path.stem}_projected.meta.json"
    _setup_logging_file(logs_dir / "projection.log")

    if not force and cfg.get("skip_existing") and out_h5ad.is_file() and out_h5ad.stat().st_size > 0:
        logger.info("Skipping existing projected output: %s", out_h5ad)
        return out_h5ad

    logger.info("Loading atlas: %s", atlas_path)
    reference = sc.read_h5ad(atlas_path)
    _ensure_gene_names(reference)
    logger.info("Atlas: %s cells x %s genes", f"{reference.n_obs:,}", f"{reference.n_vars:,}")

    logger.info("Loading query: %s", query_path)
    query = sc.read_h5ad(query_path)
    _ensure_gene_names(query)
    logger.info("Query: %s cells x %s genes", f"{query.n_obs:,}", f"{query.n_vars:,}")

    random.seed(cfg["seed"])
    np.random.seed(cfg["seed"])
    sc.settings.seed = cfg["seed"]

    reference, query = _align_genes(reference, query)
    _prepare_reference(reference, cfg)
    _normalize_if_needed(query, target_sum=cfg["normalize_target_sum"])

    label_cols = _detect_label_columns(reference, label_columns)
    if "X_umap" in query.obsm:
        query.obsm["X_umap_query"] = np.asarray(query.obsm["X_umap"])

    projected_cols: list[str] = []
    ingested = False
    for label_col in label_cols:
        projected_col = f"projected_{label_col}"
        logger.info("Projecting label column via ingest: %s -> %s", label_col, projected_col)
        if not ingested:
            sc.tl.ingest(query, reference, obs=label_col)
            ingested = True
            query.obs[projected_col] = query.obs[label_col].astype("category")
        else:
            qry = query.copy()
            sc.tl.ingest(qry, reference, obs=label_col)
            query.obs[projected_col] = qry.obs[label_col].astype("category")
        projected_cols.append(projected_col)
    query.uns["projection_atlas"] = str(atlas_path)
    query.uns["projection_query"] = str(query_path)
    query.uns["projection_label_columns"] = label_cols

    logger.info("Saving projected query: %s", out_h5ad)
    query.write_h5ad(out_h5ad)

    summary_rows: list[dict[str, Any]] = []
    if cfg.get("save_plots", True):
        ref_umap = np.asarray(reference.obsm["X_umap"])
        qry_umap = np.asarray(query.obsm["X_umap"])

        primary = projected_cols[0]
        primary_ref = label_cols[0]
        palette = _palette_for_label(primary_ref)

        overlay_path = figures_dir / f"umap_overlay_{primary_ref}.png"
        _plot_overlay(
            ref_umap,
            qry_umap,
            query.obs[primary],
            overlay_path,
            title=f"Atlas + projected query ({primary_ref})",
            dpi=cfg["dpi"],
            palette=palette,
        )
        logger.info("Saved figure: %s", overlay_path)

        query_umap_path = figures_dir / f"umap_query_{primary_ref}.png"
        _save_umap_scatter(
            qry_umap,
            query.obs[primary],
            query_umap_path,
            title=f"Projected query ({primary_ref})",
            dpi=cfg["dpi"],
            palette=palette,
            figsize=(10, 6),
        )
        logger.info("Saved figure: %s", query_umap_path)

        ref_comp = _composition_table(reference.obs[primary_ref], "atlas")
        query_comp = _composition_table(query.obs[primary], "projected_query")
        comp_csv = tables_dir / f"composition_{primary_ref}.csv"
        pd.concat([ref_comp, query_comp], ignore_index=True).to_csv(comp_csv, index=False)
        comp_plot = figures_dir / f"composition_{primary_ref}.png"
        _plot_composition(ref_comp, query_comp, comp_plot, dpi=cfg["dpi"])
        logger.info("Saved composition plot: %s", comp_plot)

        batch_key = query_batch_key if query_batch_key in query.obs.columns else None
        if batch_key:
            import matplotlib.pyplot as plt

            cats = query.obs[batch_key].astype(str).unique().tolist()
            n_cols = min(4, max(1, len(cats)))
            n_rows = int(np.ceil(len(cats) / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 4 * n_rows))
            axes = np.atleast_1d(axes).ravel()
            for ax, cat in zip(axes, cats, strict=False):
                mask = query.obs[batch_key].astype(str) == cat
                sub = qry_umap[mask]
                sub_lbl = query.obs[primary][mask]
                for label in pd.Categorical(sub_lbl).categories:
                    m = sub_lbl.astype(str) == str(label)
                    ax.scatter(sub[m, 0], sub[m, 1], s=5, alpha=0.7, rasterized=True, label=str(label))
                ax.set_title(f"{batch_key}={cat}")
            for ax in axes[len(cats):]:
                ax.axis("off")
            fig.suptitle(f"Projected query split by {batch_key}")
            split_path = figures_dir / f"umap_split_{batch_key}_{primary_ref}.png"
            fig.tight_layout()
            split_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(split_path, dpi=cfg["dpi"], bbox_inches="tight")
            plt.close(fig)
            logger.info("Saved split UMAP: %s", split_path)

        truth_col = query_label_column
        if truth_col and truth_col in query.obs.columns:
            conf_path = figures_dir / f"confusion_{truth_col}_vs_{primary_ref}.png"
            conf_df = _plot_confusion(query.obs[truth_col], query.obs[primary], conf_path, dpi=cfg["dpi"])
            conf_csv = tables_dir / f"confusion_{truth_col}_vs_{primary_ref}.csv"
            conf_df.to_csv(conf_csv)
            logger.info("Saved confusion matrix: %s", conf_path)

        for col, ref_col in zip(projected_cols, label_cols, strict=True):
            summary_rows.append(
                {
                    "projected_column": col,
                    "reference_column": ref_col,
                    "n_query_cells": int(query.n_obs),
                    "n_labels": int(query.obs[col].nunique()),
                }
            )

    summary_df = pd.DataFrame(summary_rows or [{"projected_column": projected_cols[0], "reference_column": label_cols[0]}])
    summary_csv = tables_dir / "projection_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    meta = {
        "atlas_h5ad": str(atlas_path),
        "query_h5ad": str(query_path),
        "projected_h5ad": str(out_h5ad),
        "label_columns": label_cols,
        "projected_columns": projected_cols,
        "n_shared_genes": int(reference.n_vars),
    }
    meta_json.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    logger.info("Projection completed. Outputs under: %s", out_root)
    return out_h5ad
