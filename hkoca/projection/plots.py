"""Projection figures: atlas UMAP overlay, composition, uncertainty."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from hkoca.annotation.colors import lookup_color, palette_for_obs_key

LEVEL2_ORDER = [
    "undifferentiated_pluripotent_cells",
    "Off_target_cells",
    "Nephron_progenitor_cells",
    "Proliferative_tubule",
    "Podocytes",
    "Proximal_tubule",
    "LOH",
    "Distal_tubule",
    "Ureteric_Bud",
    "Collecting_duct",
    "Stromal_progenitors",
    "Proliferative_stroma",
    "Mesangial_cells",
    "Pericytes",
    "Interstitial_Fibroblasts",
    "Activated Fibroblasts",
    "Early_endothelium",
    "Late_endothelium",
]

SIM_CMAP = LinearSegmentedColormap.from_list(
    "atlas_query_sim", ["#C62828", "#F57C00", "#FDD835", "#42A5F5", "#0D47A1"]
)


def _palette(labels: pd.Series, obs_key: str) -> dict[str, str]:
    base = palette_for_obs_key(obs_key) or palette_for_obs_key("Level_3") or {}
    out = dict(base)
    for lab in pd.Series(labels).astype(str).unique():
        if lab not in out:
            out[lab] = lookup_color(lab, base)
    return out


def _style_umap(ax, title: str) -> None:
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=11, pad=8)


def knn_project_onto_atlas_umap(
    atlas_lat: np.ndarray,
    atlas_umap: np.ndarray,
    query_lat: np.ndarray,
    k: int = 30,
) -> tuple[np.ndarray, np.ndarray, float]:
    from sklearn.neighbors import NearestNeighbors

    k = min(int(k), atlas_lat.shape[0] - 1)
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
    nn.fit(atlas_lat)
    atlas_d, _ = nn.kneighbors(atlas_lat)
    d_ref = float(np.percentile(atlas_d[:, min(4, k - 1)], 95))
    d_ref = max(d_ref, 1e-6)
    q_d, q_idx = nn.kneighbors(query_lat)
    w = 1.0 / (q_d + 1e-6)
    w = w / w.sum(axis=1, keepdims=True)
    q_umap = (w[:, :, None] * atlas_umap[q_idx]).sum(axis=1)
    sim = np.clip(1.0 - q_d[:, 0] / d_ref, 0.0, 1.0).astype(np.float32)
    return q_umap.astype(np.float32), sim, d_ref


def subsample_rows(n: int, n_max: int, seed: int) -> np.ndarray:
    take = min(int(n_max), int(n))
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, size=take, replace=False))


def composition_table(labels: pd.Series, name: str) -> pd.DataFrame:
    counts = labels.astype(str).value_counts()
    total = max(int(counts.sum()), 1)
    return pd.DataFrame(
        {
            "dataset": name,
            "label": counts.index.astype(str),
            "n_cells": counts.values,
            "fraction": (counts.values / total).round(4),
        }
    )


def plot_overlay(
    ref_umap: np.ndarray,
    query_umap: np.ndarray,
    query_labels: pd.Series,
    out_path: Path,
    *,
    title: str,
    dpi: int,
    palette_key: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["axes.grid"] = False
    palette = _palette(query_labels, palette_key)
    fig, ax = plt.subplots(figsize=(12.5, 7))
    ax.grid(False)
    ax.scatter(ref_umap[:, 0], ref_umap[:, 1], s=2, c="#DDDDDD", rasterized=True, alpha=0.35, label="Atlas")
    lbl = query_labels.astype(str)
    cats = [c for c in pd.Categorical(lbl).categories if (lbl == c).any()]
    for cat in cats:
        mask = lbl == cat
        ax.scatter(
            query_umap[mask, 0],
            query_umap[mask, 1],
            s=6,
            c=[palette.get(cat, "#999999")],
            rasterized=True,
            alpha=0.75,
            label=cat,
        )
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0, 0.82, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_query_categorical(
    atlas_umap_bg: np.ndarray,
    q_umap: np.ndarray,
    labels: np.ndarray,
    out_path: Path,
    *,
    title: str,
    dpi: int,
    palette_key: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    palette = _palette(pd.Series(labels), palette_key)
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.grid(False)
    ax.scatter(atlas_umap_bg[:, 0], atlas_umap_bg[:, 1], s=0.12, c="#DDDDDD", alpha=0.35, rasterized=True, linewidths=0)
    order = [x for x in LEVEL2_ORDER if x in set(labels)] + sorted(set(labels) - set(LEVEL2_ORDER))
    for lv in order:
        m = labels == lv
        if not m.any():
            continue
        ax.scatter(
            q_umap[m, 0],
            q_umap[m, 1],
            s=1.5,
            c=[palette.get(lv, "#999999")],
            alpha=0.95,
            rasterized=True,
            linewidths=0,
            label=str(lv).replace("_", " "),
        )
    _style_umap(ax, title)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, markerscale=4)
    fig.tight_layout(rect=(0, 0, 0.78, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_similarity(
    atlas_umap_bg: np.ndarray,
    q_umap: np.ndarray,
    sim: np.ndarray,
    d_ref: float,
    out_path: Path,
    *,
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.grid(False)
    ax.scatter(atlas_umap_bg[:, 0], atlas_umap_bg[:, 1], s=0.12, c="#DDDDDD", alpha=0.35, rasterized=True, linewidths=0)
    sc = ax.scatter(
        q_umap[:, 0],
        q_umap[:, 1],
        s=1.2,
        c=sim,
        cmap=SIM_CMAP,
        vmin=0,
        vmax=1,
        alpha=0.9,
        rasterized=True,
        linewidths=0,
    )
    _style_umap(
        ax,
        "Query on atlas UMAP · latent similarity\n"
        f"(blue = close to atlas, red = farther; 95% atlas NN dist = {d_ref:.3f})",
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Atlas similarity score", fontsize=9)
    txt = (
        f"Mean similarity: {float(sim.mean()) * 100:.1f}%\n"
        f"Median: {float(np.median(sim)) * 100:.1f}%\n"
        f">=50%: {float((sim >= 0.50).mean()) * 100:.1f}% of query cells\n"
        f">=70%: {float((sim >= 0.70).mean()) * 100:.1f}%\n"
        f">=85%: {float((sim >= 0.85).mean()) * 100:.1f}%"
    )
    ax.text(
        0.02,
        0.02,
        txt,
        transform=ax.transAxes,
        fontsize=8.5,
        va="bottom",
        ha="left",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.92, "edgecolor": "#888888"},
    )
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#DDDDDD", markersize=7, label="Atlas (background)"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#0D47A1", markersize=7, label="Query (high sim)"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor="#C62828", markersize=7, label="Query (low sim)"),
        ],
        loc="upper right",
        frameon=False,
        fontsize=8,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_uncertainty(
    atlas_umap_bg: np.ndarray,
    q_umap: np.ndarray,
    uncert: np.ndarray,
    out_path: Path,
    *,
    dpi: int,
    title: str = "Query on atlas UMAP · Level_3_uncert",
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.grid(False)
    ax.scatter(atlas_umap_bg[:, 0], atlas_umap_bg[:, 1], s=0.12, c="#DDDDDD", alpha=0.35, rasterized=True, linewidths=0)
    vmax = max(0.5, float(np.nanpercentile(uncert, 99)))
    sc = ax.scatter(
        q_umap[:, 0],
        q_umap[:, 1],
        s=1.2,
        c=uncert,
        cmap="magma",
        vmin=0,
        vmax=vmax,
        alpha=0.9,
        rasterized=True,
        linewidths=0,
    )
    _style_umap(ax, title)
    fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02, label="Level_3_uncert")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_composition(ref_comp: pd.DataFrame, query_comp: pd.DataFrame, out_path: Path, *, dpi: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    merged = pd.concat([ref_comp, query_comp], ignore_index=True)
    labels = sorted(merged["label"].unique(), key=str)
    datasets = merged["dataset"].unique().tolist()
    x = np.arange(len(labels))
    width = 0.35 if len(datasets) == 2 else 0.8 / max(len(datasets), 1)
    fig, ax = plt.subplots(figsize=(max(10, len(labels) * 0.45), 5))
    ax.grid(False)
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


def plot_composition_by_sample(
    obs: pd.DataFrame,
    out_path: Path,
    *,
    dpi: int,
    batch_key: str = "sample_id",
    label_col: str = "Level_2_pred",
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if batch_key not in obs.columns or label_col not in obs.columns:
        return
    palette = _palette(obs[label_col], label_col)
    comp = obs.groupby([batch_key, label_col], observed=True).size().unstack(fill_value=0).astype(float)
    if comp.empty:
        return
    comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100.0
    order = [x for x in LEVEL2_ORDER if x in comp_pct.columns] + [c for c in comp_pct.columns if c not in LEVEL2_ORDER]
    comp_pct = comp_pct[order]
    bottom = np.zeros(len(comp_pct))
    x = np.arange(len(comp_pct))
    fig, ax = plt.subplots(figsize=(max(8, len(comp_pct) * 0.7), 5))
    ax.grid(False)
    for col in comp_pct.columns:
        vals = comp_pct[col].to_numpy()
        ax.bar(x, vals, bottom=bottom, color=palette.get(str(col), "#999999"), width=0.75, label=str(col).replace("_", " "))
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(comp_pct.index.astype(str), rotation=30, ha="right")
    ax.set_ylabel("Fraction of cells (%)")
    ax.set_title(f"Query {label_col} composition by {batch_key}")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7)
    fig.tight_layout(rect=(0, 0, 0.78, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_confusion(truth: pd.Series, predicted: pd.Series, out_path: Path, *, dpi: int) -> pd.DataFrame:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.crosstab(truth.astype(str), predicted.astype(str))
    fig, ax = plt.subplots(figsize=(max(6, df.shape[1] * 0.5), max(5, df.shape[0] * 0.5)))
    ax.grid(False)
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


def plot_joint_umap_panels(lat, query_uncert: np.ndarray, out_path: Path, *, dpi: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    lat = lat.copy()
    lat.obs["uncert_plot"] = np.nan
    q_mask = lat.obs["dataset_role"].astype(str).to_numpy() == "query"
    lat.obs.loc[q_mask, "uncert_plot"] = np.asarray(query_uncert, dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax in axes:
        ax.grid(False)
    sc.pl.umap(lat, color="dataset_role", ax=axes[0], show=False, frameon=False, title="Atlas vs query")
    sc.pl.umap(
        lat,
        color="Level_2_plot",
        ax=axes[1],
        show=False,
        frameon=False,
        title="Level_2 (atlas / pred)",
        legend_loc="right margin",
    )
    sc.pl.umap(
        lat,
        color="uncert_plot",
        ax=axes[2],
        show=False,
        frameon=False,
        title="Query label uncertainty",
        cmap="magma",
        vmin=0,
        vmax=1,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
