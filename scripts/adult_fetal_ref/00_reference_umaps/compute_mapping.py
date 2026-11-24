#!/usr/bin/env python
"""Compute wkNN organoid→reference mapping (heavy). Figures: use plot_mapping_umaps.py.

Builds wkNN, label transfer, and soft UMAP coords under
00_reference_umaps/full_projection/{adult,fetal}/tables/.

Presence scores / UMAP live in 01_presence_coverage/plot_presence.py.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mapping_data import load_adult, load_fetal  # noqa: E402
from wknn_mapping import get_wknn, transfer_labels  # noqa: E402
from lib_config import ORG_L2_COLORS, ROOT, apply_fig_style  # noqa: E402

apply_fig_style()

OUT_ROOT = ROOT / "00_reference_umaps" / "full_projection"

# wkNN defaults
WKNN_K = 100
WEIGHTING = "jaccard_square"


def style_umap(ax, title: str) -> None:
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel("UMAP1", fontsize=14)
    ax.set_ylabel("UMAP2", fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    for spn in ax.spines.values():
        spn.set_visible(False)


def freeze_limits(ax, coords: np.ndarray, pad: float = 0.8) -> None:
    ax.set_xlim(float(coords[:, 0].min()) - pad, float(coords[:, 0].max()) + pad)
    ax.set_ylim(float(coords[:, 1].min()) - pad, float(coords[:, 1].max()) + pad)
    ax.set_aspect("equal", adjustable="box")


def keep_topn_neighbors(wknn: sp.spmatrix, top_n: int) -> sp.csr_matrix:
    """Keep only the heaviest ``top_n`` wkNN edges per query cell (CSR)."""
    w = wknn.tocsr(copy=True)
    if top_n <= 0:
        return w
    data = w.data.copy()
    for i in range(w.shape[0]):
        start, end = w.indptr[i], w.indptr[i + 1]
        length = end - start
        if length <= top_n:
            continue
        row = data[start:end]
        keep = np.argpartition(row, -top_n)[-top_n:]
        mask = np.zeros(length, dtype=bool)
        mask[keep] = True
        row[~mask] = 0.0
        data[start:end] = row
    out = sp.csr_matrix((data, w.indices, w.indptr), shape=w.shape)
    out.eliminate_zeros()
    return out


def soft_umap_from_wknn(
    wknn: sp.spmatrix,
    ref_umap: np.ndarray,
    *,
    top_n: int | None = None,
    power: float = 1.0,
) -> np.ndarray:
    """Soft UMAP landing: row-normalized wkNN @ ref UMAP.

    Parameters (viz knobs; wkNN graph itself unchanged)
    ----------
    top_n
        If set, average only the heaviest ``top_n`` neighbors (wkNN weights).
        ``None`` = all k neighbors (exact vignette / matched-transcriptome softmax).
    power
        Raise edge weights to ``power`` before normalizing (``>1`` sharpens toward
        the strongest neighbors; ``1`` = plain softmax).
    """
    w = wknn.tocsr()
    if top_n is not None:
        w = keep_topn_neighbors(w, top_n=int(top_n))
    if power != 1.0:
        w = w.copy()
        w.data = np.power(np.maximum(w.data, 0.0), float(power))
    row_sum = np.asarray(w.sum(axis=1)).ravel()
    row_sum = np.maximum(row_sum, 1e-12)
    soft = sp.diags(1.0 / row_sum) @ w
    return np.asarray(soft @ ref_umap, dtype=np.float32)


# Figure landing: soft average, but sharper than plain k=100
# (plain k=100 is muddy). Label transfer / soft UMAP still use full wkNN.
VIZ_SOFT_TOP_N = 10  # average top-10 neighbors (None = all 100)
VIZ_SOFT_POWER = 2.0  # sharpen weights before softmax
MUDDY_L2 = {"Off_target_cells", "undifferentiated_pluripotent_cells"}
MUDDY_MAX_DRAW = 25_000


def plot_full_mapping(
    ref_umap: np.ndarray,
    org_umap: np.ndarray,
    org_l2: np.ndarray,
    conf: np.ndarray,
    transferred: np.ndarray,
    out_png: Path,
    ref_name: str,
    *,
    l2_colors: dict[str, str] | None = None,
    level2_title_suffix: str = "",
    write_transferred: bool = True,
) -> None:
    """Level_2 overlay + transfer confidence (sharp landing, readable rendering)."""
    colors = l2_colors or ORG_L2_COLORS
    n_all = len(org_l2)
    rng = np.random.default_rng(42)
    conf = np.asarray(conf, dtype=np.float64)
    vmax_conf = float(np.nanpercentile(conf, 98))
    vmax_conf = max(vmax_conf, 1e-6)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8.5))

    def _ref(ax):
        # Density backdrop (avoids gray-point fog); organoids stay as scatter
        ax.hexbin(
            ref_umap[:, 0],
            ref_umap[:, 1],
            gridsize=220,
            cmap="Greys",
            mincnt=1,
            linewidths=0,
            alpha=0.45,
            zorder=1,
            rasterized=True,
        )

    ax = axes[0]
    _ref(ax)
    vc = pd.Series(org_l2).value_counts()
    # Muddy types under (fainter); kidney types on top — same style as full_projection
    draw_order = [lab for lab in vc.index if lab in MUDDY_L2] + [
        lab for lab in vc.index if lab not in MUDDY_L2
    ]
    for lab in draw_order:
        m = org_l2 == lab
        if not np.any(m):
            continue
        muddy = lab in MUDDY_L2
        ax.scatter(
            org_umap[m, 0],
            org_umap[m, 1],
            s=0.8 if muddy else 1.6,
            c=colors.get(lab, "#888888"),
            alpha=0.35 if muddy else 0.7,
            linewidths=0,
            rasterized=True,
            zorder=2 if muddy else 3,
        )
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=colors.get(lab, "#888888"),
            markersize=8,
            label=f"{lab} ({int(vc[lab]):,})",
        )
        for lab in vc.index
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        fontsize=10,
    )
    suffix = f" {level2_title_suffix}" if level2_title_suffix else ""
    style_umap(
        ax,
        f"Organoids on {ref_name} — Level_2"
        f"{suffix} "
        f"(top-{VIZ_SOFT_TOP_N}, power={VIZ_SOFT_POWER:g}; n = {n_all:,})",
    )
    freeze_limits(ax, ref_umap)

    ax = axes[1]
    _ref(ax)
    # Subsample lowest half of scores so high-confidence structure is visible
    ord_all = np.argsort(conf)
    n_low = n_all // 2
    low_idx = ord_all[:n_low]
    if len(low_idx) > 60_000:
        low_idx = rng.choice(low_idx, 60_000, replace=False)
    high_idx = ord_all[n_low:]
    plot_idx = np.concatenate([np.sort(low_idx), high_idx])
    sc = ax.scatter(
        org_umap[plot_idx, 0],
        org_umap[plot_idx, 1],
        s=np.where(conf[plot_idx] >= np.median(conf), 1.8, 0.5),
        c=np.clip(conf[plot_idx], 0, vmax_conf),
        cmap="viridis",
        vmin=0.0,
        vmax=vmax_conf,
        alpha=0.75,
        linewidths=0,
        rasterized=True,
        zorder=2,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.038, pad=0.02)
    cbar.set_label(f"Transfer confidence (vmax = p98 = {vmax_conf:.2f})", fontsize=13)
    cbar.ax.tick_params(labelsize=12)
    style_umap(
        ax,
        f"Mapping confidence on {ref_name} "
        f"(soft landing top-{VIZ_SOFT_TOP_N}; n = {n_all:,})",
    )
    freeze_limits(ax, ref_umap)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)

    if not write_transferred:
        return

    # Second figure: transferred reference labels
    fig, ax = plt.subplots(figsize=(12, 9))
    _ref(ax)
    tv = pd.Series(transferred).value_counts()
    cmap = plt.get_cmap("tab20")
    lab_cols = {lab: cmap(i % 20) for i, lab in enumerate(tv.index)}
    # Draw rare labels last so they stay visible
    for lab in reversed(list(tv.index)):
        m = transferred == lab
        idx = np.flatnonzero(m)
        if len(idx) > 80_000:
            idx = rng.choice(idx, 80_000, replace=False)
        ax.scatter(
            org_umap[idx, 0],
            org_umap[idx, 1],
            s=1.2,
            c=[lab_cols[lab]],
            alpha=0.7,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=lab_cols[lab],
            markersize=8,
            label=f"{lab} ({int(tv[lab]):,})",
        )
        for lab in list(tv.index)[:25]
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        fontsize=10,
        title="Transferred (top 25)",
    )
    style_umap(ax, f"Transferred {ref_name} labels (n = {n_all:,})")
    freeze_limits(ax, ref_umap)
    out2 = out_png.parent / f"{ref_name}_transferred_labels.png"
    fig.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out2}", flush=True)


def run_atlas(
    ref_name: str,
    ref_lat: np.ndarray,
    ref_umap: np.ndarray,
    ref_labels: np.ndarray,
    org_lat: np.ndarray,
    org_l2: np.ndarray,
    label_key: str,
) -> None:
    out_dir = OUT_ROOT / ref_name
    fig_dir = out_dir / "figures"
    tab_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"=== Mapping: {ref_name} | ref={ref_lat.shape[0]:,} query={org_lat.shape[0]:,} ===",
        flush=True,
    )
    print(
        f"  get_wknn(k={WKNN_K}, weighting={WEIGHTING}, query2ref=True, ref2query=False)",
        flush=True,
    )
    wknn = get_wknn(
        ref=np.asarray(ref_lat, dtype=np.float32),
        query=np.asarray(org_lat, dtype=np.float32),
        k=WKNN_K,
        query2ref=True,
        ref2query=False,
        weighting_scheme=WEIGHTING,
    )
    # Ensure (n_query, n_ref)
    if wknn.shape[0] != org_lat.shape[0]:
        wknn = wknn.T.tocsr()
    print(f"  wkNN shape={wknn.shape} nnz={wknn.nnz:,}", flush=True)
    sp.save_npz(tab_dir / "wknn.npz", wknn)

    # Build minimal AnnData shells for transfer_labels
    import anndata as ad

    ref_ad = ad.AnnData(X=sp.csr_matrix((ref_lat.shape[0], 1)))
    ref_ad.obs[label_key] = pd.Categorical(ref_labels.astype(str))
    ref_ad.obs_names = [f"ref_{i}" for i in range(ref_lat.shape[0])]
    ref_ad.obsm["X_latent"] = np.asarray(ref_lat, dtype=np.float32)
    ref_ad.obsm["X_umap"] = np.asarray(ref_umap, dtype=np.float32)

    qry_ad = ad.AnnData(X=sp.csr_matrix((org_lat.shape[0], 1)))
    qry_ad.obs["Level_2"] = pd.Categorical(org_l2.astype(str))
    qry_ad.obs_names = [f"qry_{i}" for i in range(org_lat.shape[0])]
    qry_ad.obsm["X_latent"] = np.asarray(org_lat, dtype=np.float32)

    print("  transfer_labels…", flush=True)
    scores = transfer_labels(ref_ad, qry_ad, wknn, label_key=label_key)
    scores.insert(0, "Level_2", org_l2)
    scores.to_csv(tab_dir / "transferred_labels.csv")
    best_label = scores["best_label"].astype(str).to_numpy()
    best_score = scores["best_score"].to_numpy(dtype=np.float32)

    print("  soft UMAP landing (full k=100 softmax; saved)…", flush=True)
    org_umap_soft = soft_umap_from_wknn(wknn, ref_umap, top_n=None, power=1.0)
    np.save(tab_dir / "organoid_umap_soft.npy", org_umap_soft)
    print(
        f"  figure soft landing (wkNN weights; top_n={VIZ_SOFT_TOP_N}, "
        f"power={VIZ_SOFT_POWER})…",
        flush=True,
    )
    org_umap = soft_umap_from_wknn(
        wknn, ref_umap, top_n=VIZ_SOFT_TOP_N, power=VIZ_SOFT_POWER
    )
    np.save(tab_dir / "organoid_umap_soft_viz.npy", org_umap)

    plot_full_mapping(
        ref_umap,
        org_umap,
        org_l2,
        best_score,
        best_label,
        fig_dir / f"{ref_name}_Level2_confidence.png",
        ref_name,
    )

    # Summary tables
    land = (
        pd.Series(org_l2)
        .value_counts()
        .rename_axis("Level_2")
        .reset_index(name="n_cells")
    )
    land["pct"] = 100.0 * land["n_cells"] / len(org_l2)
    land.to_csv(tab_dir / "organoid_Level2_counts.csv", index=False)

    xfer = (
        pd.Series(best_label)
        .value_counts()
        .rename_axis(label_key)
        .reset_index(name="n_cells")
    )
    xfer["pct"] = 100.0 * xfer["n_cells"] / len(best_label)
    xfer["median_best_score"] = [
        float(np.median(best_score[best_label == lab])) for lab in xfer[label_key]
    ]
    xfer.to_csv(tab_dir / "transferred_label_counts.csv", index=False)

    meta = {
        "method": "wkNN soft mapping (vendored)",
        "note": "weighted kNN + soft UMAP landing + label transfer",
        "reference": ref_name,
        "n_ref": int(ref_lat.shape[0]),
        "n_query": int(org_lat.shape[0]),
        "latent": "X_scANVI (precomputed; no AtlasMapper retrain)",
        "wknn": {
            "k": WKNN_K,
            "weighting_scheme": WEIGHTING,
            "query2ref": True,
            "ref2query": False,
        },
        "label_key": label_key,
        "umap_landing_saved": "full-k softmax (organoid_umap_soft.npy)",
        "umap_landing_figures": {
            "method": "soft average of wkNN neighbors",
            "top_n": VIZ_SOFT_TOP_N,
            "power": VIZ_SOFT_POWER,
            "file": "organoid_umap_soft_viz.npy",
        },
        "plotted": "all query cells (muddy L2 fainter; confidence low-half subsampled)",
        "confidence_vmax": "98th percentile of best_score",
        "note": "presence scores are computed in 01_presence_coverage/plot_presence.py",
    }
    (tab_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"=== Done {ref_name} → {out_dir} ===", flush=True)


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description="Compute wkNN mapping tables (figures: plot_mapping_umaps.py)")
    ap.add_argument("--which", choices=["adult", "fetal", "both"], default="both")
    ap.add_argument("--viz-top-n", type=int, default=None)
    ap.add_argument("--viz-power", type=float, default=None)
    args = ap.parse_args()
    global VIZ_SOFT_TOP_N, VIZ_SOFT_POWER
    if args.viz_top_n is not None:
        VIZ_SOFT_TOP_N = args.viz_top_n
    if args.viz_power is not None:
        VIZ_SOFT_POWER = args.viz_power

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    if args.which in ("adult", "both"):
        print("Loading adult…", flush=True)
        ref_lat, ref_umap, ref_lab, org_lat, org_l2, key = load_adult()
        run_atlas("adult", ref_lat, ref_umap, ref_lab, org_lat, org_l2, key)
    if args.which in ("fetal", "both"):
        print("Loading fetal…", flush=True)
        ref_lat, ref_umap, ref_lab, org_lat, org_l2, key = load_fetal()
        run_atlas("fetal", ref_lat, ref_umap, ref_lab, org_lat, org_l2, key)
    print("Mapping complete. Replot figures with plot_mapping_umaps.py", flush=True)


if __name__ == "__main__":
    main()
