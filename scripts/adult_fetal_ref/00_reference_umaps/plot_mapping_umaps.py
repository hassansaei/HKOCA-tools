#!/usr/bin/env python
"""Plot organoid-on-reference UMAPs (3 figures per adult/fetal).

Level_2, transfer confidence, transferred labels.
Presence UMAP lives in 01_presence_coverage/plot_presence.py.
"""
from __future__ import annotations

import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib_config import ORG_L2_COLORS, PATHS, ROOT, WORKSPACE, apply_fig_style, read_h5_obs_col  # noqa: E402

apply_fig_style()

OUT = ROOT / "00_reference_umaps" / "full_projection"
HCKA_FULL = WORKSPACE / "human_kidney_reference_projection" / "Single_cell_human_kidney_cell_atlas.h5ad"
L2_ORDER = list(ORG_L2_COLORS.keys())


def style(ax, title: str) -> None:
    ax.set_title(title, fontsize=16, pad=8)
    ax.set_xlabel("UMAP1", fontsize=14)
    ax.set_ylabel("UMAP2", fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def freeze(ax, coords: np.ndarray, pad: float = 0.8) -> None:
    ax.set_xlim(float(coords[:, 0].min()) - pad, float(coords[:, 0].max()) + pad)
    ax.set_ylim(float(coords[:, 1].min()) - pad, float(coords[:, 1].max()) + pad)
    ax.set_aspect("equal", adjustable="box")


def draw_ref(ax, ref_umap: np.ndarray) -> None:
    ax.hexbin(
        ref_umap[:, 0], ref_umap[:, 1],
        gridsize=220, cmap="Greys", mincnt=1,
        linewidths=0, alpha=0.45, zorder=1, rasterized=True,
    )


def load_ref_umap(ref: str) -> np.ndarray:
    if ref == "adult":
        with h5py.File(HCKA_FULL, "r") as f:
            barcodes = read_h5_obs_col(f, "_index")
            disease = read_h5_obs_col(f, "disease")
            site = np.array(
                [str(x).strip().lower() for x in read_h5_obs_col(f, "sampled_site_condition")],
                dtype=object,
            )
            umap_full = np.asarray(f["obsm"]["X_scANVI_umap"][:], dtype=np.float32)
        keep = (disease == "normal") & np.isin(site, ["healthy", "adjacent"])
        ref_umap, ref_bc = umap_full[keep], barcodes[keep]
        with h5py.File(PATHS["adult_meta_broad"], "r") as f:
            lat_bc = read_h5_obs_col(f, "_index")
        pos = {b: i for i, b in enumerate(ref_bc)}
        return ref_umap[np.array([pos[b] for b in lat_bc], dtype=np.int64)]

    with h5py.File(PATHS["fetal_harmony"], "r") as fh:
        harm_bc = read_h5_obs_col(fh, "_index")
        harm_umap = np.asarray(fh["obsm"]["X_umap"][:], dtype=np.float32)
    with h5py.File(PATHS["fetal_umap"], "r") as fs:
        scanvi_ibc = read_h5_obs_col(fs, "integrated_barcode")
    harm_pos = {b: i for i, b in enumerate(harm_bc)}
    return harm_umap[np.array([harm_pos[b] for b in scanvi_ibc], dtype=np.int64)]


def plot_level2(ref_umap, org_umap, org_l2, out_png, ref: str) -> None:
    fig, ax = plt.subplots(figsize=(11, 8.5))
    draw_ref(ax, ref_umap)
    vc = pd.Series(org_l2).value_counts()
    order = [l for l in L2_ORDER if l in vc.index] + [l for l in vc.index if l not in L2_ORDER]
    for lab in order:
        m = org_l2 == lab
        ax.scatter(
            org_umap[m, 0], org_umap[m, 1],
            s=1.5, c=ORG_L2_COLORS.get(lab, "#888888"),
            alpha=0.75, linewidths=0, rasterized=True, zorder=2,
        )
    handles = [
        plt.Line2D(
            [0], [0], marker="o", color="w",
            markerfacecolor=ORG_L2_COLORS.get(lab, "#888888"), markersize=9,
            label=f"{lab.replace('_', ' ')} ({int(vc[lab]):,})",
        )
        for lab in order
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=11)
    style(ax, f"Organoids on {ref} — Level_2")
    freeze(ax, ref_umap)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def plot_confidence(ref_umap, org_umap, conf, out_png, ref: str) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    draw_ref(ax, ref_umap)
    conf = np.asarray(conf, dtype=np.float64)
    vmax = max(float(np.nanpercentile(conf, 98)), 1e-6)
    order = np.argsort(conf)
    sc = ax.scatter(
        org_umap[order, 0], org_umap[order, 1],
        s=1.2, c=np.clip(conf[order], 0, vmax),
        cmap="RdBu_r", vmin=0.0, vmax=vmax,
        alpha=0.7, linewidths=0, rasterized=True, zorder=2,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Transfer confidence", fontsize=13)
    style(ax, f"Mapping confidence on {ref}")
    freeze(ax, ref_umap)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def plot_transferred(ref_umap, org_umap, labels, out_png, ref: str) -> None:
    fig, ax = plt.subplots(figsize=(11, 8.5))
    draw_ref(ax, ref_umap)
    tv = pd.Series(labels).value_counts()
    base = [
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
        "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77",
        "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
    ]
    cols = {lab: base[i % len(base)] for i, lab in enumerate(tv.index)}
    for lab in reversed(list(tv.index)):
        m = labels == lab
        ax.scatter(
            org_umap[m, 0], org_umap[m, 1],
            s=1.2, c=cols[lab], alpha=0.7, linewidths=0, rasterized=True, zorder=2,
        )
    handles = [
        plt.Line2D(
            [0], [0], marker="o", color="w", markerfacecolor=cols[lab], markersize=8,
            label=f"{str(lab).replace('_', ' ')} ({int(tv[lab]):,})",
        )
        for lab in list(tv.index)[:22]
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=11)
    style(ax, f"Transferred labels on {ref}")
    freeze(ax, ref_umap)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def replot(ref: str) -> None:
    tab = OUT / ref / "tables"
    fig = OUT / ref / "figures"
    fig.mkdir(parents=True, exist_ok=True)
    # strip any leftover presence PNGs (owned by 01_presence_coverage)
    for stale in fig.glob(f"{ref}_*presence*.png"):
        stale.unlink()
        print(f"  removed leftover {stale.name}", flush=True)
    print(f"=== Mapping {ref} (3 figures) ===", flush=True)

    ref_umap = load_ref_umap(ref)
    org_umap = np.load(tab / "organoid_umap_soft_viz.npy")
    scores = pd.read_csv(tab / "transferred_labels.csv", index_col=0)
    org_l2 = scores["Level_2"].astype(str).to_numpy()
    conf = scores["best_score"].to_numpy(float)
    labels = scores["best_label"].astype(str).to_numpy()

    plot_level2(ref_umap, org_umap, org_l2, fig / f"{ref}_Level2.png", ref)
    plot_confidence(ref_umap, org_umap, conf, fig / f"{ref}_confidence.png", ref)
    plot_transferred(ref_umap, org_umap, labels, fig / f"{ref}_transferred_labels.png", ref)


def main() -> None:
    replot("adult")
    replot("fetal")
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
