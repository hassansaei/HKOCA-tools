#!/usr/bin/env python
"""00 — Reference annotation UMAPs (after adult/fetal processing scripts).

Adult  — frozen HCKA X_scANVI_umap (healthy/adjacent), broad_annotation
Fetal  — frozen Harmony X_umap (weeks 9–17), harmonized_celltype

Upstream first: upstream/process_adult.py → upstream/process_fetal.py → this script.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib_config import (  # noqa: E402
    DIRS,
    PATHS,
    WORKSPACE,
    apply_fig_style,
    ensure_dirs,
    read_h5_obs_col,
)

apply_fig_style()

HCKA_FULL = (
    WORKSPACE
    / "human_kidney_reference_projection"
    / "Single_cell_human_kidney_cell_atlas.h5ad"
)

# High-contrast categorical colors (avoid near-shades within lineages).
ADULT_BROAD_COLORS = {
    "PT": "#E41A1C",
    "DCT": "#377EB8",
    "TAL": "#4DAF4A",
    "ATL": "#FF7F00",
    "DTL": "#984EA3",
    "CNT": "#A65628",
    "PC": "#F781BF",
    "IC-A": "#999999",
    "IC-B": "#66C2A5",
    "MD": "#FC8D62",
    "Podocytes": "#E78AC3",
    "PEC": "#8DA0CB",
    "Papillary epithelium": "#FFD92F",
    "Stroma": "#1B9E77",
    "Mesangial": "#D95F02",
    "EC": "#000000",
    "Neuron": "#E6AB02",
    "CD8 T cells": "#7570B3",
    "CD4 T cells": "#66A61E",
    "T cells": "#A6761D",
    "B cells": "#E7298A",
    "NK cells": "#1F78B4",
    "Macrophages": "#B2DF8A",
    "Monocytes": "#33A02C",
    "cDCs": "#FB9A99",
    "pDCs": "#CAB2D6",
    "Plasma cells": "#6A3D9A",
    "Mast cells": "#B15928",
    "Neutrophils": "#FFFF99",
    "RBC": "#FDBF6F",
    "Platelets": "#BDBDBD",
    "NA": "#D9D9D9",
}

FETAL_HARM_COLORS = {
    "NPC": "#E41A1C",
    "PTA": "#377EB8",
    "Renal_vesicle_CSB": "#4DAF4A",
    "SSB": "#FF7F00",
    "Podocytes": "#984EA3",
    "Early_podocyte": "#F781BF",
    "PEC": "#A6CEE3",
    "Proximal_tubule": "#A65628",
    "Distal_tubule_LOH": "#1F78B4",
    "Connecting_tubule": "#6A3D9A",
    "Ureteric_epithelium": "#66C2A5",
    "Fibroblast": "#E7298A",
    "Stroma_progenitor": "#66A61E",
    "Myofibroblast": "#E6AB02",
    "Mesangial": "#D95F02",
    "Endothelial": "#000000",
    "Immune": "#7570B3",
    "Megakaryocyte": "#BDBDBD",
    "Erythroid": "#FB9A99",
    "Neuron": "#FFFF33",
    "Unassigned": "#D9D9D9",
}

DISPLAY = {
    "PT": "Proximal tubule",
    "DCT": "Distal convoluted tubule",
    "TAL": "Thick ascending limb",
    "ATL": "Thin ascending limb",
    "DTL": "Thin descending limb",
    "CNT": "Connecting tubule",
    "PC": "Principal cells",
    "IC-A": "Type A intercalated cells",
    "IC-B": "Type B intercalated cells",
    "MD": "Macula densa",
    "PEC": "Parietal epithelial cells",
    "EC": "Endothelial cells",
    "NPC": "Nephron progenitor cells",
    "PTA": "Pretubular aggregate",
    "Renal_vesicle_CSB": "Renal vesicle",
    "SSB": "S-shaped body",
    "Early_podocyte": "Early podocyte",
    "Distal_tubule_LOH": "Distal tubule / LOH",
    "Connecting_tubule": "Connecting tubule",
    "Ureteric_epithelium": "Ureteric epithelium",
    "Proximal_tubule": "Proximal tubule",
    "Stroma_progenitor": "Stroma progenitor",
}


def pretty(lab: str) -> str:
    return DISPLAY.get(str(lab), str(lab).replace("_", " "))


def palette_for(labels: list[str], fixed: dict[str, str]) -> dict[str, str]:
    """Fixed distinct colors; leftover labels get spaced tab20 (rarely needed)."""
    out = {}
    leftovers = []
    for lab in labels:
        if lab in fixed:
            out[lab] = fixed[lab]
        else:
            leftovers.append(lab)
    if leftovers:
        cmap = plt.get_cmap("tab20")
        for i, lab in enumerate(leftovers):
            out[lab] = cmap((i * 3) % 20 / 20.0)
    for na in ("NA", "Unassigned"):
        if na in out:
            out[na] = "#D0D0D0"
    return out


def plot_broad_umap(
    coords: np.ndarray,
    labels: np.ndarray,
    out_png: Path,
    title: str,
    color_map: dict[str, str],
    *,
    point_size: float,
    alpha: float,
) -> pd.DataFrame:
    labels = np.asarray(labels, dtype=object)
    vc = pd.Series(labels).value_counts()
    order = [lab for lab in vc.index if lab not in ("NA", "Unassigned")]
    for na in ("NA", "Unassigned"):
        if na in vc.index:
            order.append(na)
    cmap = palette_for(order, color_map)

    fig, ax = plt.subplots(figsize=(12.5, 9.0))
    # draw large / grey last; rare types on top for visibility
    draw_order = list(reversed(order))
    for lab in draw_order:
        m = labels == lab
        ax.scatter(
            coords[m, 0],
            coords[m, 1],
            s=point_size,
            c=cmap[lab],
            alpha=0.25 if lab in ("NA", "Unassigned") else alpha,
            linewidths=0,
            rasterized=True,
        )
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel("UMAP1", fontsize=14)
    ax.set_ylabel("UMAP2", fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    pad = 0.8
    ax.set_xlim(float(coords[:, 0].min()) - pad, float(coords[:, 0].max()) + pad)
    ax.set_ylim(float(coords[:, 1].min()) - pad, float(coords[:, 1].max()) + pad)
    ax.set_aspect("equal", adjustable="box")

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=cmap[lab],
            markersize=9,
            label=f"{pretty(lab)} ({int(vc[lab]):,})",
        )
        for lab in order
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=11,
        ncol=1 if len(order) <= 22 else 2,
        handletextpad=0.4,
    )
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white", pad_inches=0.2)
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)

    tab = vc.rename_axis("broad_annotation").reset_index(name="n_cells")
    tab["pct_reference"] = 100.0 * tab["n_cells"] / len(labels)
    tab["display_name"] = tab["broad_annotation"].map(pretty)
    return tab


def run_adult() -> None:
    print("=== Adult broad annotation UMAP (X_scANVI_umap) ===", flush=True)
    with h5py.File(HCKA_FULL, "r") as f:
        barcodes = read_h5_obs_col(f, "_index")
        disease = read_h5_obs_col(f, "disease")
        site = np.array(
            [str(x).strip().lower() for x in read_h5_obs_col(f, "sampled_site_condition")],
            dtype=object,
        )
        broad = read_h5_obs_col(f, "broad_annotation")
        umap_full = np.asarray(f["obsm"]["X_scANVI_umap"][:], dtype=np.float32)

    keep = (disease == "normal") & np.isin(site, ["healthy", "adjacent"])
    ref_umap = umap_full[keep]
    ref_broad = broad[keep]
    ref_bc = barcodes[keep]

    with h5py.File(PATHS["adult_meta_broad"], "r") as f:
        lat_bc = read_h5_obs_col(f, "_index")
    pos = {b: i for i, b in enumerate(ref_bc)}
    order = np.array([pos[b] for b in lat_bc], dtype=np.int64)
    ref_umap = ref_umap[order]
    ref_broad = ref_broad[order]
    print(f"  n = {len(ref_broad):,} healthy/adjacent cells", flush=True)

    tab = plot_broad_umap(
        ref_umap,
        ref_broad,
        DIRS["00_fig"] / "adult_UMAP_broad_annotation.png",
        title=f"Adult kidney — broad annotation (scANVI UMAP, n = {len(ref_broad):,})",
        color_map=ADULT_BROAD_COLORS,
        point_size=0.55,
        alpha=0.40,
    )
    tab.to_csv(DIRS["00_tab"] / "adult_broad_annotation_counts.csv", index=False)
    (DIRS["00_tab"] / "adult_umap_source.json").write_text(
        json.dumps(
            {
                "umap_source": str(HCKA_FULL),
                "umap_key": "X_scANVI_umap",
                "label_key": "broad_annotation",
                "filter": "disease==normal & sampled_site_condition in {healthy,adjacent}",
                "n_cells": int(len(ref_broad)),
            },
            indent=2,
        )
    )


def run_fetal() -> None:
    print("=== Fetal harmonized annotation UMAP (Harmony) ===", flush=True)
    with h5py.File(PATHS["fetal_harmony"], "r") as fh:
        harm_umap = np.asarray(fh["obsm"]["X_umap"][:], dtype=np.float32)
        harm_broad = read_h5_obs_col(fh, "harmonized_celltype")
    print(f"  n = {len(harm_broad):,} cells", flush=True)

    tab = plot_broad_umap(
        harm_umap,
        harm_broad,
        DIRS["00_fig"] / "fetal_UMAP_broad_annotation.png",
        title=f"Fetal kidney — harmonized cell type (Harmony UMAP, n = {len(harm_broad):,})",
        color_map=FETAL_HARM_COLORS,
        point_size=2.2,
        alpha=0.55,
    )
    tab.to_csv(DIRS["00_tab"] / "fetal_broad_annotation_counts.csv", index=False)
    (DIRS["00_tab"] / "fetal_umap_source.json").write_text(
        json.dumps(
            {
                "umap_source": str(PATHS["fetal_harmony"]),
                "umap_key": "X_umap",
                "embedding": "X_pca_harmony",
                "label_key": "harmonized_celltype",
                "n_cells": int(len(harm_broad)),
            },
            indent=2,
        )
    )


def main() -> None:
    ensure_dirs()
    DIRS["00_fig"].mkdir(parents=True, exist_ok=True)
    ap = argparse.ArgumentParser(description="Regenerate adult/fetal broad UMAPs only")
    ap.add_argument("--which", choices=["adult", "fetal", "both"], default="both")
    args = ap.parse_args()
    if args.which in ("adult", "both"):
        run_adult()
    if args.which in ("fetal", "both"):
        run_fetal()
    print("Reference broad UMAPs complete.", flush=True)


if __name__ == "__main__":
    main()
