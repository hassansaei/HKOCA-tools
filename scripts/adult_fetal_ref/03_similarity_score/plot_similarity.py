#!/usr/bin/env python
"""Organoid vs primary transcriptional similarity (Spearman on log2-CPM).

Pseudobulk log2-CPM; HVG + TF gene sets. Matched niches use equal-weight
means across primary subtypes. Writes HVG/TF heatmaps, protocol bars,
and fetal-vs-adult matched fidelity.
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
import seaborn as sns
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib_config import (  # noqa: E402
    DIRS,
    PATHS,
    PROTOCOL_COLORS,
    PROTOCOLS_3,
    REF_COLORS,
    apply_fig_style,
    build_adult_fine_to_broad,
    ensure_dirs,
    load_human_tfs,
    merge_protocol,
    read_h5_obs_col,
)

apply_fig_style()

# All organoid Level_2 (ordered); Off-target / UPC included on heatmaps only
ORG_L2_ALL = [
    "Nephron_progenitor_cells",
    "Proliferative_nephron",
    "Podocytes",
    "PEC-like",
    "Proximal_tubule",
    "Distal_tubule",
    "LOH",
    "Collecting_duct",
    "Proliferative_stroma",
    "S1",
    "S2",
    "S3",
    "Endothelial_cells",
    "Off_target_cells",
    "undifferentiated_pluripotent_cells",
]

# Organoid Level_2 → adult broad_annotation (from reference)
ADULT_PAIRS_BROAD = {
    "Distal_tubule": ["DCT"],
    "LOH": ["TAL", "ATL", "DTL"],
    "Collecting_duct": ["PC", "IC-A", "IC-B", "CNT", "Papillary epithelium"],
    "PEC-like": ["PEC"],
    "Podocytes": ["Podocytes"],
    "Proliferative_stroma": ["Stroma", "Mesangial"],
    "Proximal_tubule": ["PT"],
    "S1": ["Stroma", "Mesangial"],
    "S2": ["Stroma", "Mesangial"],
    "S3": ["Stroma", "Mesangial"],
    "Endothelial_cells": ["EC"],
}

# Organoid Level_2 → fetal harmonized_celltype
FETAL_PAIRS = {
    "Distal_tubule": ["Distal_tubule_LOH", "Connecting_tubule"],
    "LOH": ["Distal_tubule_LOH"],
    "Collecting_duct": ["Ureteric_epithelium"],
    "Nephron_progenitor_cells": ["NPC"],
    "PEC-like": ["PEC"],
    "Podocytes": ["Podocytes", "Early_podocyte"],
    "Proliferative_nephron": ["NPC", "PTA", "Renal_vesicle_CSB", "SSB"],
    "Proliferative_stroma": ["Fibroblast", "Stroma_progenitor", "Myofibroblast"],
    "Proximal_tubule": ["Proximal_tubule"],
    "S1": ["Fibroblast", "Stroma_progenitor", "Myofibroblast"],
    "S2": ["Fibroblast", "Stroma_progenitor", "Myofibroblast"],
    "S3": ["Fibroblast", "Stroma_progenitor", "Myofibroblast"],
    "Endothelial_cells": ["Endothelial"],
}


def read_var_genes(f: h5py.File) -> np.ndarray:
    k = "_index" if "_index" in f["var"] else "gene_symbol"
    arr = f["var"][k][:]
    return np.array([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr])


def sum_counts_by_group(
    h5_path: Path,
    groups: np.ndarray,
    keep_mask: np.ndarray,
    chunk_size: int = 30000,
) -> tuple[pd.DataFrame, pd.Series]:
    """Return (genes × groups counts, n_cells per group)."""
    keep_idx = np.flatnonzero(keep_mask)
    unique_g, inv = np.unique(groups[keep_idx], return_inverse=True)
    n_per = pd.Series(np.bincount(inv), index=[str(x) for x in unique_g], dtype=int)
    with h5py.File(h5_path, "r") as f:
        genes = read_var_genes(f)
        n_genes = len(genes)
        g = f["X"]
        indptr = g["indptr"][:]
        n_obs = len(indptr) - 1
        cell_to_grp = np.full(n_obs, -1, dtype=np.int32)
        cell_to_grp[keep_idx] = inv.astype(np.int32)
        mat = np.zeros((len(unique_g), n_genes), dtype=np.float64)
        ind_ds, dat_ds = g["indices"], g["data"]
        for start in range(0, n_obs, chunk_size):
            end = min(n_obs, start + chunk_size)
            gc = cell_to_grp[start:end]
            loc_keep = np.flatnonzero(gc >= 0)
            if len(loc_keep) == 0:
                continue
            a, b = int(indptr[start]), int(indptr[end])
            ind_c, dat_c = ind_ds[a:b], dat_ds[a:b]
            ip_c = indptr[start : end + 1] - a
            for loc in loc_keep:
                r0, r1 = int(ip_c[loc]), int(ip_c[loc + 1])
                if r1 > r0:
                    np.add.at(mat[gc[loc]], ind_c[r0:r1], dat_c[r0:r1])
    return pd.DataFrame(mat.T, index=genes, columns=[str(x) for x in unique_g]), n_per


def counts_to_log2_cpm(df: pd.DataFrame) -> pd.DataFrame:
    lib = df.sum(axis=0).replace(0, np.nan)
    return np.log2(df.div(lib, axis=1) * 1e6 + 1.0)


def select_hvgs(ref_logcpm: pd.DataFrame, n_top: int = 2000) -> list[str]:
    var = ref_logcpm.var(axis=1)
    mean = ref_logcpm.mean(axis=1)
    ok = (ref_logcpm.max(axis=1) > 0.5) & np.isfinite(var) & (mean > 0)
    return var[ok].sort_values(ascending=False).head(n_top).index.astype(str).tolist()


def spearman(a: pd.Series, b: pd.Series, genes: list[str]) -> float:
    va, vb = a.loc[genes].to_numpy(float), b.loc[genes].to_numpy(float)
    m = np.isfinite(va) & np.isfinite(vb)
    if m.sum() < 20:
        return np.nan
    r, _ = spearmanr(va[m], vb[m])
    return float(r) if np.isfinite(r) else np.nan


# Display labels: expand abbreviations; data keys stay unchanged
DISPLAY_NAMES = {
    # adult broad
    "PT": "Proximal tubule",
    "DCT": "Distal convoluted tubule",
    "TAL": "Thick ascending limb",
    "ATL": "Thin ascending limb",
    "DTL": "Thin descending limb",
    "PC": "Principal cells",
    "IC-A": "Type A intercalated cells",
    "IC-B": "Type B intercalated cells",
    "CNT": "Connecting tubule",
    "PEC": "Parietal epithelial cells",
    "EC": "Endothelial cells",
    "Papillary epithelium": "Papillary epithelium",
    "Stroma": "Stroma",
    "Mesangial": "Mesangial",
    "Podocytes": "Podocytes",
    # fetal harmonized
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
    "Myofibroblast": "Myofibroblast",
    "Fibroblast": "Fibroblast",
    "Endothelial": "Endothelial",
}


def pretty(name: str) -> str:
    s = str(name)
    if s in DISPLAY_NAMES:
        return DISPLAY_NAMES[s]
    return s.replace("_", " ")


def matched_profile(ref_lc: pd.DataFrame, labels: list[str]) -> pd.Series | None:
    """Equal-weight mean of log2-CPM across primary subtypes (avoids size dominance)."""
    present = [l for l in labels if l in ref_lc.columns]
    if not present:
        return None
    return ref_lc[present].mean(axis=1)


def ref_labels_for_plot(ref_name: str, ref_lc: pd.DataFrame) -> list[str]:
    """Primary labels for heatmap columns: adult broad; fetal harmonized."""
    pairs = ADULT_PAIRS_BROAD if ref_name == "adult" else FETAL_PAIRS
    labs = []
    for ct in ORG_L2_ALL:
        for l in pairs.get(ct, []):
            if l in ref_lc.columns and l not in labs:
                labs.append(l)
    return labs


def matched_labels(ref_name: str, cell_type: str) -> list[str]:
    if ref_name == "adult":
        return ADULT_PAIRS_BROAD.get(cell_type, [])
    return FETAL_PAIRS.get(cell_type, [])


def plot_heatmap(mat: pd.DataFrame, title: str, xlab: str, out_png: Path) -> None:
    show = mat.copy()
    show.index = [pretty(i) for i in show.index]
    show.columns = [pretty(c) for c in show.columns]
    show = show.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if show.empty:
        return
    n_x, n_y = show.shape[1], show.shape[0]
    fig_w = max(9.0, 0.85 * n_x + 4.0)
    fig_h = max(6.0, 0.48 * n_y + 2.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        show.astype(float),
        ax=ax,
        cmap="RdBu_r",
        center=0.55,
        vmin=0.2,
        vmax=0.95,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 9},
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Spearman ρ", "shrink": 0.75},
        xticklabels=True,
        yticklabels=True,
    )
    ax.set_title(title, fontsize=16, pad=10)
    ax.set_xlabel(xlab, fontsize=14, labelpad=8)
    ax.set_ylabel("Organoid Level_2", fontsize=14)
    # Tick at column center; vertical labels sit directly under each column
    ax.set_xticks(np.arange(n_x) + 0.5)
    ax.set_xticklabels(list(show.columns), rotation=90, ha="center", va="top", fontsize=11)
    ax.set_yticks(np.arange(n_y) + 0.5)
    ax.set_yticklabels(list(show.index), rotation=0, va="center", fontsize=12)
    ax.tick_params(axis="both", which="major", length=0, pad=4)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white", pad_inches=0.25)
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def plot_matched_by_protocol(sim: pd.DataFrame, ref_name: str, out_png: Path) -> None:
    """Matched niche ρ (HVG) by protocol for each organoid type."""
    cts = (
        sim[sim["protocol"] == "All organoids"]
        .sort_values("spearman_HVG")
        ["cell_type"]
        .tolist()
    )
    if not cts:
        return
    y = np.arange(len(cts))
    fig, ax = plt.subplots(figsize=(9.5, max(4.5, 0.5 * len(cts) + 1.8)))
    # All organoids as light bars
    all_v = (
        sim[sim["protocol"] == "All organoids"]
        .set_index("cell_type")
        .reindex(cts)["spearman_HVG"]
    )
    ax.barh(y, all_v, height=0.55, color="#E0E0E0", edgecolor="#888888", label="All organoids")
    offsets = {"Morizane": 0.18, "Takasato": 0.0, "Freedman": -0.18}
    for proto in PROTOCOLS_3:
        sub = sim[sim["protocol"] == proto].set_index("cell_type").reindex(cts)
        ax.scatter(
            sub["spearman_HVG"],
            y + offsets[proto],
            s=90,
            color=PROTOCOL_COLORS[proto],
            edgecolors="white",
            linewidths=0.7,
            zorder=4,
            label=proto,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([pretty(c) for c in cts], fontsize=12)
    ax.set_xlabel("Spearman ρ (HVG vs matched primary)", fontsize=14)
    ax.set_title(f"Matched fidelity by protocol vs {ref_name}", fontsize=16)
    ax.set_xlim(0.3, 1.0)
    ax.legend(loc="lower right", frameon=False, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def plot_fetal_vs_adult(fetal: pd.DataFrame, adult: pd.DataFrame, out_png: Path) -> None:
    f = fetal[fetal["protocol"] == "All organoids"].set_index("cell_type")
    a = adult[adult["protocol"] == "All organoids"].set_index("cell_type")
    shared = [c for c in ORG_L2_ALL if c in f.index and c in a.index]
    if not shared:
        return
    y = np.arange(len(shared))
    fig, ax = plt.subplots(figsize=(9, max(4.5, 0.48 * len(shared) + 1.6)))
    ax.barh(
        y - 0.18,
        [f.loc[c, "spearman_HVG"] for c in shared],
        height=0.35,
        color=REF_COLORS["fetal"],
        label="vs fetal",
    )
    ax.barh(
        y + 0.18,
        [a.loc[c, "spearman_HVG"] for c in shared],
        height=0.35,
        color=REF_COLORS["adult"],
        label="vs adult",
    )
    ax.set_yticks(y)
    ax.set_yticklabels([pretty(c) for c in shared], fontsize=12)
    ax.set_xlabel("Spearman ρ (HVG, matched niche)", fontsize=14)
    ax.set_title("Matched fidelity: fetal vs adult", fontsize=16)
    ax.set_xlim(0.3, 1.0)
    ax.legend(loc="lower right", frameon=False, fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def run_for_reference(
    ref_name: str,
    org_all_lc: pd.DataFrame,
    org_proto_lc: dict[str, pd.DataFrame],
    org_n: pd.Series,
    ref_counts: pd.DataFrame,
    ref_n: pd.Series,
    tfs: list[str],
) -> pd.DataFrame:
    print(f"=== Similarity vs {ref_name} ===", flush=True)
    ref_lc = counts_to_log2_cpm(ref_counts)
    common = ref_lc.index.intersection(org_all_lc.index)
    hvgs = [g for g in select_hvgs(ref_lc.loc[common]) if g in common]
    tf_set = set(tfs)
    tf_genes = [
        g
        for g in common
        if g in tf_set and (ref_lc.loc[g].max() > 0.2 or org_all_lc.loc[g].max() > 0.2)
    ]
    print(f"  {ref_name}: {len(hvgs)} HVGs, {len(tf_genes)} TFs", flush=True)

    # --- A. Cross heatmap: all organoid Level_2 × primary labels ---
    xlab_labs = ref_labels_for_plot(ref_name, ref_lc)
    org_cts = [c for c in ORG_L2_ALL if c in org_all_lc.columns]
    mat_hvg = pd.DataFrame(index=org_cts, columns=xlab_labs, dtype=float)
    mat_tf = pd.DataFrame(index=org_cts, columns=xlab_labs, dtype=float)
    for oc in org_cts:
        for lab in xlab_labs:
            mat_hvg.loc[oc, lab] = spearman(org_all_lc[oc], ref_lc[lab], hvgs)
            mat_tf.loc[oc, lab] = spearman(org_all_lc[oc], ref_lc[lab], tf_genes)

    mat_hvg.to_csv(DIRS["03_tab"] / f"{ref_name}_similarity_matrix_HVG_actual.csv")
    mat_tf.to_csv(DIRS["03_tab"] / f"{ref_name}_similarity_matrix_TF_actual.csv")
    xlab = (
        f"{ref_name.capitalize()} broad annotation"
        if ref_name == "adult"
        else "Fetal harmonized cell type"
    )
    plot_heatmap(
        mat_hvg,
        f"HVG similarity vs {ref_name}",
        xlab,
        DIRS["03_fig"] / f"{ref_name}_similarity_heatmap_HVG.png",
    )
    plot_heatmap(
        mat_tf,
        f"TF similarity vs {ref_name}",
        xlab,
        DIRS["03_fig"] / f"{ref_name}_similarity_heatmap_TF.png",
    )

    # --- B. Matched niche by protocol (types with a defined primary match) ---
    rows = []
    for ct in ORG_L2_ALL:
        labs = matched_labels(ref_name, ct)
        matched = matched_profile(ref_lc, labs)
        if matched is None:
            continue
        n_ref = int(sum(int(ref_n.get(l, 0)) for l in labs if l in ref_n.index))
        for pname, olc in [("All organoids", org_all_lc)] + [
            (p, org_proto_lc[p]) for p in PROTOCOLS_3 if p in org_proto_lc
        ]:
            if ct not in olc.columns:
                continue
            rows.append(
                {
                    "reference": ref_name,
                    "cell_type": ct,
                    "protocol": pname,
                    "matched_primary_labels": ";".join([l for l in labs if l in ref_lc.columns]),
                    "n_org_cells": int(org_n.get(ct, 0)) if pname == "All organoids" else np.nan,
                    "n_ref_cells_matched": n_ref,
                    "spearman_HVG": spearman(olc[ct], matched, hvgs),
                    "spearman_TF": spearman(olc[ct], matched, tf_genes),
                    "n_hvgs": len(hvgs),
                    "n_tfs": len(tf_genes),
                }
            )
    sim = pd.DataFrame(rows)
    sim.to_csv(DIRS["03_tab"] / f"{ref_name}_similarity_matched_by_protocol.csv", index=False)
    plot_matched_by_protocol(
        sim, ref_name, DIRS["03_fig"] / f"{ref_name}_similarity_matched_by_protocol.png"
    )

    # cell counts for transparency
    n_tab = pd.DataFrame(
        {
            "label": list(ref_n.index) + [f"organoid::{c}" for c in org_cts],
            "n_cells": list(ref_n.values) + [int(org_n.get(c, 0)) for c in org_cts],
            "compartment": ["primary"] * len(ref_n) + ["organoid"] * len(org_cts),
        }
    )
    n_tab.to_csv(DIRS["03_tab"] / f"{ref_name}_profile_cell_counts.csv", index=False)
    return sim


def main() -> None:
    ensure_dirs()
    for p in DIRS["03_fig"].glob("*.png"):
        p.unlink()
    for p in DIRS["03_tab"].glob("*.csv"):
        p.unlink()

    tfs = load_human_tfs()
    with h5py.File(PATHS["organoid_counts"], "r") as f:
        org_l2 = read_h5_obs_col(f, "Level_2")
        org_proto = np.array(
            [merge_protocol(p) for p in read_h5_obs_col(f, "diff_protocol")], dtype=object
        )

    keep = np.isin(org_l2, ORG_L2_ALL)
    print(f"Organoid Level_2 included: {[c for c in ORG_L2_ALL if c in set(org_l2)]}", flush=True)
    org_all_c, org_n = sum_counts_by_group(PATHS["organoid_counts"], org_l2, keep)
    org_all_lc = counts_to_log2_cpm(org_all_c)

    keys = np.array([f"{p}||{c}" for p, c in zip(org_proto, org_l2)], dtype=object)
    org_p_c, _ = sum_counts_by_group(
        PATHS["organoid_counts"], keys, keep & np.isin(org_proto, PROTOCOLS_3)
    )
    org_p_lc = counts_to_log2_cpm(org_p_c)
    org_proto_lc = {}
    for p in PROTOCOLS_3:
        cols = [c for c in org_p_lc.columns if c.startswith(f"{p}||")]
        sub = org_p_lc[cols].copy()
        sub.columns = [c.split("||", 1)[1] for c in cols]
        org_proto_lc[p] = sub

    with h5py.File(PATHS["fetal_counts"], "r") as f:
        f_harm = read_h5_obs_col(f, "harmonized_celltype")
    f_c, f_n = sum_counts_by_group(
        PATHS["fetal_counts"], f_harm, np.ones(len(f_harm), dtype=bool)
    )

    # Adult: group by broad_annotation (not fine PTS1/2/3 …)
    fine_to_broad = build_adult_fine_to_broad()
    with h5py.File(PATHS["adult_counts"], "r") as f:
        a_fine = read_h5_obs_col(f, "fine_annotation")
    a_broad = np.array([fine_to_broad.get(str(x), "NA") for x in a_fine], dtype=object)
    a_keep = ~np.isin(a_broad, ["NA", "Unassigned"])
    a_c, a_n = sum_counts_by_group(PATHS["adult_counts"], a_broad, a_keep)
    print(f"Adult broad types used: {sorted(a_c.columns)}", flush=True)

    fetal_sim = run_for_reference(
        "fetal", org_all_lc, org_proto_lc, org_n, f_c, f_n, tfs
    )
    adult_sim = run_for_reference(
        "adult", org_all_lc, org_proto_lc, org_n, a_c, a_n, tfs
    )
    plot_fetal_vs_adult(fetal_sim, adult_sim, DIRS["03_fig"] / "fetal_vs_adult_matched_fidelity.png")
    print("Similarity complete.", flush=True)


if __name__ == "__main__":
    main()
