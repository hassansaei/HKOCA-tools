#!/usr/bin/env python
"""Presence coverage — compute + plot 3 figures per adult/fetal reference.

1. Composition heatmap (confident wkNN labels × protocol)
2. Size-normalized presence by cell type (bars)
3. Presence on reference UMAP (blue→red)

Per-cell presence scores are **computed here** (not in 00):
  estimate_presence_score on the saved wkNN from
  00_reference_umaps/full_projection/{ref}/tables/wknn.npz
and written to
  01_presence_coverage/tables/{ref}_presence_on_reference.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mapping_data import load_adult, load_fetal  # noqa: E402
from wknn_mapping import estimate_presence_score  # noqa: E402
from lib_config import (  # noqa: E402
    ADULT_IMMUNE_TYPES,
    CONF_THRESHOLD,
    DIRS,
    FETAL_IMMUNE_TYPES,
    PATHS,
    PRESENCE_MIN_CELLS,
    PRESENCE_MIN_FRAC,
    PROTOCOLS_3,
    REF_MIN_CELLS,
    ROOT,
    apply_fig_style,
    build_adult_fine_to_broad,
    build_fetal_author_to_broad,
    ensure_dirs,
    merge_protocol,
    read_h5_obs_col,
)

apply_fig_style()

MAP_ROOT = ROOT / "00_reference_umaps" / "full_projection"
N0_SIZE = 100
PRESENT_MIN = 0.25
ABSENT_MAX = 0.05
RWR_ALPHA = 0.1
RWR_ROUNDS = 100

DISPLAY = {
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
    "MD": "Macula densa",
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
    "immune cells": "Immune cells",
}


def pretty(name: str) -> str:
    s = str(name)
    return DISPLAY.get(s, s.replace("_", " "))


# ── 1) composition heatmap ──────────────────────────────────────────────────


def classify_detection(n_conf, frac_org, median_conf, n_ref) -> str:
    if n_ref < REF_MIN_CELLS:
        return "Ref_rare"
    if (
        n_conf >= PRESENCE_MIN_CELLS
        and frac_org >= PRESENCE_MIN_FRAC
        and np.isfinite(median_conf)
        and median_conf >= CONF_THRESHOLD
    ):
        return "Present"
    if n_conf >= 10:
        return "Trace"
    return "Absent"


def compute_composition(ref_name, ref_broad, org_broad, org_conf, org_proto) -> pd.DataFrame:
    ref_ok = ~np.isin(ref_broad, ["NA", "Unassigned"])
    ref_vc = pd.Series(ref_broad[ref_ok]).value_counts()
    n_ref_tot = int(ref_ok.sum())
    rows = []
    for group in ["All organoids"] + PROTOCOLS_3:
        mask = np.ones(len(org_conf), dtype=bool) if group == "All organoids" else (org_proto == group)
        n_org = int(mask.sum())
        if n_org == 0:
            continue
        lab, conf = org_broad[mask], org_conf[mask]
        conf_ok = conf >= CONF_THRESHOLD
        for ct, n_ref in ref_vc.items():
            m = (lab == ct) & conf_ok
            n_c = int(m.sum())
            med = float(np.median(conf[m])) if n_c else np.nan
            rows.append(
                {
                    "reference": ref_name,
                    "protocol_group": group,
                    "cell_type": ct,
                    "n_ref": int(n_ref),
                    "pct_ref": 100.0 * n_ref / n_ref_tot,
                    "n_org_total_group": n_org,
                    "n_org_mapped_conf": n_c,
                    "pct_org_representation": 100.0 * n_c / n_org,
                    "median_conf_confident": med,
                    "detection_status": classify_detection(n_c, n_c / n_org, med, int(n_ref)),
                }
            )
    return pd.DataFrame(rows)


def collapse_immune_composition(df: pd.DataFrame, immune: set[str]) -> pd.DataFrame:
    if not immune:
        return df
    keep = df[~df["cell_type"].isin(immune)].copy()
    im = df[df["cell_type"].isin(immune)]
    if im.empty:
        return df
    agg = im.groupby(["reference", "protocol_group"], as_index=False).agg(
        n_ref=("n_ref", "sum"),
        pct_ref=("pct_ref", "sum"),
        n_org_total_group=("n_org_total_group", "first"),
        n_org_mapped_conf=("n_org_mapped_conf", "sum"),
        pct_org_representation=("pct_org_representation", "sum"),
    )
    agg["cell_type"] = "immune cells"
    agg["median_conf_confident"] = np.nan
    agg["detection_status"] = np.where(
        agg["n_org_mapped_conf"] >= PRESENCE_MIN_CELLS, "Present", "Trace"
    )
    return pd.concat([keep, agg], ignore_index=True)


def plot_composition_heatmap(df: pd.DataFrame, ref_name: str, out_png: Path) -> None:
    cols = ["All organoids"] + [p for p in PROTOCOLS_3 if p in set(df["protocol_group"])]
    all_df = df[df["protocol_group"] == "All organoids"].copy()
    order = all_df.sort_values("pct_org_representation", ascending=False)["cell_type"].tolist()
    piv_pct = (
        df.pivot(index="cell_type", columns="protocol_group", values="pct_org_representation")
        .reindex(index=order, columns=cols).fillna(0.0)
    )
    piv_n = (
        df.pivot(index="cell_type", columns="protocol_group", values="n_org_mapped_conf")
        .reindex(index=order, columns=cols).fillna(0.0)
    )
    unmapped_pct, unmapped_n = {}, {}
    for c in cols:
        sub = df[df["protocol_group"] == c]
        n_tot = int(sub["n_org_total_group"].iloc[0])
        n_map = int(sub["n_org_mapped_conf"].sum())
        n_u = max(0, n_tot - n_map)
        unmapped_n[c] = float(n_u)
        unmapped_pct[c] = 100.0 * n_u / max(n_tot, 1)
    piv_pct.loc["Unmapped"] = pd.Series(unmapped_pct)
    piv_n.loc["Unmapped"] = pd.Series(unmapped_n)
    row_order = order + ["Unmapped"]
    piv_pct, piv_n = piv_pct.reindex(row_order), piv_n.reindex(row_order)

    annot = pd.DataFrame(index=row_order, columns=cols, dtype=object)
    for r in row_order:
        for c in cols:
            pct, n = float(piv_pct.loc[r, c]), float(piv_n.loc[r, c])
            n_str = f"{n/1000:.1f}k" if n >= 1000 else f"{int(n)}"
            annot.loc[r, c] = "—" if n == 0 else f"{pct:.1f}\n({n_str})"

    meta = all_df.set_index("cell_type").reindex(order)
    ref_pct, n_ref = meta["pct_ref"].fillna(0.0), meta["n_ref"].fillna(0).astype(int)
    n_rows = len(row_order)
    fig = plt.figure(figsize=(15.0, max(7.0, 0.38 * n_rows + 1.8)))
    gs = fig.add_gridspec(
        1, 3, width_ratios=[4.0, 1.45, 0.12], wspace=0.45,
        left=0.20, right=0.92, top=0.90, bottom=0.14,
    )
    ax, ax_ref, cax = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1], sharey=None), fig.add_subplot(gs[0, 2])
    ax_ref.sharey(ax)
    vmax = max(10.0, float(np.nanpercentile(piv_pct.loc[order].to_numpy(), 90)))
    sns.heatmap(
        piv_pct, ax=ax, cmap="Blues", annot=annot, fmt="",
        annot_kws={"size": 10, "linespacing": 0.9}, vmin=0.0, vmax=vmax,
        cbar_ax=cax, cbar_kws={"label": ""}, linewidths=0.4, linecolor="white",
    )
    cax.yaxis.set_label_position("right")
    cax.yaxis.tick_right()
    cax.set_ylabel("% organoid", fontsize=12, labelpad=8)
    ax.axhline(len(order), color="#333333", linewidth=1.0)
    ax.set_title(f"Organoid composition on {ref_name}", fontsize=16, pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=25, labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    y = np.arange(len(order)) + 0.5
    ax_ref.barh(y, ref_pct.reindex(order).to_numpy(), height=0.7, color="#9E9E9E")
    for yi, ct in zip(y, order):
        ax_ref.text(
            float(ref_pct.get(ct, 0)) + 0.3, yi, f"n={int(n_ref.get(ct, 0)):,}",
            va="center", ha="left", fontsize=10, color="#333333",
        )
    ax_ref.set_xlim(0, max(12.0, float(ref_pct.max()) * 1.45))
    ax_ref.set_xlabel(f"{ref_name} ref %", fontsize=12)
    ax_ref.set_title("Primary composition", fontsize=13)
    ax_ref.set_ylim(n_rows, 0)
    ax_ref.tick_params(axis="y", labelleft=False)
    ax_ref.spines[["top", "right"]].set_visible(False)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white", pad_inches=0.25)
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


# ── 2) presence by cell type ─────────────────────────────────────────────────


def collapse_immune_presence(df: pd.DataFrame, immune: set[str]) -> pd.DataFrame:
    if not immune:
        return df
    keep = df[~df["cell_type"].isin(immune)].copy()
    im = df[df["cell_type"].isin(immune)]
    if im.empty:
        return df
    n = im["n_ref"].to_numpy(float)
    n_tot = int(im["n_ref"].sum())
    mean_p = float(np.average(im["mean_presence"], weights=n)) if n.sum() else 0.0
    row = {
        "cell_type": "immune cells",
        "n_ref": n_tot,
        "pct_ref": float(im["pct_ref"].sum()),
        "mean_presence": mean_p,
        "median_presence": float(np.average(im["median_presence"], weights=n)) if n.sum() else 0.0,
        "sum_presence": float(im["sum_presence"].sum()),
        "presence_norm": mean_p * n_tot / (n_tot + N0_SIZE),
    }
    return pd.concat([keep, pd.DataFrame([row])], ignore_index=True)


def aggregate_presence(presence: np.ndarray, labels: np.ndarray, immune: set[str]) -> pd.DataFrame:
    ok = ~np.isin(labels, ["NA", "Unassigned", "unassigned"])
    presence, labels = presence[ok], labels[ok]
    n_tot = len(labels)
    rows = []
    for ct in sorted(set(labels.tolist())):
        vals = presence[labels == ct]
        n = len(vals)
        mean_p = float(vals.mean())
        rows.append(
            {
                "cell_type": ct,
                "n_ref": n,
                "pct_ref": 100.0 * n / max(n_tot, 1),
                "mean_presence": mean_p,
                "median_presence": float(np.median(vals)),
                "sum_presence": float(vals.sum()),
                "presence_norm": mean_p * n / (n + N0_SIZE),
            }
        )
    df = collapse_immune_presence(pd.DataFrame(rows), immune)
    df["status"] = np.where(
        df["presence_norm"] >= PRESENT_MIN, "Present",
        np.where(df["presence_norm"] <= ABSENT_MAX, "Absent", "Low"),
    )
    df["display_name"] = df["cell_type"].map(pretty)
    return df.sort_values("presence_norm", ascending=False).reset_index(drop=True)


def plot_presence_bars(df: pd.DataFrame, ref_name: str, out_png: Path) -> None:
    show = df.sort_values("presence_norm", ascending=True)
    colors = {"Present": "#1565C0", "Low": "#90CAF9", "Absent": "#BDBDBD"}
    fig, ax = plt.subplots(figsize=(11.0, max(6.0, 0.40 * len(show) + 2.2)))
    y = np.arange(len(show))
    vals = show["presence_norm"].to_numpy(float)
    ax.barh(y, vals, color=[colors[s] for s in show["status"]], height=0.75)
    ax.axvline(PRESENT_MIN, color="#1565C0", ls="--", lw=1.0, alpha=0.85)
    ax.axvline(ABSENT_MAX, color="#757575", ls=":", lw=1.0, alpha=0.85)
    xmax = min(1.15, float(vals.max()) * 1.45 + 0.08)
    for yi, (_, r) in enumerate(show.iterrows()):
        ax.text(
            float(r["presence_norm"]) + 0.012, yi,
            f"{r['presence_norm']:.2f}  (mean={r['mean_presence']:.2f}, n={int(r['n_ref']):,})  {r['status']}",
            va="center", ha="left", fontsize=10, color="#333333",
        )
    ax.set_yticks(y)
    ax.set_yticklabels(show["display_name"].tolist(), fontsize=12)
    ax.set_xlabel("Size-normalized presence score", fontsize=14)
    ax.set_xlim(0, xmax)
    ax.set_title(
        f"Organoid presence on {ref_name} primary cell types\n"
        "Score normalized for cell-type size",
        fontsize=14, pad=10,
    )
    ax.legend(
        handles=[
            Patch(facecolor=colors["Present"], label=f"Present (≥{PRESENT_MIN:g})"),
            Patch(facecolor=colors["Low"], label="Low"),
            Patch(facecolor=colors["Absent"], label=f"Absent (≤{ABSENT_MAX:g})"),
        ],
        loc="lower right", fontsize=11, frameon=False,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white", pad_inches=0.25)
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


# ── 3) presence scores + UMAP ────────────────────────────────────────────────


def compute_presence_scores(ref_name: str, tab_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run estimate_presence_score; save CSV under 01 tables."""
    wknn_path = MAP_ROOT / ref_name / "tables" / "wknn.npz"
    if not wknn_path.exists():
        raise FileNotFoundError(
            f"Missing {wknn_path}\n"
            "Run first: python scripts/00_reference_umaps/compute_mapping.py"
        )
    print(f"  loading wkNN {wknn_path.name}…", flush=True)
    wknn = sp.load_npz(wknn_path).tocsr()

    if ref_name == "adult":
        ref_lat, ref_umap, labels, org_lat, org_l2, _key = load_adult()
    else:
        ref_lat, ref_umap, labels, org_lat, org_l2, _key = load_fetal()
    labels = np.asarray(labels, dtype=object)

    if wknn.shape[0] != org_lat.shape[0]:
        wknn = wknn.T.tocsr()
    if wknn.shape != (org_lat.shape[0], ref_lat.shape[0]):
        raise ValueError(
            f"wkNN shape {wknn.shape} != query×ref "
            f"({org_lat.shape[0]}, {ref_lat.shape[0]})"
        )

    ref_ad = ad.AnnData(X=sp.csr_matrix((ref_lat.shape[0], 1)))
    ref_ad.obs_names = [f"ref_{i}" for i in range(ref_lat.shape[0])]
    ref_ad.obsm["X_latent"] = np.asarray(ref_lat, dtype=np.float32)

    qry_ad = ad.AnnData(X=sp.csr_matrix((org_lat.shape[0], 1)))
    qry_ad.obs["Level_2"] = pd.Categorical(np.asarray(org_l2, dtype=str))
    qry_ad.obs_names = [f"qry_{i}" for i in range(org_lat.shape[0])]
    qry_ad.obsm["X_latent"] = np.asarray(org_lat, dtype=np.float32)

    print(
        f"  estimate_presence_score (RWR α={RWR_ALPHA}, rounds={RWR_ROUNDS})…",
        flush=True,
    )
    pres = estimate_presence_score(
        ref_ad,
        qry_ad,
        wknn=wknn,
        use_rep_ref_wknn="X_latent",
        use_rep_query_wknn="X_latent",
        do_random_walk=True,
        alpha_random_walk=RWR_ALPHA,
        num_rounds_random_walk=RWR_ROUNDS,
        log=True,
    )
    presence = np.asarray(pres["max"], dtype=np.float32)
    out_csv = tab_dir / f"{ref_name}_presence_on_reference.csv"
    pd.DataFrame({"presence_max": presence}).to_csv(out_csv, index=False)
    print(f"  Wrote {out_csv}", flush=True)
    return presence, np.asarray(ref_umap, dtype=np.float32), labels


def plot_presence_umap(ref_umap, presence, out_png, ref: str) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    presence = np.asarray(presence, dtype=np.float64)
    order = np.argsort(presence)
    sc = ax.scatter(
        ref_umap[order, 0], ref_umap[order, 1],
        s=1.4, c=presence[order], cmap="RdBu_r", vmin=0, vmax=1,
        alpha=0.85, linewidths=0, rasterized=True,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Presence score", fontsize=13)
    ax.set_title(f"Organoid presence on {ref} reference", fontsize=16, pad=8)
    ax.set_xlabel("UMAP1", fontsize=14)
    ax.set_ylabel("UMAP2", fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    pad = 0.8
    ax.set_xlim(float(ref_umap[:, 0].min()) - pad, float(ref_umap[:, 0].max()) + pad)
    ax.set_ylim(float(ref_umap[:, 1].min()) - pad, float(ref_umap[:, 1].max()) + pad)
    ax.set_aspect("equal", adjustable="box")
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


# ── driver ───────────────────────────────────────────────────────────────────


def run_ref(ref_name: str) -> None:
    print(f"=== Presence {ref_name} (3 figures) ===", flush=True)
    fig_dir, tab_dir = DIRS["01_fig"], DIRS["01_tab"]

    # (1) composition
    if ref_name == "adult":
        with h5py.File(PATHS["adult_meta_broad"], "r") as fm:
            ref_broad = read_h5_obs_col(fm, "broad_annotation")
        fine_to_broad = build_adult_fine_to_broad()
        with h5py.File(PATHS["org_on_adult"], "r") as fo:
            wknn_lab = read_h5_obs_col(fo, "wkNN_label")
            org_conf = np.asarray(fo["obs"]["wkNN_confidence"][:], dtype=np.float32)
            org_proto = read_h5_obs_col(fo, "diff_protocol")
        org_broad = np.array([fine_to_broad.get(str(x), "NA") for x in wknn_lab], dtype=object)
        immune = ADULT_IMMUNE_TYPES
    else:
        with h5py.File(PATHS["fetal_umap"], "r") as fr:
            ref_broad = read_h5_obs_col(fr, "harmonized_celltype")
        auth_to_broad = build_fetal_author_to_broad()
        with h5py.File(PATHS["org_on_fetal"], "r") as fo:
            wknn_lab = read_h5_obs_col(fo, "wkNN_label")
            org_conf = np.asarray(fo["obs"]["wkNN_confidence"][:], dtype=np.float32)
            org_proto = read_h5_obs_col(fo, "diff_protocol")
        org_broad = np.array([auth_to_broad.get(str(x), "NA") for x in wknn_lab], dtype=object)
        immune = FETAL_IMMUNE_TYPES

    org_proto = np.array([merge_protocol(p) for p in org_proto], dtype=object)
    by_proto = compute_composition(ref_name, ref_broad, org_broad, org_conf, org_proto)
    by_proto.to_csv(tab_dir / f"{ref_name}_presence_by_protocol.csv", index=False)
    collapsed = collapse_immune_composition(by_proto, immune)
    collapsed.to_csv(tab_dir / f"{ref_name}_presence_by_protocol_grouped.csv", index=False)
    plot_composition_heatmap(
        collapsed, ref_name, fig_dir / f"{ref_name}_presence_by_protocol_heatmap.png"
    )

    # (2)+(3) compute presence here, then bars + UMAP
    presence, ref_umap, labels = compute_presence_scores(ref_name, tab_dir)

    pres_df = aggregate_presence(presence, labels, immune)
    pres_df.to_csv(tab_dir / f"{ref_name}_presence_by_celltype.csv", index=False)
    plot_presence_bars(pres_df, ref_name, fig_dir / f"{ref_name}_presence_by_celltype.png")
    plot_presence_umap(
        ref_umap, presence, fig_dir / f"{ref_name}_presence_umap.png", ref_name
    )


def main() -> None:
    ensure_dirs()
    DIRS["01_fig"].mkdir(parents=True, exist_ok=True)
    run_ref("adult")
    run_ref("fetal")
    print("Presence complete (3 figures × adult/fetal).", flush=True)


if __name__ == "__main__":
    main()
