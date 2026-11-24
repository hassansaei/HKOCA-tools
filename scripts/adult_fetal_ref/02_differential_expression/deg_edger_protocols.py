#!/usr/bin/env python
"""Stage 2: Pseudobulk Differential Expression (edgeR QLF) per cell-type x protocol
vs matched Fetal and Adult reference cells.

- Uses ONLY the 3 merged protocols: Morizane (+ modified), Takasato (+ modified), Freedman
- Enforces >= 3 replicates per group (sample_id for org/adult, specimen for fetal)
- Includes age as a within-group centered covariate (+ QC covariates)
- Tests the 9 target cell types:
    Distal_tubule, Nephron_progenitor_cells, PEC-like, Podocytes,
    Proliferative_nephron, Proliferative_stroma, Proximal_tubule, S1, S2
- Quantifies Shared Up, Shared Down, Inconsistent, and Protocol-specific DEGs
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib_config import (  # noqa: E402
    DEG_FDR,
    DEG_LFC,
    DIRS,
    EDGER_R,
    MIN_CELLS_PB,
    MIN_COUNTS_PB,
    MIN_REPLICATES,
    PAIRS_L2,
    PATHS,
    PROTOCOL_COLORS,
    PROTOCOLS_3,
    RSCRIPT,
    TARGET_CELLTYPES_9,
    apply_fig_style,
    ensure_dirs,
    merge_protocol,
    parse_age_numeric,
    read_h5_obs_col,
)

apply_fig_style()


def read_var_genes(f: h5py.File) -> np.ndarray:
    k = "_index" if "_index" in f["var"] else "gene_symbol"
    arr = f["var"][k][:]
    return np.array([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr])


def aggregate_pseudobulks_csr(
    h5_path: Path,
    group_keys: np.ndarray,
    keep_mask: np.ndarray,
    age_vals: np.ndarray,
    chunk_size: int = 25000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stream CSR matrix X in chunks and sum raw counts per group_key.
    Also computes n_cells, med_genes (detected > 0), and median age_val per group.
    """
    keep_idx = np.flatnonzero(keep_mask)
    if len(keep_idx) == 0:
        return pd.DataFrame(), pd.DataFrame()

    unique_keys, inv = np.unique(group_keys[keep_idx], return_inverse=True)
    n_groups = len(unique_keys)
    # Map cell index -> group id (-1 if not kept)
    with h5py.File(h5_path, "r") as f:
        genes = read_var_genes(f)
        n_genes = len(genes)
        g = f["X"]
        indptr = g["indptr"][:]
        n_obs = len(indptr) - 1

        cell_to_grp = np.full(n_obs, -1, dtype=np.int32)
        cell_to_grp[keep_idx] = inv.astype(np.int32)

        counts_mat = np.zeros((n_groups, n_genes), dtype=np.float64)
        n_genes_per_cell = np.zeros(len(keep_idx), dtype=np.int32)

        # Track position in keep_idx for detected genes per cell
        pos_in_keep = np.full(n_obs, -1, dtype=np.int32)
        pos_in_keep[keep_idx] = np.arange(len(keep_idx), dtype=np.int32)

        indices_ds = g["indices"]
        data_ds = g["data"]

        for start in range(0, n_obs, chunk_size):
            end = min(n_obs, start + chunk_size)
            grp_chunk = cell_to_grp[start:end]
            local_kept = np.flatnonzero(grp_chunk >= 0)
            if len(local_kept) == 0:
                continue

            a = int(indptr[start])
            b = int(indptr[end])
            ind_c = indices_ds[a:b]
            dat_c = data_ds[a:b]
            ip_c = indptr[start : end + 1] - a

            for loc in local_kept:
                r0 = int(ip_c[loc])
                r1 = int(ip_c[loc + 1])
                if r1 > r0:
                    cols = ind_c[r0:r1]
                    vals = dat_c[r0:r1]
                    np.add.at(counts_mat[grp_chunk[loc]], cols, vals)
                    n_genes_per_cell[pos_in_keep[start + loc]] = int((vals > 0).sum())

    # Build metadata per group
    meta_rows = []
    keep_ages = age_vals[keep_idx]
    for gi, key in enumerate(unique_keys):
        m = inv == gi
        n_c = int(m.sum())
        tot = float(counts_mat[gi].sum())
        med_g = float(np.median(n_genes_per_cell[m])) if n_c > 0 else 0.0
        a_sub = keep_ages[m]
        a_fin = a_sub[np.isfinite(a_sub)]
        med_age = float(np.median(a_fin)) if len(a_fin) > 0 else np.nan
        meta_rows.append(
            {
                "key": str(key),
                "n_cells": n_c,
                "total_counts": tot,
                "med_genes": med_g,
                "age_val": med_age,
            }
        )
    meta_df = pd.DataFrame(meta_rows)
    counts_df = pd.DataFrame(
        np.rint(counts_mat.T).astype(np.int64),
        index=genes,
        columns=[str(k) for k in unique_keys],
    )
    return counts_df, meta_df


def build_all_pseudobulks() -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    print("=== Building pseudobulks (Organoid 3 protocols, Fetal, Adult) ===", flush=True)

    # 1. Organoid: key = Level_2 || protocol || sample_id
    with h5py.File(PATHS["organoid_counts"], "r") as f:
        l2 = read_h5_obs_col(f, "Level_2")
        proto_raw = read_h5_obs_col(f, "diff_protocol")
        sid = read_h5_obs_col(f, "sample_id")
        age_raw = read_h5_obs_col(f, "Age")
    proto = np.array([merge_protocol(p) for p in proto_raw], dtype=object)
    org_age = np.array([parse_age_numeric(a) for a in age_raw], dtype=float)
    org_keep = np.isin(l2, TARGET_CELLTYPES_9) & np.isin(proto, PROTOCOLS_3)
    org_keys = np.array([f"{a}||{b}||{c}" for a, b, c in zip(l2, proto, sid)], dtype=object)
    org_c, org_m = aggregate_pseudobulks_csr(PATHS["organoid_counts"], org_keys, org_keep, org_age)
    parts = org_m["key"].str.split(r"\|\|", expand=True)
    org_m["cell_type"] = parts[0]
    org_m["protocol"] = parts[1]
    org_m["sample_id"] = parts[2]
    print(f"  Organoid pseudobulks: {len(org_m)} raw", flush=True)

    # 2. Fetal: aggregate per (target_celltype || specimen) since fetal pairs can span multiple harmonized types
    with h5py.File(PATHS["fetal_counts"], "r") as f:
        f_harm = read_h5_obs_col(f, "harmonized_celltype")
        f_spec = read_h5_obs_col(f, "specimen")
        f_wk = np.asarray(f["obs"]["week"][:], dtype=float) * 7.0  # days
    # Aggregate by harmonized_celltype || specimen first (exact, no double reading!)
    f_keys = np.array([f"{a}||{b}" for a, b in zip(f_harm, f_spec)], dtype=object)
    f_c, f_m = aggregate_pseudobulks_csr(
        PATHS["fetal_counts"], f_keys, np.ones(len(f_harm), dtype=bool), f_wk
    )
    f_parts = f_m["key"].str.split(r"\|\|", expand=True)
    f_m["ref_label"] = f_parts[0]
    f_m["sample_id"] = f_parts[1]
    print(f"  Fetal pseudobulks: {len(f_m)} raw", flush=True)

    # 3. Adult: aggregate by fine_annotation || sample_id (only needed fine labels!)
    needed_adult = sorted({lab for _, (_, alabs) in PAIRS_L2.items() for lab in alabs})
    with h5py.File(PATHS["adult_counts"], "r") as f, h5py.File(
        PATHS["adult_meta_broad"], "r"
    ) as fm:
        a_fine = read_h5_obs_col(f, "fine_annotation")
        a_sid = read_h5_obs_col(f, "sample_id")
        a_age_raw = read_h5_obs_col(fm, "age_range")
    a_age = np.array([parse_age_numeric(a) for a in a_age_raw], dtype=float)
    a_keep = np.isin(a_fine, needed_adult)
    a_keys = np.array([f"{a}||{b}" for a, b in zip(a_fine, a_sid)], dtype=object)
    a_c, a_m = aggregate_pseudobulks_csr(PATHS["adult_counts"], a_keys, a_keep, a_age)
    a_parts = a_m["key"].str.split(r"\|\|", expand=True)
    a_m["ref_label"] = a_parts[0]
    a_m["sample_id"] = a_parts[1]
    print(f"  Adult pseudobulks: {len(a_m)} raw", flush=True)

    return {"org": (org_c, org_m), "fetal": (f_c, f_m), "adult": (a_c, a_m)}


def collapse_primary_for_celltype(
    ref_c: pd.DataFrame,
    ref_m: pd.DataFrame,
    ref_labels: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum reference pseudobulks across matching ref_labels per sample_id."""
    if not ref_labels:
        return pd.DataFrame(), pd.DataFrame()
    sub_m = ref_m[ref_m["ref_label"].isin(ref_labels)].copy()
    if sub_m.empty:
        return pd.DataFrame(), pd.DataFrame()

    cols_dict = {}
    meta_rows = []
    for sid, g in sub_m.groupby("sample_id"):
        keys = g["key"].tolist()
        summed = ref_c[keys].sum(axis=1)
        n_c = int(g["n_cells"].sum())
        tot = float(summed.sum())
        if n_c < MIN_CELLS_PB or tot < MIN_COUNTS_PB:
            continue
        w = g["n_cells"].to_numpy(float)
        med_g = float(np.average(g["med_genes"].to_numpy(float), weights=w))
        ages = g["age_val"].to_numpy(float)
        fin = np.isfinite(ages)
        age_v = float(np.average(ages[fin], weights=w[fin])) if np.any(fin) else np.nan
        col_name = f"ref__{sid}"
        cols_dict[col_name] = summed
        meta_rows.append(
            {
                "sample": col_name,
                "group": "Control",
                "n_cells": n_c,
                "total_counts": tot,
                "med_genes": med_g,
                "age_val": age_v,
            }
        )
    if not cols_dict:
        return pd.DataFrame(), pd.DataFrame()
    return pd.DataFrame(cols_dict, index=ref_c.index), pd.DataFrame(meta_rows)


def run_edger_contrast(
    ref_c: pd.DataFrame,
    ref_m: pd.DataFrame,
    org_c: pd.DataFrame,
    org_m: pd.DataFrame,
    outdir: Path,
) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    common_genes = ref_c.index.intersection(org_c.index)
    comb_c = pd.concat([ref_c.loc[common_genes], org_c.loc[common_genes]], axis=1)
    comb_m = pd.concat([ref_m, org_m], ignore_index=True)

    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        c_path = td_p / "counts.tsv"
        m_path = td_p / "meta.tsv"
        comb_c.to_csv(c_path, sep="\t")
        comb_m.to_csv(m_path, sep="\t", index=False)

        proc = subprocess.run(
            [str(RSCRIPT), str(EDGER_R), str(c_path), str(m_path), str(outdir)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"edgeR failed for {outdir}:\n{proc.stderr}\n{proc.stdout}")

    return json.loads((outdir / "design.json").read_text())


def run_deg_for_reference(
    ref_name: str,
    pbs: dict[str, tuple[pd.DataFrame, pd.DataFrame]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print(f"\n=== Stage 2: edgeR DEG vs {ref_name} ===", flush=True)
    org_c, org_m = pbs["org"]
    ref_c_raw, ref_m_raw = pbs[ref_name]
    out_base = DIRS["02_tab"] / f"{ref_name}_deg_per_celltype_protocol"
    out_base.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    consensus_summary_rows = []
    consensus_gene_rows = []

    for ct in TARGET_CELLTYPES_9:
        flabs, alabs = PAIRS_L2[ct]
        ref_labs = alabs if ref_name == "adult" else flabs
        ref_c, ref_m = collapse_primary_for_celltype(ref_c_raw, ref_m_raw, ref_labs)
        n_ref = len(ref_m)

        proto_tables: dict[str, pd.DataFrame] = {}

        for proto in PROTOCOLS_3:
            sub_om = org_m[
                (org_m["cell_type"] == ct)
                & (org_m["protocol"] == proto)
                & (org_m["n_cells"] >= MIN_CELLS_PB)
                & (org_m["total_counts"] >= MIN_COUNTS_PB)
            ].copy()
            n_org = len(sub_om)

            if n_ref < MIN_REPLICATES or n_org < MIN_REPLICATES:
                reason = (
                    "no_reference_counterpart"
                    if len(ref_labs) == 0
                    else f"insufficient_replicates(ref={n_ref},org={n_org})"
                )
                summary_rows.append(
                    {
                        "reference": ref_name,
                        "cell_type": ct,
                        "protocol": proto,
                        "status": "skipped",
                        "reason": reason,
                        "n_ref_replicates": n_ref,
                        "n_org_replicates": n_org,
                        "n_genes_tested": 0,
                        "n_sig": 0,
                        "n_sig_up": 0,
                        "n_sig_down": 0,
                    }
                )
                continue

            om_clean = pd.DataFrame(
                {
                    "sample": [f"org__{s}" for s in sub_om["sample_id"]],
                    "group": "Query",
                    "n_cells": sub_om["n_cells"].to_numpy(),
                    "total_counts": sub_om["total_counts"].to_numpy(),
                    "med_genes": sub_om["med_genes"].to_numpy(),
                    "age_val": sub_om["age_val"].to_numpy(),
                }
            )
            oc_clean = org_c[sub_om["key"].tolist()].copy()
            oc_clean.columns = om_clean["sample"].tolist()

            outdir = out_base / f"{ct}__{proto}"
            info = run_edger_contrast(ref_c, ref_m, oc_clean, om_clean, outdir)
            deg_df = pd.read_csv(outdir / "deg_full.tsv", sep="\t")
            proto_tables[proto] = deg_df

            summary_rows.append(
                {
                    "reference": ref_name,
                    "cell_type": ct,
                    "protocol": proto,
                    "status": "tested",
                    "reason": "ok",
                    "n_ref_replicates": info["n_control"],
                    "n_org_replicates": info["n_query"],
                    "n_genes_tested": info["n_genes_tested"],
                    "n_sig": info["n_sig"],
                    "n_sig_up": info["n_sig_up"],
                    "n_sig_down": info["n_sig_down"],
                    "formula": info["formula"],
                }
            )
            print(
                f"  {ref_name} | {ct:26s} | {proto:8s} | ref={info['n_control']:3d}, org={info['n_query']:2d} "
                f"| sig={info['n_sig']:4d} (up={info['n_sig_up']:4d}, down={info['n_sig_down']:4d})",
                flush=True,
            )

        # Protocol consensus across the 3 protocols
        if len(proto_tables) == 3:
            all_genes = sorted(set.union(*[set(df["gene"]) for df in proto_tables.values()]))
            # Build lookup dicts
            dir_map = {}
            lfc_map = {}
            fdr_map = {}
            for p, df in proto_tables.items():
                sig_m = (df["FDR"] < DEG_FDR) & (df["logFC"].abs() > DEG_LFC)
                d = np.where(~sig_m, 0, np.where(df["logFC"] > 0, 1, -1))
                dir_map[p] = dict(zip(df["gene"], d))
                lfc_map[p] = dict(zip(df["gene"], df["logFC"]))
                fdr_map[p] = dict(zip(df["gene"], df["FDR"]))

            counts_cat = {
                "Shared_Up_3of3": 0,
                "Shared_Down_3of3": 0,
                "Consensus_Up_2of3": 0,
                "Consensus_Down_2of3": 0,
                "Inconsistent_direction": 0,
                "Protocol_specific_1of3": 0,
            }
            for g in all_genes:
                ds = [dir_map[p].get(g, 0) for p in PROTOCOLS_3]
                n_up = sum(x == 1 for x in ds)
                n_dn = sum(x == -1 for x in ds)
                if n_up == 0 and n_dn == 0:
                    continue
                if n_up > 0 and n_dn > 0:
                    cat = "Inconsistent_direction"
                elif n_up == 3:
                    cat = "Shared_Up_3of3"
                elif n_dn == 3:
                    cat = "Shared_Down_3of3"
                elif n_up == 2:
                    cat = "Consensus_Up_2of3"
                elif n_dn == 2:
                    cat = "Consensus_Down_2of3"
                else:
                    cat = "Protocol_specific_1of3"
                counts_cat[cat] += 1

                lfcs = [lfc_map[p].get(g, np.nan) for p in PROTOCOLS_3]
                fdrs = [fdr_map[p].get(g, np.nan) for p in PROTOCOLS_3]
                consensus_gene_rows.append(
                    {
                        "reference": ref_name,
                        "cell_type": ct,
                        "gene": g,
                        "category": cat,
                        "n_protocols_up": n_up,
                        "n_protocols_down": n_dn,
                        "mean_logFC": float(np.nanmean(lfcs)),
                        "Morizane_logFC": lfcs[0],
                        "Takasato_logFC": lfcs[1],
                        "Freedman_logFC": lfcs[2],
                        "Morizane_FDR": fdrs[0],
                        "Takasato_FDR": fdrs[1],
                        "Freedman_FDR": fdrs[2],
                    }
                )

            total_any = sum(counts_cat.values())
            consensus_summary_rows.append(
                {
                    "reference": ref_name,
                    "cell_type": ct,
                    "n_genes_union_tested": len(all_genes),
                    "n_sig_any_protocol": total_any,
                    **counts_cat,
                    "n_shared_3of3_total": counts_cat["Shared_Up_3of3"] + counts_cat["Shared_Down_3of3"],
                    "pct_shared_3of3_of_any": 100.0
                    * (counts_cat["Shared_Up_3of3"] + counts_cat["Shared_Down_3of3"])
                    / max(total_any, 1),
                    "pct_inconsistent_of_any": 100.0 * counts_cat["Inconsistent_direction"] / max(total_any, 1),
                }
            )

    sum_df = pd.DataFrame(summary_rows)
    cons_sum_df = pd.DataFrame(consensus_summary_rows)
    cons_gene_df = pd.DataFrame(consensus_gene_rows)

    sum_df.to_csv(DIRS["02_tab"] / f"{ref_name}_deg_summary.csv", index=False)
    cons_sum_df.to_csv(DIRS["02_tab"] / f"{ref_name}_protocol_consensus_summary.csv", index=False)
    cons_gene_df.to_csv(DIRS["02_tab"] / f"{ref_name}_protocol_consensus_genes.csv", index=False)

    plot_deg_counts_by_protocol(sum_df, ref_name, DIRS["02_fig"] / f"{ref_name}_deg_counts_by_protocol.png")
    if not cons_sum_df.empty:
        plot_protocol_consensus(cons_sum_df, ref_name, DIRS["02_fig"] / f"{ref_name}_deg_protocol_consensus.png")

    return sum_df, cons_sum_df, cons_gene_df


def plot_deg_counts_by_protocol(sum_df: pd.DataFrame, ref_name: str, out_png: Path) -> None:
    df = sum_df[sum_df["status"] == "tested"].copy()
    if df.empty:
        return
    cts = [c for c in TARGET_CELLTYPES_9 if c in df["cell_type"].values]
    x = np.arange(len(cts))
    width = 0.26

    fig, ax = plt.subplots(figsize=(max(11, 1.35 * len(cts) + 2.5), 7))
    for i, proto in enumerate(PROTOCOLS_3):
        sub = df[df["protocol"] == proto].set_index("cell_type").reindex(cts)
        up = sub["n_sig_up"].fillna(0).to_numpy()
        dn = -sub["n_sig_down"].fillna(0).to_numpy()
        pos = x + (i - 1) * width
        col = PROTOCOL_COLORS[proto]
        ax.bar(pos, up, width=width * 0.92, color=col, alpha=0.92, label=f"{proto} (Up in org)")
        ax.bar(pos, dn, width=width * 0.92, color=col, alpha=0.42, hatch="//", edgecolor=col, linewidth=0.5)

    ax.axhline(0, color="#222222", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(cts, rotation=30, ha="right", fontsize=12)
    ax.set_ylabel("Number of DEGs (+ Up / − Down in organoid)", fontsize=14)
    ax.set_title(
        f"edgeR pseudobulk DEGs vs {ref_name} reference (FDR < {DEG_FDR}, |log2FC| > {DEG_LFC}, age-adjusted)",
        fontsize=15,
    )
    ax.legend(loc="upper right", frameon=False, fontsize=12, ncol=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def plot_protocol_consensus(cons_df: pd.DataFrame, ref_name: str, out_png: Path) -> None:
    """Two-panel story figure:
    (A) Stacked horizontal bar of Shared Up (3/3), Shared Down (3/3), Consensus (2/3),
        Inconsistent, and Protocol-specific (1/3) gene counts
    (B) Proportion (%) of Shared (3/3) vs Inconsistent DEGs across protocols
    """
    df = cons_df.copy().sort_values("n_shared_3of3_total", ascending=True).reset_index(drop=True)
    y = np.arange(len(df))

    fig, axes = plt.subplots(1, 2, figsize=(18, max(6.2, 0.72 * len(df) + 2.0)), sharey=True)

    # Panel A: Absolute gene counts
    ax = axes[0]
    cats = [
        ("Shared_Up_3of3", "Shared Up (3/3)", "#D62728"),
        ("Shared_Down_3of3", "Shared Down (3/3)", "#1F77B4"),
        ("Consensus_Up_2of3", "Consensus Up (2/3)", "#FF9896"),
        ("Consensus_Down_2of3", "Consensus Down (2/3)", "#AEC7E8"),
        ("Inconsistent_direction", "Inconsistent (opp. sign)", "#9467BD"),
        ("Protocol_specific_1of3", "Protocol-specific (1/3)", "#C7C7C7"),
    ]
    left = np.zeros(len(df))
    for col, lab, color in cats:
        vals = df[col].to_numpy()
        ax.barh(y, vals, left=left, height=0.65, color=color, label=lab)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(df["cell_type"], fontsize=13)
    ax.set_xlabel("Number of significant DEGs (union across 3 protocols)", fontsize=14)
    ax.set_title(f"A. Protocol concordance of DEGs vs {ref_name}", fontsize=15)
    ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.25)

    # Panel B: Breakdown of Shared 3/3 (Up vs Down) vs Inconsistent
    ax = axes[1]
    up3 = df["Shared_Up_3of3"].to_numpy()
    dn3 = df["Shared_Down_3of3"].to_numpy()
    inc = df["Inconsistent_direction"].to_numpy()
    ax.barh(y + 0.22, up3, height=0.22, color="#D62728", label="Shared Up (3/3)")
    ax.barh(y, dn3, height=0.22, color="#1F77B4", label="Shared Down (3/3)")
    ax.barh(y - 0.22, inc, height=0.22, color="#9467BD", label="Inconsistent")
    for yi, u, d, ic in zip(y, up3, dn3, inc):
        ax.text(u + 15, yi + 0.22, str(int(u)), va="center", fontsize=11)
        ax.text(d + 15, yi, str(int(d)), va="center", fontsize=11)
        ax.text(ic + 15, yi - 0.22, str(int(ic)), va="center", fontsize=11)
    ax.set_xlabel("Number of genes", fontsize=14)
    ax.set_title(f"B. Shared Up / Down (3/3) vs Inconsistent ({ref_name})", fontsize=15)
    ax.set_xlim(0, max(50, float(max(up3.max(), dn3.max(), inc.max())) * 1.22))
    ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#DDDDDD", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def main() -> None:
    ensure_dirs()
    pbs = build_all_pseudobulks()

    # Save concise pseudobulk summary
    pb_rows = []
    for k, (_, m) in pbs.items():
        pb_rows.append(
            {
                "source": k,
                "n_raw_pseudobulks": len(m),
                "n_passing_qc": int(((m["n_cells"] >= MIN_CELLS_PB) & (m["total_counts"] >= MIN_COUNTS_PB)).sum()),
            }
        )
    pd.DataFrame(pb_rows).to_csv(DIRS["02_tab"] / "pseudobulk_summary.csv", index=False)

    run_deg_for_reference("adult", pbs)
    run_deg_for_reference("fetal", pbs)
    print("Stage 2 (edgeR DEG) complete.", flush=True)


if __name__ == "__main__":
    main()
