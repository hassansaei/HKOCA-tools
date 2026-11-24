#!/usr/bin/env python
"""Preranked GSEA on organoid vs primary DEGs (complement to ORA).

For each reference × cell type we build a gene ranking from the mean logFC
across the three protocols (genes seen in ≥2 protocols), then run gseapy
prerank on Hallmark + GO Biological Process + Reactome.

Positive NES → enrichment toward organoid-high genes.
Negative NES → enrichment toward primary-high genes.

Outputs under 02_differential_expression/gsea/.
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

warnings.filterwarnings("ignore", category=UserWarning)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib_config import (  # noqa: E402
    DIRS,
    PATHS,
    PROTOCOLS_3,
    apply_fig_style,
    ensure_dirs,
)

apply_fig_style()

GMT_FILES = {
    "Hallmark": "Enrichr.MSigDB_Hallmark_2020.gmt",
    "GO_BP": "Enrichr.GO_Biological_Process_2023.gmt",
    "Reactome": "Enrichr.Reactome_2022.gmt",
}
MIN_PROTOCOLS = 2
PERM = 500
THREADS = 8
FDR_CUT = 0.25  # GSEA often uses 0.25 for exploratory FDR
TOP_N_HEAT = 25
TOP_N_BARS = 20


def load_gmt_dict() -> dict[str, set[str]]:
    """Combine selected GMTs into one gene-set dict (name prefixed by library)."""
    out: dict[str, set[str]] = {}
    for lib, fname in GMT_FILES.items():
        p = PATHS["gmt_dir"] / fname
        if not p.exists():
            print(f"  WARN missing GMT {p}", flush=True)
            continue
        with open(p, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                genes = {g.strip().upper() for g in parts[2:] if g.strip()}
                if 15 <= len(genes) <= 500:
                    out[f"[{lib}] {parts[0].strip()}"] = genes
    return out


def build_rank(ref: str, cell_type: str) -> pd.Series | None:
    """Mean logFC across protocols (≥ MIN_PROTOCOLS); higher = organoid up."""
    deg_dir = DIRS["02_tab"] / f"{ref}_deg_per_celltype_protocol"
    frames = []
    for proto in PROTOCOLS_3:
        tsv = deg_dir / f"{cell_type}__{proto}" / "deg_full.tsv"
        if not tsv.exists():
            continue
        df = pd.read_csv(tsv, sep="\t", usecols=["gene", "logFC"])
        df["gene"] = df["gene"].astype(str).str.upper()
        df = df.rename(columns={"logFC": proto})
        frames.append(df[["gene", proto]].drop_duplicates("gene"))
    if len(frames) < MIN_PROTOCOLS:
        return None
    merged = frames[0]
    for f in frames[1:]:
        merged = merged.merge(f, on="gene", how="outer")
    score_cols = [c for c in merged.columns if c != "gene"]
    n_obs = merged[score_cols].notna().sum(axis=1)
    merged = merged.loc[n_obs >= MIN_PROTOCOLS].copy()
    merged["score"] = merged[score_cols].mean(axis=1)
    # drop ties / zeros-only noise: require some absolute signal
    merged = merged.loc[np.isfinite(merged["score"])]
    if len(merged) < 200:
        return None
    rnk = merged.set_index("gene")["score"].sort_values(ascending=False)
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    return rnk


def clean_term(term: str, max_len: int = 55) -> str:
    s = re.sub(r"^\[.*?\]\s*", "", str(term))
    s = re.sub(r"\s*\(GO:\d+\)$", "", s)
    s = re.sub(r"\s*R-HSA-\d+$", "", s)
    return s[: max_len - 1] + "…" if len(s) > max_len else s


def run_prerank(rnk: pd.Series, gene_sets: dict[str, set[str]], label: str) -> pd.DataFrame:
    import gseapy as gp

    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        permutation_num=PERM,
        min_size=15,
        max_size=500,
        threads=THREADS,
        seed=42,
        no_plot=True,
        outdir=None,
        verbose=False,
        method="permutation",
    )
    res = pre.res2d.copy()
    # normalize column names across gseapy versions
    rename = {
        "Term": "term",
        "NES": "nes",
        "ES": "es",
        "NOM p-val": "pval",
        "FDR q-val": "fdr",
        "FWER p-val": "fwer",
        "Gene %": "gene_pct",
        "Lead_genes": "lead_genes",
        "Tag %": "tag_pct",
    }
    res = res.rename(columns={k: v for k, v in rename.items() if k in res.columns})
    for c in ("nes", "es", "pval", "fdr"):
        if c in res.columns:
            res[c] = pd.to_numeric(res[c], errors="coerce")
    res.insert(0, "n_genes_ranked", int(len(rnk)))
    res.insert(0, "cell_type", label.split("|", 1)[1] if "|" in label else label)
    res.insert(0, "reference", label.split("|", 1)[0] if "|" in label else "")
    return res


def plot_nes_heatmap(df: pd.DataFrame, ref: str, out_png: Path) -> None:
    sig = df[df["fdr"] <= FDR_CUT].copy()
    if sig.empty:
        sig = df.nsmallest(40, "fdr")
    # pick terms with strongest |NES| recurrence
    score = (
        sig.assign(abs_nes=sig["nes"].abs())
        .groupby("term")
        .agg(n_ct=("cell_type", "nunique"), max_abs=("abs_nes", "max"), mean_nes=("nes", "mean"))
        .sort_values(["n_ct", "max_abs"], ascending=[False, False])
    )
    terms = score.head(TOP_N_HEAT).index.tolist()
    cts = sorted(sig["cell_type"].unique())
    mat = (
        sig[sig["term"].isin(terms)]
        .pivot_table(index="term", columns="cell_type", values="nes", aggfunc="mean")
        .reindex(index=terms, columns=cts)
    )
    if mat.empty:
        print(f"  no heatmap for {ref}", flush=True)
        return
    fig_h = max(7.0, 0.38 * len(mat) + 2.5)
    fig, ax = plt.subplots(figsize=(max(8.0, 1.3 * len(cts) + 6), fig_h))
    sns.heatmap(
        mat,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        linewidths=0.4,
        linecolor="white",
        xticklabels=True,
        yticklabels=[clean_term(t) for t in mat.index],
        cbar_kws={"label": "NES", "shrink": 0.6},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(
        f"GSEA NES vs {ref} primary (FDR ≤ {FDR_CUT})\n"
        "positive = organoid-high · negative = primary-high",
        fontsize=15,
        pad=10,
    )
    ax.tick_params(axis="x", labelsize=12, rotation=35)
    ax.tick_params(axis="y", labelsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def plot_top_bars(df: pd.DataFrame, ref: str, out_png: Path) -> None:
    """Global top pathways by mean |NES| across cell types (signed mean NES)."""
    sig = df[df["fdr"] <= FDR_CUT].copy()
    if sig.empty:
        sig = df.copy()
    agg = (
        sig.groupby("term")
        .agg(mean_nes=("nes", "mean"), max_abs=("nes", lambda s: float(np.max(np.abs(s)))), n_ct=("cell_type", "nunique"))
        .reset_index()
    )
    agg["sort"] = agg["max_abs"]
    show = agg.nlargest(TOP_N_BARS, "sort").sort_values("mean_nes")
    colors = [("#C62828" if v >= 0 else "#1565C0") for v in show["mean_nes"]]
    fig, ax = plt.subplots(figsize=(11.0, max(5.5, 0.38 * len(show) + 1.8)))
    y = np.arange(len(show))
    ax.barh(y, show["mean_nes"], color=colors, height=0.72)
    ax.axvline(0, color="#222222", lw=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([clean_term(t) for t in show["term"]], fontsize=12)
    ax.set_xlabel("Mean NES across cell types", fontsize=14)
    ax.set_title(
        f"Top GSEA pathways vs {ref} primary\n"
        "red = organoid-high · blue = primary-high",
        fontsize=15,
        pad=8,
    )
    ax.legend(
        handles=[
            Patch(facecolor="#C62828", label="Organoid-high"),
            Patch(facecolor="#1565C0", label="Primary-high"),
        ],
        loc="lower right",
        frameon=False,
        fontsize=12,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def list_cell_types(ref: str) -> list[str]:
    cons = DIRS["02_tab"] / f"{ref}_protocol_consensus_genes.csv"
    if cons.exists():
        return sorted(pd.read_csv(cons)["cell_type"].astype(str).unique())
    deg_dir = DIRS["02_tab"] / f"{ref}_deg_per_celltype_protocol"
    return sorted({p.name.split("__", 1)[0] for p in deg_dir.glob("*__*") if p.is_dir()})


def run_ref(ref: str, gene_sets: dict[str, set[str]]) -> None:
    print(f"=== GSEA prerank: {ref} ===", flush=True)
    fig_dir, tab_dir = DIRS["02_gsea_fig"], DIRS["02_gsea_tab"]
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for ct in list_cell_types(ref):
        rnk = build_rank(ref, ct)
        if rnk is None:
            print(f"  skip {ct}: insufficient ranked genes", flush=True)
            continue
        print(f"  {ct}: ranking {len(rnk):,} genes…", flush=True)
        try:
            res = run_prerank(rnk, gene_sets, f"{ref}|{ct}")
        except Exception as exc:  # noqa: BLE001
            print(f"  GSEA failed for {ct}: {exc}", flush=True)
            continue
        res.to_csv(tab_dir / f"{ref}_{ct}_gsea_prerank.csv", index=False)
        rnk.to_csv(tab_dir / f"{ref}_{ct}_ranking.rnk", sep="\t", header=False)
        n_sig = int((res["fdr"] <= FDR_CUT).sum()) if "fdr" in res.columns else 0
        print(f"    terms={len(res)}  FDR≤{FDR_CUT}={n_sig}", flush=True)
        rows.append(res)

    if not rows:
        print(f"  no GSEA results for {ref}", flush=True)
        return
    all_df = pd.concat(rows, ignore_index=True)
    all_df.to_csv(tab_dir / f"{ref}_gsea_prerank_all.csv", index=False)
    plot_nes_heatmap(all_df, ref, fig_dir / f"{ref}_gsea_nes_heatmap.png")
    plot_top_bars(all_df, ref, fig_dir / f"{ref}_gsea_top_pathways.png")


def main() -> None:
    ensure_dirs()
    DIRS["02_gsea_fig"].mkdir(parents=True, exist_ok=True)
    DIRS["02_gsea_tab"].mkdir(parents=True, exist_ok=True)
    gene_sets = load_gmt_dict()
    print(f"Loaded {len(gene_sets)} gene sets from {list(GMT_FILES)}", flush=True)
    run_ref("adult", gene_sets)
    run_ref("fetal", gene_sets)
    print("GSEA complete.", flush=True)


if __name__ == "__main__":
    main()
