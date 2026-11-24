#!/usr/bin/env python
"""Stage 2b: Pathway Over-Representation Analysis (ORA) on Shared Up and Shared Down
DEGs across the 3 protocols (Morizane, Takasato, Freedman), separately for Adult and Fetal.

Uses local Enrichr GMT libraries (GO_Biological_Process_2023, MSigDB_Hallmark_2020,
Reactome_2022, KEGG_2021_Human) with exact hypergeometric test + BH FDR correction.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib_config import DIRS, PATHS, apply_fig_style, ensure_dirs  # noqa: E402

apply_fig_style()

GMT_FILES = {
    "Hallmark": "Enrichr.MSigDB_Hallmark_2020.gmt",
    "GO_BP": "Enrichr.GO_Biological_Process_2023.gmt",
    "Reactome": "Enrichr.Reactome_2022.gmt",
    "KEGG": "Enrichr.KEGG_2021_Human.gmt",
}


def load_gmts() -> dict[str, dict[str, set[str]]]:
    libs = {}
    for lib_name, fname in GMT_FILES.items():
        p = PATHS["gmt_dir"] / fname
        if not p.exists():
            continue
        gs = {}
        with open(p, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                term = parts[0].strip()
                genes = {g.strip().upper() for g in parts[2:] if g.strip()}
                if 10 <= len(genes) <= 500:
                    gs[term] = genes
        libs[lib_name] = gs
    return libs


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, ranked[i] * n / (i + 1))
        adj[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0.0, 1.0)
    return out


def run_ora(
    query_genes: set[str],
    bg_genes: set[str],
    gmt_libs: dict[str, dict[str, set[str]]],
) -> pd.DataFrame:
    q = {g.upper() for g in query_genes} & bg_genes
    M = len(bg_genes)
    N = len(q)
    if N < 10 or M < 100:
        return pd.DataFrame()

    rows = []
    for lib_name, terms in gmt_libs.items():
        lib_rows = []
        pvals = []
        for term, tgenes in terms.items():
            t_in_bg = tgenes & bg_genes
            n_t = len(t_in_bg)
            if n_t < 10:
                continue
            ov = q & t_in_bg
            k = len(ov)
            if k < 3:
                continue
            p = float(hypergeom.sf(k - 1, M, n_t, N))
            fe = (k / N) / (n_t / M)
            pvals.append(p)
            lib_rows.append(
                {
                    "library": lib_name,
                    "term": term,
                    "overlap": k,
                    "term_size_in_bg": n_t,
                    "query_size": N,
                    "bg_size": M,
                    "fold_enrichment": fe,
                    "pval": p,
                    "genes": ";".join(sorted(ov)[:20]),
                }
            )
        if lib_rows:
            fdrs = bh_adjust(np.array(pvals, dtype=float))
            for r, f in zip(lib_rows, fdrs):
                r["fdr"] = float(f)
            rows.extend(lib_rows)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["fdr", "pval"])


def clean_term_name(term: str, max_len: int = 48) -> str:
    s = re.sub(r"\s*\(GO:\d+\)$", "", str(term))
    s = re.sub(r"\s*R-HSA-\d+$", "", s)
    if len(s) > max_len:
        s = s[: max_len - 1].rstrip() + "…"
    return s


def plot_pathway_dotplot(ora_df: pd.DataFrame, ref_name: str, out_png: Path) -> None:
    """Two-panel dotplot: Left = Shared Up (3/3) pathways, Right = Shared Down (3/3) pathways."""
    # Prefer Hallmark + GO_BP for crisp biological story
    sub_all = ora_df[ora_df["library"].isin(["Hallmark", "GO_BP", "Reactome"])].copy()
    sub_all = sub_all[sub_all["fdr"] < 0.05]

    fig, axes = plt.subplots(1, 2, figsize=(21, 9.5))

    for ax, direction, cmap, title in [
        (axes[0], "Shared_Up_3of3", "Reds", f"A. Enriched in Shared Up DEGs (organoid > {ref_name})"),
        (axes[1], "Shared_Down_3of3", "Blues", f"B. Enriched in Shared Down DEGs (organoid < {ref_name})"),
    ]:
        sub = sub_all[sub_all["direction"] == direction].copy()
        if sub.empty:
            ax.text(0.5, 0.5, "No enriched pathways (FDR < 0.05)", ha="center", va="center", fontsize=14)
            ax.set_title(title, fontsize=15)
            ax.axis("off")
            continue

        # Select top pathways: recurring across cell types or strongest FDR (prioritize Hallmark/GO_BP)
        sub["neg_log10_fdr"] = -np.log10(np.clip(sub["fdr"].to_numpy(), 1e-30, 1.0))
        term_score = (
            sub.groupby("term")
            .agg(n_ct=("cell_type", "nunique"), max_sig=("neg_log10_fdr", "max"))
            .sort_values(["n_ct", "max_sig"], ascending=[False, False])
        )
        top_terms = term_score.head(16).index.tolist()
        plot_df = sub[sub["term"].isin(top_terms)].copy()
        cts = sorted(plot_df["cell_type"].unique().tolist())

        # Order terms by peak cell type then score
        term_order = list(reversed(top_terms))
        x_map = {c: i for i, c in enumerate(cts)}
        y_map = {t: i for i, t in enumerate(term_order)}

        xs = [x_map[c] for c in plot_df["cell_type"]]
        ys = [y_map[t] for t in plot_df["term"]]
        sizes = np.clip(plot_df["fold_enrichment"].to_numpy() * 35.0, 40.0, 320.0)
        colors = np.clip(plot_df["neg_log10_fdr"].to_numpy(), 1.3, 15.0)

        sc = ax.scatter(
            xs,
            ys,
            s=sizes,
            c=colors,
            cmap=cmap,
            vmin=1.3,
            vmax=15.0,
            edgecolors="#333333",
            linewidths=0.5,
            alpha=0.92,
        )
        ax.set_xticks(range(len(cts)))
        ax.set_xticklabels(cts, rotation=35, ha="right", fontsize=12)
        ax.set_yticks(range(len(term_order)))
        ax.set_yticklabels([clean_term_name(t) for t in term_order], fontsize=11.5)
        ax.set_title(title, fontsize=15, pad=10)
        ax.grid(True, linestyle=":", alpha=0.35)
        ax.spines[["top", "right"]].set_visible(False)

        cbar = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
        cbar.set_label("-log10(FDR)", fontsize=12)
        cbar.ax.tick_params(labelsize=11)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)


def run_for_reference(ref_name: str, gmt_libs: dict[str, dict[str, set[str]]]) -> None:
    print(f"=== Stage 2b: Pathway ORA for {ref_name} ===", flush=True)
    cons_path = DIRS["02_tab"] / f"{ref_name}_protocol_consensus_genes.csv"
    deg_dir = DIRS["02_tab"] / f"{ref_name}_deg_per_celltype_protocol"
    if not cons_path.exists():
        print(f"  Missing {cons_path}, skipping.", flush=True)
        return

    cons_df = pd.read_csv(cons_path)
    all_rows = []
    for ct, grp in cons_df.groupby("cell_type"):
        # Background = union of tested genes for this cell type
        bg_genes = set()
        for p_dir in deg_dir.glob(f"{ct}__*"):
            f_tsv = p_dir / "deg_full.tsv"
            if f_tsv.exists():
                bg_genes.update(pd.read_csv(f_tsv, sep="\t", usecols=["gene"])["gene"].astype(str).str.upper())

        for direction in ["Shared_Up_3of3", "Shared_Down_3of3"]:
            q_genes = set(grp.loc[grp["category"] == direction, "gene"].astype(str))
            ora = run_ora(q_genes, bg_genes, gmt_libs)
            if not ora.empty:
                ora.insert(0, "direction", direction)
                ora.insert(0, "cell_type", ct)
                ora.insert(0, "reference", ref_name)
                all_rows.append(ora)
                n_sig = int((ora["fdr"] < 0.05).sum())
                print(f"  {ref_name} | {ct:26s} | {direction:16s} | q={len(q_genes):4d} | sig_terms={n_sig:3d}", flush=True)

    if all_rows:
        out_df = pd.concat(all_rows, ignore_index=True)
    else:
        out_df = pd.DataFrame(
            columns=["reference", "cell_type", "direction", "library", "term", "overlap", "term_size_in_bg", "query_size", "bg_size", "fold_enrichment", "pval", "fdr", "genes"]
        )
    out_csv = DIRS["02_tab"] / f"{ref_name}_pathway_enrichment.csv"
    out_df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} ({len(out_df)} rows)", flush=True)
    plot_pathway_dotplot(out_df, ref_name, DIRS["02_fig"] / f"{ref_name}_pathway_up_down_dotplot.png")


def main() -> None:
    ensure_dirs()
    gmt_libs = load_gmts()
    print(f"Loaded GMT libraries: { {k: len(v) for k, v in gmt_libs.items()} }", flush=True)
    run_for_reference("adult", gmt_libs)
    run_for_reference("fetal", gmt_libs)
    print("Stage 2b (Pathway ORA) complete.", flush=True)


if __name__ == "__main__":
    main()
