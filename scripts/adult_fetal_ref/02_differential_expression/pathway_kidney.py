#!/usr/bin/env python
"""Kidney pathway ORA filter + diverging-bar figures (per cell type).

Keeps pathways that match kidney / transport / metabolism names, or belong to
kidney GO ontology descendants. Uses Shared Up/Down consensus DEGs.

Outputs under 02_differential_expression/pathways/:
  figures/{adult,fetal}_pathway_kidney_diverging_bars_{cell_type}.png
  tables/...
"""
from __future__ import annotations

import re
import sys
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import obonet
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import hypergeom

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib_config import DIRS, PATHS, ROOT, apply_fig_style  # noqa: E402

apply_fig_style()

FIG = ROOT / "02_differential_expression" / "pathways" / "figures"
TAB = ROOT / "02_differential_expression" / "pathways" / "tables"
OBO_PATH = TAB / "go-basic.obo"
GMT_PATH = PATHS["gmt_dir"] / "Enrichr.GO_Biological_Process_2023.gmt"
P_THRESH = 0.05
UP, DN = "#C62828", "#1565C0"
UP_NS, DN_NS = "#EF9A9A", "#90CAF9"

KIDNEY_GO_PARENTS = [
    "GO:0001822",  # kidney development
    "GO:0072001",  # renal system development
    "GO:0001655",  # urogenital system development
    "GO:0003094",  # glomerular filtration
    "GO:0097205",  # renal filtration
    "GO:0003097",  # renal water transport
    "GO:0097254",  # renal tubular secretion
    "GO:0003093",  # regulation of glomerular filtration
    "GO:0032835",  # glomerulus development
    "GO:0072080",  # nephron tubule development
    "GO:0072073",  # kidney epithelium development
    "GO:0072012",  # glomerulus vasculature development
    "GO:0072310",  # glomerular epithelial cell development
    "GO:0072311",  # glomerular epithelial cell differentiation
    "GO:0090521",  # podocyte cell migration
    "GO:0072112",  # podocyte differentiation
    "GO:0070293",  # renal absorption
]

NAME_RE = re.compile(
    r"kidney|renal|nephron|glomerul|podocyt|ureteric|collecting duct|"
    r"proximal tubul|distal tubul|loop of henle|thick ascending|"
    r"juxtaglomerul|renin|aldosterone|aquaporin|metaneph|"
    r"organic cation|organic anion|zwitterion transport|"
    r"bicarbonate reclamation|sodium reabsorb|potassium secretion|"
    r"water reabsorb|nephropathy|nephritis|polycystic kidney|"
    r"vitamin d biosynthesis|vitamin d metabolism|calcium reabsorb|"
    r"slc transporter|slc[- ]mediated|amino acid reabsorb|phosphate reabsorb|"
    r"inorganic (?:cation|anion) (?:import|transport|transmembrane)|"
    r"transport of inorganic cations|"
    r"transmembrane transporter|carboxylic acid transmembrane|"
    r"modified amino acid transport|nitrogen compound transport|"
    r"amino acid (?:catabolic|metabolic|transport)|"
    r"fatty acid (?:metabol|degradation|beta)|"
    r"oxidative phosphoryl|xenobiotic metabol|bile acid metabol|"
    r"citric acid|tca cycle|biological oxidat|drug metabol|"
    r"monocarboxylic acid metabol|pyruvate metabol|propanoate metabol|"
    r"tryptophan metabol|arginine and proline|histidine metabol|"
    r"beta-alanine metabol|glycine, serine|glyoxylate|"
    r"alanine, aspartate|metabolism of amino acids|"
    r"metabolism of vitamins|peroxisomal lipid|respiratory electron|"
    r"disorders of transmembrane",
    re.I,
)


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


def ensure_obo(path: Path) -> None:
    if path.exists() and path.stat().st_size > 1_000_000:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading go-basic.obo → {path}", flush=True)
    urls = [
        "http://current.geneontology.org/ontology/go-basic.obo",
        "https://purl.obolibrary.org/obo/go/go-basic.obo",
    ]
    last_err: Exception | None = None
    for url in urls:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=180) as resp, open(path, "wb") as out:
                out.write(resp.read())
            if path.stat().st_size > 1_000_000:
                return
        except Exception as exc:  # noqa: BLE001
            last_err = exc
    raise RuntimeError(f"Failed to download go-basic.obo: {last_err}")


def kidney_go_ids(obo_path: Path) -> set[str]:
    graph = obonet.read_obo(obo_path)
    union: set[str] = set()
    for parent in KIDNEY_GO_PARENTS:
        if parent not in graph:
            continue
        union |= set(nx.ancestors(graph, parent)) | {parent}
    return union


def load_go_gmt(gmt_path: Path) -> dict[str, tuple[str, set[str]]]:
    out: dict[str, tuple[str, set[str]]] = {}
    with open(gmt_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            m = re.search(r"(GO:\d+)", parts[0])
            if not m:
                continue
            genes = {g.strip().upper() for g in parts[2:] if g.strip()}
            if 10 <= len(genes) <= 500:
                out[m.group(1)] = (parts[0], genes)
    return out


def clean_term(term: str, max_len: int = 52) -> str:
    s = re.sub(r"\s*\(GO:\d+\)$", "", str(term))
    s = re.sub(r"\s*R-HSA-\d+$", "", s)
    return s[: max_len - 1] + "…" if len(s) > max_len else s


def load_query_and_bg(ref: str, cell_type: str, direction: str) -> tuple[set[str], set[str]]:
    cons = pd.read_csv(DIRS["02_tab"] / f"{ref}_protocol_consensus_genes.csv")
    q = set(
        cons.loc[
            (cons["cell_type"] == cell_type) & (cons["category"] == direction),
            "gene",
        ]
        .astype(str)
        .str.upper()
    )
    bg: set[str] = set()
    deg_dir = DIRS["02_tab"] / f"{ref}_deg_per_celltype_protocol"
    for p in deg_dir.glob(f"{cell_type}__*/deg_full.tsv"):
        bg.update(pd.read_csv(p, sep="\t", usecols=["gene"])["gene"].astype(str).str.upper())
    return q, bg


def run_ora(
    query: set[str],
    background: set[str],
    gene_sets: dict[str, tuple[str, set[str]]],
) -> pd.DataFrame:
    q = {g.upper() for g in query} & background
    M, N = len(background), len(q)
    if N < 5 or M < 50:
        return pd.DataFrame()
    rows, pvals = [], []
    for gid, (term, genes) in gene_sets.items():
        t_in = genes & background
        n_t = len(t_in)
        if n_t < 5:
            continue
        k = len(q & t_in)
        p = float(hypergeom.sf(k - 1, M, n_t, N)) if k > 0 else 1.0
        pvals.append(p)
        rows.append(
            {
                "goid": gid,
                "term": term,
                "overlap": k,
                "term_size_in_bg": n_t,
                "query_size": N,
                "bg_size": M,
                "fold_enrichment": (k / N) / (n_t / M) if N and n_t else 0.0,
                "pval": p,
                "genes": ";".join(sorted(q & t_in)[:25]),
            }
        )
    if not rows:
        return pd.DataFrame()
    for r, f in zip(rows, bh_adjust(np.array(pvals, float))):
        r["fdr"] = float(f)
    return pd.DataFrame(rows).sort_values(["pval", "fdr"])


def filter_kidney(ora_df: pd.DataFrame, ontology_goids: set[str]) -> pd.DataFrame:
    term = ora_df["term"].astype(str)
    name_keep = term.str.contains(NAME_RE)
    if ontology_goids:
        ont_pat = "|".join(re.escape(g) for g in sorted(ontology_goids))
        ont_keep = term.str.contains(ont_pat, regex=True)
    else:
        ont_keep = pd.Series(False, index=ora_df.index)
    return ora_df.loc[name_keep | ont_keep].copy()


def plot_diverging(ora_df: pd.DataFrame, ref: str, cell_type: str, out_png: Path) -> pd.DataFrame:
    kid = ora_df[ora_df["cell_type"] == cell_type].copy()
    if kid.empty:
        return pd.DataFrame()
    kid["nlp"] = -np.log10(np.clip(kid["pval"].to_numpy(float), 1e-300, 1.0))
    kid["label"] = [f"[{lib}] {clean_term(t)}" for lib, t in zip(kid["library"], kid["term"])]

    rows = []
    for (lib, term, lab), g in kid.groupby(["library", "term", "label"]):
        up, dn = g[g["direction"] == "Shared_Up_3of3"], g[g["direction"] == "Shared_Down_3of3"]
        up_p = float(up["pval"].min()) if len(up) else np.nan
        dn_p = float(dn["pval"].min()) if len(dn) else np.nan
        rows.append(
            {
                "library": lib,
                "term": term,
                "label": lab,
                "cell_type": cell_type,
                "up_neglog10_p": float(up["nlp"].max()) if len(up) else 0.0,
                "down_neglog10_p": float(dn["nlp"].max()) if len(dn) else 0.0,
                "up_sig": bool(up_p < P_THRESH) if np.isfinite(up_p) else False,
                "down_sig": bool(dn_p < P_THRESH) if np.isfinite(dn_p) else False,
            }
        )
    bar = pd.DataFrame(rows)
    if bar.empty:
        return bar
    bar["sort_key"] = bar[["up_neglog10_p", "down_neglog10_p"]].max(axis=1)
    bar = bar.sort_values("sort_key").reset_index(drop=True)
    if len(bar) > 45:
        bar = bar.nlargest(45, "sort_key").sort_values("sort_key").reset_index(drop=True)

    y = np.arange(len(bar))
    fig, ax = plt.subplots(figsize=(11.5, max(5.0, 0.38 * len(bar) + 2.0)))
    for i, r in bar.iterrows():
        if r["down_neglog10_p"] > 0:
            ax.barh(i, -r["down_neglog10_p"], height=0.7, color=DN if r["down_sig"] else DN_NS, edgecolor="none")
        if r["up_neglog10_p"] > 0:
            ax.barh(i, r["up_neglog10_p"], height=0.7, color=UP if r["up_sig"] else UP_NS, edgecolor="none")
    ax.axvline(0, color="#222222", lw=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(bar["label"], fontsize=12)
    ax.set_xlabel("−log10(p)  (left = Shared Down, right = Shared Up)", fontsize=14)
    ax.set_title(f"{cell_type.replace('_', ' ')} vs {ref} primary", fontsize=16, pad=8)
    vmax = max(1.5, float(bar[["up_neglog10_p", "down_neglog10_p"]].to_numpy().max()) * 1.12)
    ax.set_xlim(-vmax, vmax)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        handles=[
            Patch(facecolor=UP, label="Shared Up"),
            Patch(facecolor=DN, label="Shared Down"),
            Patch(facecolor=UP_NS, label="n.s."),
        ],
        loc="lower right",
        frameon=False,
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png} ({len(bar)} pathways)", flush=True)
    return bar


def load_full_ora(ref: str) -> pd.DataFrame:
    for cand in (
        TAB / f"{ref}_pathway_ora_all.csv",
        DIRS["02_tab"] / f"{ref}_pathway_enrichment.csv",
    ):
        if cand.exists():
            df = pd.read_csv(cand)
            if cand.name.endswith("pathway_enrichment.csv"):
                df.to_csv(TAB / f"{ref}_pathway_ora_all.csv", index=False)
            return df
    raise FileNotFoundError(f"Missing pathway ORA table for {ref}")


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    ensure_obo(OBO_PATH)

    print("Kidney GO descendants…", flush=True)
    union = kidney_go_ids(OBO_PATH)
    gmt = load_go_gmt(GMT_PATH)
    kidney_sets = {gid: gmt[gid] for gid in union if gid in gmt}
    ontology_goids = set(kidney_sets)
    print(f"  ontology={len(union)}  with GMT genes={len(kidney_sets)}", flush=True)

    all_rows = []
    for ref in ("adult", "fetal"):
        cons = pd.read_csv(DIRS["02_tab"] / f"{ref}_protocol_consensus_genes.csv")
        cell_types = sorted(cons["cell_type"].astype(str).unique())
        print(f"\n=== {ref}: {cell_types} ===", flush=True)

        # Dedicated kidney-GO ORA (keeps sparse terms missing from genome-wide ORA)
        ont_parts = []
        for ct in cell_types:
            for direction in ("Shared_Up_3of3", "Shared_Down_3of3"):
                q, bg = load_query_and_bg(ref, ct, direction)
                ora = run_ora(q, bg, kidney_sets)
                if ora.empty:
                    continue
                ora.insert(0, "direction", direction)
                ora.insert(0, "cell_type", ct)
                ora.insert(0, "reference", ref)
                ora["library"] = "GO_BP"
                ont_parts.append(ora)
        from_ont = pd.concat(ont_parts, ignore_index=True) if ont_parts else pd.DataFrame()

        ora_full = load_full_ora(ref)
        from_full = filter_kidney(ora_full, ontology_goids)

        if not from_ont.empty and not from_full.empty:
            for c in from_full.columns:
                if c not in from_ont.columns:
                    from_ont[c] = pd.NA
            from_ont = from_ont.reindex(columns=from_full.columns)
            keys = set(
                zip(
                    from_full["cell_type"].astype(str),
                    from_full["direction"].astype(str),
                    from_full["library"].astype(str),
                    from_full["term"].astype(str),
                )
            )
            only_ont = from_ont[
                ~from_ont.apply(
                    lambda r: (str(r["cell_type"]), str(r["direction"]), str(r["library"]), str(r["term"]))
                    in keys,
                    axis=1,
                )
            ]
            merged = pd.concat([from_full, only_ont], ignore_index=True)
        elif not from_full.empty:
            merged = from_full
        else:
            merged = from_ont

        merged.to_csv(TAB / f"{ref}_pathway_kidney_ora_rows.csv", index=False)
        all_rows.append(merged)
        print(f"  kidney rows={len(merged)}", flush=True)

        for p in FIG.glob(f"{ref}_pathway_kidney_diverging_bars_*.png"):
            p.unlink()
        for p in TAB.glob(f"{ref}_pathway_kidney_diverging_bars_*.csv"):
            p.unlink()
        for ct in sorted(merged["cell_type"].astype(str).unique()):
            bar = plot_diverging(
                merged, ref, ct, FIG / f"{ref}_pathway_kidney_diverging_bars_{ct}.png"
            )
            if not bar.empty:
                bar.to_csv(TAB / f"{ref}_pathway_kidney_diverging_bars_{ct}.csv", index=False)

    if all_rows:
        pd.concat(all_rows, ignore_index=True).to_csv(TAB / "pathway_kidney_ora_all.csv", index=False)
    print("Kidney pathways complete.", flush=True)


if __name__ == "__main__":
    main()
