#!/usr/bin/env python3
"""Snapseed v5 hierarchy on a clustered atlas; flag clusters with score < 0.6.

If leiden_res_1.0 is missing, runs a short Leiden res=1.0 prelude first.
"""
from __future__ import annotations

import argparse
import gc
import os
import sys
import traceback
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import snapseed
import yaml

SCORE_FLAG_THRESHOLD = 0.6
CLUSTER_KEY = "leiden_res_1.0"
LEVEL_LABELS = ["Level_1", "Level_2", "Level_3"]


def log(msg: str = "") -> None:
    print(msg, flush=True)


def load_markers(path: Path) -> dict:
    with path.open() as f:
        raw = yaml.safe_load(f)
    # drop non-hierarchy helpers if present
    return {k: v for k, v in raw.items() if isinstance(v, dict) and "marker_genes" in v}


def collect_marker_genes(tree: dict) -> set[str]:
    genes: set[str] = set()

    def walk(node: dict):
        for _, v in node.items():
            if not isinstance(v, dict):
                continue
            for g in v.get("marker_genes", []) or []:
                genes.add(str(g))
            if "subtypes" in v and isinstance(v["subtypes"], dict):
                walk(v["subtypes"])

    walk(tree)
    return genes


def cluster_leiden(adata: ad.AnnData, key: str = CLUSTER_KEY, resolution: float = 1.0) -> None:
    """Normalize → HVG → PCA → neighbors → Leiden → UMAP (keeps all genes for Snapseed)."""
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    batch = "sample_id" if "sample_id" in adata.obs else None
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key=batch, flavor="seurat")
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=50)
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    del hvg
    gc.collect()
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
    sc.tl.leiden(adata, resolution=resolution, key_added=key, random_state=42)
    sc.tl.umap(adata, random_state=42)
    adata.obs[key] = adata.obs[key].astype(str)


def metrics_to_score_tables(metrics: dict, cluster_key: str) -> dict[str, pd.DataFrame]:
    out = {}
    for level, df in metrics.items():
        d = df.copy()
        # snapseed uses group_name column = cluster_key
        if cluster_key in d.columns:
            d = d.rename(columns={cluster_key: "cluster"})
        elif "group" in d.columns:
            d = d.rename(columns={"group": "cluster"})
        d["cluster"] = d["cluster"].astype(str)
        d["level"] = level
        out[level] = d
    return out


def build_cluster_summary(
    assignments: pd.DataFrame,
    metrics: dict,
    cluster_key: str,
    adata: ad.AnnData,
    threshold: float,
) -> pd.DataFrame:
    """One row / cluster with labels, per-level scores, overall_score, flag."""
    assign = assignments.copy()
    assign.index = assign.index.astype(str)
    level_cols = list(assign.columns)  # level_1, level_2, ...

    score_map = {}  # level -> cluster -> score
    class_map = {}
    expr_map = {}
    auc_map = {}
    for level, df in metrics.items():
        d = df.copy()
        ccol = cluster_key if cluster_key in d.columns else (
            "group" if "group" in d.columns else d.columns[0]
        )
        d[ccol] = d[ccol].astype(str)
        score_map[level] = d.set_index(ccol)["score"].to_dict() if "score" in d.columns else {}
        class_map[level] = d.set_index(ccol)["class"].to_dict() if "class" in d.columns else {}
        expr_map[level] = d.set_index(ccol)["expr"].to_dict() if "expr" in d.columns else {}
        auc_map[level] = d.set_index(ccol)["auc"].to_dict() if "auc" in d.columns else {}

    sizes = adata.obs[cluster_key].astype(str).value_counts()

    rows = []
    for cl in sorted(assign.index, key=lambda x: int(x) if str(x).isdigit() else str(x)):
        labels = {LEVEL_LABELS[i] if i < len(LEVEL_LABELS) else f"Level_{i+1}": assign.loc[cl, c]
                  for i, c in enumerate(level_cols)}
        # fill Level_3_latest style
        latest = "Unassigned"
        latest_level = None
        for i, c in enumerate(level_cols):
            lab = assign.loc[cl, c]
            if pd.notna(lab) and str(lab) not in {"", "Unassigned", "nan", "None"}:
                latest = str(lab)
                latest_level = c  # e.g. level_2

        # scores per level
        scores = {}
        for i, c in enumerate(level_cols):
            lvl_key = c  # level_1
            scv = score_map.get(lvl_key, {}).get(cl, np.nan)
            scores[f"score_{LEVEL_LABELS[i] if i < len(LEVEL_LABELS) else c}"] = scv
            scores[f"auc_{LEVEL_LABELS[i] if i < len(LEVEL_LABELS) else c}"] = auc_map.get(lvl_key, {}).get(cl, np.nan)
            scores[f"expr_{LEVEL_LABELS[i] if i < len(LEVEL_LABELS) else c}"] = expr_map.get(lvl_key, {}).get(cl, np.nan)

        overall = score_map.get(latest_level, {}).get(cl, np.nan) if latest_level else np.nan
        # fallback: level_1 score
        if pd.isna(overall):
            overall = score_map.get(level_cols[0], {}).get(cl, np.nan)

        # also keep min score across levels that have a real assignment
        avail = []
        for c in level_cols:
            lab = assign.loc[cl, c]
            if pd.notna(lab) and str(lab) not in {"", "Unassigned", "nan", "None"}:
                s = score_map.get(c, {}).get(cl, np.nan)
                if pd.notna(s):
                    avail.append(float(s))
        min_score = float(np.min(avail)) if avail else np.nan

        flag = bool(pd.notna(overall) and float(overall) < threshold)
        rows.append(
            {
                "cluster": cl,
                "n_cells": int(sizes.get(cl, 0)),
                "pct_cells": 100.0 * float(sizes.get(cl, 0)) / adata.n_obs,
                **labels,
                "Level_3_latest": latest,
                "overall_score": float(overall) if pd.notna(overall) else np.nan,
                "min_assigned_score": min_score,
                "low_confidence_flag": flag,
                "flag_reason": (
                    f"overall_score={overall:.3f} < {threshold}"
                    if flag
                    else ("missing_score" if pd.isna(overall) else "")
                ),
                **scores,
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["low_confidence_flag", "overall_score", "cluster"],
        ascending=[False, True, True],
    )


def plot_annotation_umap(adata: ad.AnnData, out_dir: Path, max_points: int = 150_000) -> None:
    if "X_umap" not in adata.obsm:
        log("  no X_umap — skip annotation UMAP plots")
        return
    n = adata.n_obs
    if n > max_points:
        rng = np.random.default_rng(42)
        idx = np.sort(rng.choice(n, size=max_points, replace=False))
        ad_plot = adata[idx].copy()
        log(f"  plotting subsample {max_points:,}/{n:,}")
    else:
        ad_plot = adata

    colors = [c for c in [CLUSTER_KEY, "Level_1", "Level_2", "Level_3_latest", "snapseed_low_confidence"] if c in ad_plot.obs]
    fig = sc.pl.umap(ad_plot, color=colors, ncols=2, wspace=0.4, show=False, return_fig=True)
    fig.savefig(out_dir / "umap_snapseed_annotation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    del ad_plot
    gc.collect()


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--input",
        default="/data-master/pure-workspace/labss/hmami/atlas_shortcut/concatenated_filtered_singlets_important_columns.h5ad",
    )
    p.add_argument(
        "--markers",
        default="/data-master/pure-workspace/labss/hmami/atlas_shortcut_new_start/Organization/annotation/markers/snapseed_markers_v5.yaml",
    )
    p.add_argument(
        "--outdir",
        default="/data-master/pure-workspace/labss/hmami/atlas_shortcut/annotation_snapseed_v5",
    )
    p.add_argument("--cluster-key", default=CLUSTER_KEY)
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--score-threshold", type=float, default=SCORE_FLAG_THRESHOLD)
    p.add_argument("--method", default="auroc", choices=["auroc", "trinatize"])
    p.add_argument("--cluster-only", action="store_true", help="Stop after Leiden (write clustered h5ad).")
    p.add_argument(
        "--clustered-out",
        default="",
        help="Where to write clustered object (default: outdir/clustered_leiden_res1.0.h5ad)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    started = datetime.now()
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    clustered_out = Path(args.clustered_out) if args.clustered_out else out_dir / "clustered_leiden_res1.0.h5ad"

    os.environ.setdefault("OMP_NUM_THREADS", "8")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
    os.environ.setdefault("MKL_NUM_THREADS", "8")
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")

    cluster_key = args.cluster_key
    thr = args.score_threshold

    log("=" * 72)
    log("Snapseed annotation (markers v5)")
    log(f"input   : {args.input}")
    log(f"markers : {args.markers}")
    log(f"outdir  : {out_dir}")
    log(f"cluster : {cluster_key}")
    log(f"flag if overall_score < {thr}")
    log("=" * 72)

    marker_tree = load_markers(Path(args.markers))
    marker_genes = collect_marker_genes(marker_tree)
    log(f"Hierarchy roots: {list(marker_tree.keys())}")
    log(f"Unique marker genes in YAML: {len(marker_genes)}")

    log("Loading AnnData ...")
    adata = sc.read_h5ad(args.input)
    log(f"  shape: {adata.n_obs:,} × {adata.n_vars:,}")

    if cluster_key not in adata.obs:
        log(f"  {cluster_key} missing — running Leiden res={args.resolution}")
        cluster_leiden(adata, key=cluster_key, resolution=args.resolution)
        log(f"  clusters: {adata.obs[cluster_key].nunique()} → {clustered_out}")
        adata.write_h5ad(clustered_out, compression="gzip")
        log(f"  wrote clustered checkpoint ({clustered_out.stat().st_size/1024**3:.2f} GiB)")
    else:
        adata.obs[cluster_key] = adata.obs[cluster_key].astype(str).astype("category")
        log(f"  clusters: {adata.obs[cluster_key].nunique()}")

    if args.cluster_only:
        log("cluster-only done")
        return 0

    present = sorted(g for g in marker_genes if g in adata.var_names)
    missing = sorted(g for g in marker_genes if g not in adata.var_names)
    log(f"  markers present: {len(present)}/{len(marker_genes)}")
    if missing:
        log(f"  ⚠ missing markers (will be ignored by Snapseed): {missing}")
    pd.DataFrame({"gene": present, "status": "present"}).to_csv(
        out_dir / "marker_gene_presence.csv", index=False
    )
    if missing:
        pd.DataFrame({"gene": missing, "status": "missing"}).to_csv(
            out_dir / "marker_genes_missing.csv", index=False
        )

    log("Running snapseed.annotate_hierarchy ...")
    results = snapseed.annotate_hierarchy(
        adata,
        marker_tree,
        group_name=cluster_key,
        method=args.method,
    )
    assignments = results["assignments"].copy()
    metrics = results["metrics"]
    assignments.index = assignments.index.astype(str)
    assignments = assignments.fillna("Unassigned").astype(str)
    log(f"  assignment columns: {list(assignments.columns)}")
    log("\nAssignments preview:")
    log(assignments.to_string())

    # Map to obs
    level_obs_cols = []
    for i, col in enumerate(assignments.columns):
        obs_col = LEVEL_LABELS[i] if i < len(LEVEL_LABELS) else f"Level_{i+1}"
        mapped = adata.obs[cluster_key].astype(str).map(assignments[col].to_dict())
        adata.obs[obs_col] = pd.Categorical(mapped.fillna("Unassigned").astype(str))
        level_obs_cols.append(obs_col)

    # Level_3_latest = deepest non-Unassigned
    latest = adata.obs[level_obs_cols[0]].astype(str).copy()
    for col in level_obs_cols[1:]:
        cur = adata.obs[col].astype(str)
        latest = latest.where(~cur.isin(["Unassigned", "nan", "None", ""]), cur)
    # prefer fillna chain from fine to coarse
    fill = adata.obs[level_obs_cols[-1]].astype(str)
    for col in reversed(level_obs_cols[:-1]):
        fill = fill.where(~fill.isin(["Unassigned", "nan", "None", ""]), adata.obs[col].astype(str))
    adata.obs["Level_3_latest"] = pd.Categorical(fill)

    # Cluster-level summary + flags
    summary = build_cluster_summary(assignments, metrics, cluster_key, adata, thr)
    summary.to_csv(out_dir / "snapseed_cluster_annotation_scores.csv", index=False)

    flagged = summary[summary["low_confidence_flag"]].copy()
    flagged.to_csv(out_dir / "snapseed_flagged_clusters_score_lt_0.6.csv", index=False)
    log(f"\nFlagged clusters (overall_score < {thr}): {len(flagged)} / {len(summary)}")
    if len(flagged):
        show = ["cluster", "n_cells", "Level_1", "Level_2", "Level_3_latest", "overall_score", "flag_reason"]
        show = [c for c in show if c in flagged.columns]
        log(flagged[show].to_string(index=False))

    # Map flags + overall score to cells
    score_by_cl = summary.set_index("cluster")["overall_score"].to_dict()
    flag_by_cl = summary.set_index("cluster")["low_confidence_flag"].to_dict()
    cl = adata.obs[cluster_key].astype(str)
    adata.obs["snapseed_overall_score"] = cl.map(score_by_cl).astype(float)
    adata.obs["snapseed_low_confidence"] = cl.map(flag_by_cl).fillna(False).astype(bool)
    adata.obs["snapseed_flag_reason"] = cl.map(
        summary.set_index("cluster")["flag_reason"].to_dict()
    ).astype(str)

    # Store uns
    adata.uns["snapseed_assignments"] = assignments
    # metrics may contain dataframes — store as dict of records for h5ad safety via uns string tables
    adata.uns["snapseed_metrics"] = metrics
    adata.uns["snapseed_params"] = {
        "markers_yaml": str(args.markers),
        "cluster_key": cluster_key,
        "method": args.method,
        "score_flag_threshold": thr,
        "n_marker_genes_yaml": len(marker_genes),
        "n_marker_genes_present": len(present),
        "missing_markers": missing,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }

    # Per-level metrics CSVs
    for level, df in metrics.items():
        d = df.copy()
        d.to_csv(out_dir / f"snapseed_metrics_{level}.csv", index=False)

    assignments.to_csv(out_dir / "snapseed_assignments.csv")

    # Label distribution
    for col in level_obs_cols + ["Level_3_latest"]:
        vc = adata.obs[col].astype(str).value_counts()
        vc.rename_axis(col).reset_index(name="n_cells").to_csv(
            out_dir / f"cells_per_{col}.csv", index=False
        )
        log(f"\n{col} distribution (top):")
        log(vc.head(15).to_string())

    log("\nPlotting annotation UMAPs ...")
    try:
        plot_annotation_umap(adata, out_dir)
    except Exception as e:
        log(f"  ⚠ UMAP plot failed: {e}")

    out_h5ad = out_dir / "annotated_snapseed_v5_leiden_res1.0.h5ad"
    log(f"\nWriting annotated h5ad → {out_h5ad}")
    adata.write_h5ad(out_h5ad, compression="gzip")
    log(f"  wrote {out_h5ad.name} ({out_h5ad.stat().st_size/1024**3:.2f} GiB)")

    # Markdown report
    lines = [
        f"# Snapseed annotation report (markers v5)",
        "",
        f"- Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"- Started: {started.isoformat(timespec='seconds')}",
        f"- Input: `{args.input}`",
        f"- Markers: `{args.markers}`",
        f"- Cells × genes: {adata.n_obs:,} × {adata.n_vars:,}",
        f"- Cluster key: `{cluster_key}` ({adata.obs[cluster_key].nunique()} clusters)",
        f"- Method: `{args.method}`",
        f"- Flag threshold: overall_score < **{thr}**",
        f"- Flagged clusters: **{len(flagged)} / {len(summary)}**",
        f"- Markers present: {len(present)}/{len(marker_genes)}"
        + (f" (missing: {', '.join(missing)})" if missing else ""),
        "",
        "## Flagged clusters (discuss)",
        "",
    ]
    if len(flagged) == 0:
        lines.append("_None_")
    else:
        cols = ["cluster", "n_cells", "pct_cells", "Level_1", "Level_2", "Level_3_latest", "overall_score"]
        cols = [c for c in cols if c in flagged.columns]
        lines.append("| " + " | ".join(cols) + " |")
        lines.append("| " + " | ".join(["---"] * len(cols)) + " |")
        for r in flagged.itertuples(index=False):
            vals = []
            for c in cols:
                v = getattr(r, c)
                if isinstance(v, float):
                    vals.append(f"{v:.3f}")
                else:
                    vals.append(str(v))
            lines.append("| " + " | ".join(vals) + " |")

    lines += [
        "",
        "## All clusters (sorted flagged first)",
        "",
        "| cluster | n_cells | Level_1 | Level_3_latest | overall_score | flagged |",
        "| ---: | ---: | :--- | :--- | ---: | :---: |",
    ]
    for r in summary.sort_values(["low_confidence_flag", "overall_score"], ascending=[False, True]).itertuples(index=False):
        lines.append(
            f"| {r.cluster} | {r.n_cells:,} | {getattr(r,'Level_1', '')} | {r.Level_3_latest} | "
            f"{r.overall_score:.3f} | {'YES' if r.low_confidence_flag else ''} |"
        )

    lines += [
        "",
        "## Outputs",
        "",
        f"- `{out_h5ad}`",
        f"- `{out_dir / 'snapseed_cluster_annotation_scores.csv'}`",
        f"- `{out_dir / 'snapseed_flagged_clusters_score_lt_0.6.csv'}`",
        f"- `{out_dir / 'snapseed_assignments.csv'}`",
        f"- `{out_dir / 'umap_snapseed_annotation.png'}`",
        "",
    ]
    (out_dir / "annotation_report.md").write_text("\n".join(lines))
    log(f"Wrote report: {out_dir / 'annotation_report.md'}")

    log("=" * 72)
    log("DONE")
    log(f"Flagged for discussion: {sorted(flagged['cluster'].tolist())}")
    log("=" * 72)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(1)
