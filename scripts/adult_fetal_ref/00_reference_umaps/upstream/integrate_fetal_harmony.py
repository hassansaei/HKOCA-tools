#!/usr/bin/env python
"""Build the fetal weeks 9–17 Harmony reference entirely under 00_reference_umaps/.

Reads local inputs:
  organoid_vs_primary_atlas/00_reference_umaps/inputs/fetal/{stewart,hochane,tran}.h5ad

Writes local outputs:
  objects/fetal_weeks_9_17_harmony_integrated.h5ad
  objects/fetal_stewart_hochane_tran_harmonized_preHarmony.h5ad
  figures/fetal_integration/umap_*.png
  tables/fetal_integration/*.csv + INTEGRATION_SUMMARY.md

Only cells with a numeric post-fertilization week in 9–17 are kept.
Embryonic stage, Carnegie stage, and other non-week labels are dropped.
"""
from __future__ import annotations

import gc
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import anndata as ad
import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib_config import DIRS, ROOT, apply_fig_style, ensure_dirs  # noqa: E402

apply_fig_style()

INPUT_DIR = ROOT / "00_reference_umaps/inputs/fetal"
OBJ_DIR = ROOT / "00_reference_umaps/objects"
FIG_DIR = DIRS["00_fig"] / "fetal_integration"
TAB_DIR = DIRS["00_tab"] / "fetal_integration"

# author_annotation → harmonized_celltype
CELLTYPE_HARMONIZATION = {
    # Stewart
    "Cap mesenchyme": "NPC",
    "Proliferating cap mesenchyme": "NPC",
    "Proximal renal vesicle": "Renal_vesicle_CSB",
    "Distal renal vesicle": "Renal_vesicle_CSB",
    "Proliferating distal renal vesicle": "Renal_vesicle_CSB",
    "Proximal S shaped body": "SSB",
    "Medial S shaped body": "SSB",
    "Distal S shaped body": "SSB",
    "Podocyte": "Podocytes",
    "Proximal tubule": "Proximal_tubule",
    "Loop of Henle": "Distal_tubule_LOH",
    "Proximal UB": "Ureteric_epithelium",
    "CNT/PC - proximal UB": "Ureteric_epithelium",
    "Pelvic epithelium - distal UB": "Ureteric_epithelium",
    "Stroma progenitor": "Stroma_progenitor",
    "Proliferating stroma progenitor": "Stroma_progenitor",
    "Fibroblast 1": "Fibroblast",
    "Fibroblast 2": "Fibroblast",
    "Proliferating fibroblast": "Fibroblast",
    "Myofibroblast 1": "Myofibroblast",
    "Myofibroblast 2": "Myofibroblast",
    "Proliferating myofibroblast": "Myofibroblast",
    "Endothelium": "Endothelial",
    "Macrophage 1": "Immune",
    "Macrophage 2": "Immune",
    "Proliferating macrophage": "Immune",
    "Monocyte": "Immune",
    "Proliferating monocyte": "Immune",
    "cDC1": "Immune",
    "cDC2": "Immune",
    "Proliferating cDC2": "Immune",
    "pDC": "Immune",
    "NK cell": "Immune",
    "Proliferating NK cell": "Immune",
    "CD4 T cell": "Immune",
    "CD8 T cell": "Immune",
    "B cell": "Immune",
    "Proliferating B cell": "Immune",
    "Innate like lymphocyte": "Immune",
    "Neutrophil": "Immune",
    "Mast cells": "Immune",
    "Megakaryocyte": "Megakaryocyte",
    "Erythroid": "Erythroid",
    "Neuron": "Neuron",
    # Hochane
    "NPCa: nephron progenitor cells a": "NPC",
    "NPCb: nephron progenitor cells b": "NPC",
    "NPCc: nephron progenitor cells c": "NPC",
    "NPCd: nephron progenitor cells d": "NPC",
    "PTA: pretubular aggregate": "PTA",
    "RVCSBa: renal vesicle/comma-shapedbody a": "Renal_vesicle_CSB",
    "RVCSB b: renal vesicle/comma-shapedbody b": "Renal_vesicle_CSB",
    "SSBpr: s-shaped body proximal precursor cells": "SSB",
    "SSBm/d: s-shaped body medial/distal": "SSB",
    "SSBpod: s-shaped body podocyte precursor cells": "Early_podocyte",
    "Pod: podocyte": "Podocytes",
    "ErPrT: early proximal tubule": "Proximal_tubule",
    "DTLH: distal tubule/loop of Henle": "Distal_tubule_LOH",
    "CnT: connecting tubule;": "Connecting_tubule",
    "UBCD: ureteric bud/collecting duct": "Ureteric_epithelium",
    "IPC: interstitial progenitor cell": "Stroma_progenitor",
    "ICa: interstitialcells a": "Fibroblast",
    "ICa, interstitialcells b": "Fibroblast",
    "Mes: mesangial cells": "Mesangial",
    "End: endothelial cells": "Endothelial",
    "Leu: leukocyte": "Immune",
    "Prolif: proliferating cells": "Unassigned",
    # Tran DevKidCC
    "NPC": "NPC",
    "NPC-like": "NPC",
    "EN": "Renal_vesicle_CSB",
    "EPod": "Early_podocyte",
    "Pod": "Podocytes",
    "PEC": "PEC",
    "EPT": "Proximal_tubule",
    "PT": "Proximal_tubule",
    "EDT": "Distal_tubule_LOH",
    "DT": "Distal_tubule_LOH",
    "LOH": "Distal_tubule_LOH",
    "UTip": "Ureteric_epithelium",
    "UOS": "Ureteric_epithelium",
    "UIS": "Ureteric_epithelium",
    "SPC": "Stroma_progenitor",
    "MS": "Fibroblast",
    "MesS": "Fibroblast",
    "CS": "Myofibroblast",
    "Endo": "Endothelial",
    "Leukocytes": "Immune",
    "unassigned": "Unassigned",
}


def _to_symbols(adata: ad.AnnData, dataset: str) -> ad.AnnData:
    if dataset == "stewart":
        sym = adata.var["feature_name"].astype(str)
    elif dataset == "hochane":
        if "gene_symbol" in adata.var.columns:
            sym = adata.var["gene_symbol"].astype(str)
        elif "HUGO" in adata.var.columns:
            sym = adata.var["HUGO"].astype(str)
        else:
            sym = adata.var["feature_name"].astype(str)
        bad = sym.str.startswith("ENSG") | (sym == "") | sym.isna()
        if "feature_name" in adata.var.columns:
            alt = adata.var["feature_name"].astype(str)
            sym = sym.where(~bad, alt)
    else:
        sym = pd.Index(adata.var_names.astype(str))

    adata.var["gene_symbol"] = sym.astype(str).values
    adata.var_names = adata.var["gene_symbol"].astype(str)
    adata.var_names_make_unique()
    return adata


def _ensure_counts(adata: ad.AnnData) -> ad.AnnData:
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    if sparse.issparse(adata.X):
        if adata.X.data.min() < 0:
            raise ValueError("Negative values in counts matrix")
    elif np.nanmin(adata.X) < 0:
        raise ValueError("Negative values in counts matrix")
    return adata


def parse_week_stewart(obs: pd.DataFrame) -> pd.Series:
    """Numeric post-fertilization weeks only; embryonic stage (F35) → NaN."""
    stage = obs["development_stage"].astype(str)
    week = stage.str.extract(r"(\d+)(?:st|nd|rd|th) week post-fertilization")[0].astype(float)
    week = week.where(~obs["donor_id"].astype(str).eq("F35"), other=float("nan"))
    return week


def parse_week_hochane(obs: pd.DataFrame) -> pd.Series:
    """Numeric weeks only; Carnegie stage → NaN."""
    stage = obs["week_post_fertilization"].astype(str)
    week = stage.str.extract(r"(\d+)(?:st|nd|rd|th) week post-fertilization")[0].astype(float)
    carnegie = obs["development_stage"].astype(str).str.contains("Carnegie", case=False, na=False)
    week = week.where(~carnegie, other=float("nan"))
    return week


def parse_week_tran(obs: pd.DataFrame) -> pd.Series:
    return obs["week"].astype(str).str.extract(r"H(\d+)w")[0].astype(float)


def prepare_dataset(name: str) -> ad.AnnData:
    path = INPUT_DIR / f"{name}.h5ad"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing local input: {path}\n"
            f"Place stewart.h5ad, hochane.h5ad, and tran.h5ad under {INPUT_DIR}"
        )
    print(f"Loading {name} from {path} ...", flush=True)
    adata = ad.read_h5ad(path)
    adata = _to_symbols(adata, name)
    adata = _ensure_counts(adata)

    adata.obs["dataset"] = {"stewart": "Stewart", "hochane": "Hochane", "tran": "Tran"}[name]
    adata.obs["author_annotation"] = adata.obs["author_annotation"].astype(str)
    mapped = adata.obs["author_annotation"].map(CELLTYPE_HARMONIZATION)
    missing = mapped.isna()
    if missing.any():
        miss_labs = sorted(adata.obs.loc[missing, "author_annotation"].unique())
        raise KeyError(f"{name}: labels missing from dictionary: {miss_labs}")
    adata.obs["harmonized_celltype"] = mapped.astype(str)
    adata.obs["author_annotation_original"] = adata.obs["author_annotation"].astype(str)

    if name == "stewart":
        adata.obs["week"] = parse_week_stewart(adata.obs)
    elif name == "hochane":
        adata.obs["week"] = parse_week_hochane(adata.obs)
    else:
        adata.obs["week"] = parse_week_tran(adata.obs)

    adata.obs["week"] = pd.to_numeric(adata.obs["week"], errors="coerce")
    adata.obs["batch"] = adata.obs["dataset"].astype(str)

    keep_obs = [
        "dataset",
        "batch",
        "week",
        "author_annotation",
        "author_annotation_original",
        "harmonized_celltype",
    ]
    for extra in ["donor_id", "sample_id", "zone", "LineageID", "DKCC"]:
        if extra in adata.obs.columns:
            keep_obs.append(extra)
    adata.obs = adata.obs[keep_obs].copy()
    adata.var = (
        adata.var[["gene_symbol"]].copy()
        if "gene_symbol" in adata.var
        else pd.DataFrame(index=adata.var_names)
    )
    adata.obsm = {}
    adata.varm = {}
    adata.uns = {}
    adata.obsp = {}
    adata.layers = {}
    print(
        f"  {name}: {adata.n_obs} cells, "
        f"weeks={sorted(adata.obs['week'].dropna().unique().tolist())}, "
        f"n_types={adata.obs['harmonized_celltype'].nunique()}",
        flush=True,
    )
    return adata


def main() -> None:
    ensure_dirs()
    OBJ_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TAB_DIR.mkdir(parents=True, exist_ok=True)

    sc.settings.verbosity = 2
    datasets = [prepare_dataset(n) for n in ("stewart", "hochane", "tran")]

    # Keep numeric post-fertilization weeks 9–17 only
    filtered = []
    for a in datasets:
        n0 = a.n_obs
        keep = a.obs["week"].between(9, 17, inclusive="both")
        a = a[keep].copy()
        print(f"Week filter {a.obs['dataset'].iloc[0]}: {n0} -> {a.n_obs}", flush=True)
        filtered.append(a)
    datasets = filtered

    print("Concatenating on gene symbols (inner join)...", flush=True)
    for a in datasets:
        ds = a.obs["dataset"].iloc[0]
        a.obs_names = [f"{ds}_{x}" for x in a.obs_names.astype(str)]
        a.obs_names_make_unique()

    combined = ad.concat(datasets, join="inner", merge="same", index_unique=None)
    combined.obs_names_make_unique()
    if sparse.issparse(combined.X):
        combined.X.data = np.nan_to_num(combined.X.data, nan=0.0)
    print(f"Combined: {combined.n_obs} cells × {combined.n_vars} genes", flush=True)

    pre_path = OBJ_DIR / "fetal_stewart_hochane_tran_harmonized_preHarmony.h5ad"
    combined.write_h5ad(pre_path)
    print("Wrote", pre_path, flush=True)

    print("Normalize / HVG / PCA / Harmony...", flush=True)
    combined.layers["counts"] = combined.X.copy()
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    sc.pp.highly_variable_genes(combined, n_top_genes=3000, batch_key="batch", flavor="seurat")
    combined.raw = None
    combined = combined[:, combined.var["highly_variable"]].copy()
    sc.pp.scale(combined, max_value=10)
    sc.tl.pca(combined, n_comps=50, svd_solver="arpack")

    ho = hm.run_harmony(combined.obsm["X_pca"], combined.obs, ["batch"])
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != combined.n_obs and Z.shape[1] == combined.n_obs:
        Z = Z.T
    if Z.shape != (combined.n_obs, combined.obsm["X_pca"].shape[1]):
        raise ValueError(f"Unexpected Harmony output shape {Z.shape}")
    combined.obsm["X_pca_harmony"] = Z

    sc.pp.neighbors(combined, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=40)
    sc.tl.umap(combined)
    sc.tl.leiden(combined, resolution=0.6, key_added="leiden_res0.6")

    out_path = OBJ_DIR / "fetal_weeks_9_17_harmony_integrated.h5ad"
    combined.write_h5ad(out_path)
    print("Wrote", out_path, flush=True)

    sc.settings.figdir = str(FIG_DIR)
    for color in ["dataset", "week", "harmonized_celltype", "leiden_res0.6"]:
        sc.pl.umap(combined, color=color, save=f"_{color}.png", show=False, frameon=False)

    combined.obs.assign(n=1).groupby(["week", "dataset"], observed=False)["n"].count().unstack(
        fill_value=0
    ).to_csv(TAB_DIR / "cells_by_week_dataset.csv")
    combined.obs["harmonized_celltype"].value_counts().to_csv(
        TAB_DIR / "harmonized_celltype_counts.csv"
    )
    pd.crosstab(combined.obs["harmonized_celltype"], combined.obs["dataset"]).to_csv(
        TAB_DIR / "harmonized_celltype_by_dataset.csv"
    )
    pd.crosstab(combined.obs["harmonized_celltype"], combined.obs["week"]).to_csv(
        TAB_DIR / "harmonized_celltype_by_week.csv"
    )
    pd.Series(CELLTYPE_HARMONIZATION, name="harmonized_celltype").rename_axis(
        "author_annotation"
    ).reset_index().to_csv(TAB_DIR / "celltype_harmonization_dictionary.csv", index=False)

    summary = f"""# Fetal weeks 9–17 Harmony integration

## Inputs (local)
- `{INPUT_DIR}/stewart.h5ad`
- `{INPUT_DIR}/hochane.h5ad`
- `{INPUT_DIR}/tran.h5ad`
- Author labels remapped → `harmonized_celltype`
- Weeks kept: **9–17** (numeric post-fertilization weeks only; embryonic / Carnegie / other non-week labels excluded)

## Cells
- Total: **{combined.n_obs:,}**
- Genes (inner join symbols → HVG): **{combined.n_vars:,}** (3000 HVGs for PCA/Harmony)
- Batch key: `dataset` (Stewart / Hochane / Tran)
- Embedding: `X_pca_harmony` → neighbors → UMAP

## Outputs (local)
- `{out_path}`
- `{pre_path}`
- `{FIG_DIR}/umap_*.png`
- `{TAB_DIR}/*.csv`
"""
    (TAB_DIR / "INTEGRATION_SUMMARY.md").write_text(summary)
    print(summary, flush=True)

    prov = {
        "inputs": {
            "Stewart": str(INPUT_DIR / "stewart.h5ad"),
            "Hochane": str(INPUT_DIR / "hochane.h5ad"),
            "Tran": str(INPUT_DIR / "tran.h5ad"),
        },
        "weeks_kept": "9–17 (numeric post-fertilization weeks only)",
        "excluded_stages": "embryonic stage / Carnegie stage / any cells without a week number in 9–17",
        "integration": "PCA → Harmony (batch=dataset) → UMAP + Leiden 0.6",
        "n_cells": int(combined.n_obs),
        "n_genes_hvg": int(combined.n_vars),
        "object": str(out_path),
        "figures": str(FIG_DIR),
        "tables": str(TAB_DIR),
        "driver": str(Path(__file__).resolve()),
    }
    (DIRS["00_tab"] / "fetal_reference_provenance.json").write_text(json.dumps(prov, indent=2))
    print("Wrote", DIRS["00_tab"] / "fetal_reference_provenance.json", flush=True)
    gc.collect()
    print("Fetal Harmony integration done.", flush=True)


if __name__ == "__main__":
    main()
