#!/usr/bin/env python3
"""scIB metrics for enabled methods, scored against Level_2."""
from __future__ import annotations

import gc
import logging
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

try:
    import scib
except ImportError:
    import subprocess

    subprocess.check_call([sys.executable, "-m", "pip", "install", "scib", "-q"])
    import scib

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCH_DIR = os.path.dirname(SCRIPT_DIR)
with open(os.path.join(BENCH_DIR, "config.yaml"), encoding="utf-8") as _fh:
    _CFG = yaml.safe_load(_fh)

CHECKPOINTS = os.path.join(BENCH_DIR, "checkpoints")
BATCH_KEY = str(_CFG["benchmark"]["batch_key"])
LABEL_KEY = str(_CFG["benchmark"]["label_key"])
RANDOM_STATE = int(_CFG["benchmark"]["random_state"])
_METHOD_FLAGS = {str(k): bool(v) for k, v in _CFG["methods"].items()}
ALL_METHODS = [
    "pca", "seurat_v3", "seurat_v4", "scvi", "scanvi", "scpoli",
    "bbknn", "harmony", "combat", "scanorama", "css_pearson", "css_spearman",
]
METHODS = [m for m in ALL_METHODS if _METHOD_FLAGS.get(m, False)]
BIO_KEYS = ["ARI", "NMI", "ASW_celltype", "iso_label_F1", "iso_label_ASW", "cLISI"]
BATCH_KEYS = ["ASW_batch", "graph_connectivity", "iLISI", "kBET", "PCR"]


def method_enabled(name: str) -> bool:
    return bool(_METHOD_FLAGS.get(name, False))


def iter_id() -> int:
    return int(os.environ.get("BENCHMARK_ITER_ID", "1"))


def iter_dir(i: int | None = None) -> str:
    i = iter_id() if i is None else i
    return os.path.join(CHECKPOINTS, f"subsample_iter_{i}")


def seed_for_iter(i: int | None = None) -> int:
    i = iter_id() if i is None else i
    return 42 + i


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("benchmark_metrics")

ITER_ID = iter_id()
CHECKPOINT_DIR = iter_dir(ITER_ID)


def checkpoint_path(filename: str) -> str:
    return os.path.join(CHECKPOINT_DIR, filename)


def overall_score(row: dict) -> float:
    bio = np.nanmean([row.get(k, np.nan) for k in BIO_KEYS])
    batch = np.nanmean([row.get(k, np.nan) for k in BATCH_KEYS])
    return float(0.6 * bio + 0.4 * batch)


def load_obs_reference() -> tuple[pd.DataFrame, ad.AnnData]:
    """Return (obs, adata_pre) for R import + PCR. Prefer pca, else reference.h5ad."""
    for name in ("checkpoint_pca.h5ad", "reference.h5ad"):
        path = checkpoint_path(name)
        if os.path.isfile(path):
            adata = ad.read_h5ad(path)
            logger.info("Using obs/PCR reference: %s", name)
            return adata.obs.copy(), adata
    raise FileNotFoundError(
        f"Need checkpoint_pca.h5ad or reference.h5ad in {CHECKPOINT_DIR}"
    )


logger.info("Importing R embeddings (if present) | iter=%s", ITER_ID)
r_methods = {
    "seurat_v3": "seurat_v3_embeddings.csv",
    "seurat_v4": "seurat_v4_embeddings.csv",
    "css_pearson": "css_pearson_embeddings.csv",
    "css_spearman": "css_spearman_embeddings.csv",
}

master_obs, adata_pre = load_obs_reference()

for method_name, csv_file in r_methods.items():
    if not method_enabled(method_name):
        logger.info("SKIP %s — disabled in config.yaml", method_name)
        continue
    csv_path = checkpoint_path(csv_file)
    if not os.path.exists(csv_path):
        logger.warning("SKIP %s — missing %s", method_name, csv_file)
        continue
    logger.info("[%s] loading %s", method_name, csv_file)
    emb_df = pd.read_csv(csv_path, index_col=0)
    ckpt = ad.AnnData(obs=master_obs.copy())
    ckpt.obsm["X_emb"] = emb_df.reindex(ckpt.obs_names).values
    sc.pp.neighbors(ckpt, use_rep="X_emb", n_neighbors=15, random_state=RANDOM_STATE)
    ckpt.uns["method"] = method_name
    ckpt.uns["label_key_used"] = "none"
    ckpt.uns["benchmark_variant"] = "level2_org"
    ckpt.write_h5ad(checkpoint_path(f"checkpoint_{method_name}.h5ad"), compression="gzip")
    del ckpt, emb_df
    gc.collect()

logger.info("Computing scIB metrics | label=%s | batch=%s | methods=%s", LABEL_KEY, BATCH_KEY, METHODS)
all_metrics = []

for method in METHODS:
    ckpt_file = checkpoint_path(f"checkpoint_{method}.h5ad")
    if not os.path.exists(ckpt_file):
        logger.warning("SKIP %s — missing checkpoint", method)
        continue

    logger.info("Computing metrics for %s", method)
    t0 = time.time()
    ckpt = ad.read_h5ad(ckpt_file)
    if LABEL_KEY not in ckpt.obs.columns:
        raise KeyError(f"{LABEL_KEY} missing in {ckpt_file}")

    r = {"method": method, "seed": seed_for_iter(ITER_ID), "label_key": LABEL_KEY}
    embed = "X_emb"
    has_emb = embed in ckpt.obsm

    scib.metrics.cluster_optimal_resolution(
        ckpt,
        label_key=LABEL_KEY,
        cluster_key="leiden_scib",
        cluster_function=sc.tl.leiden,
        metric=scib.metrics.nmi,
    )
    r["ARI"] = scib.me.ari(ckpt, cluster_key="leiden_scib", label_key=LABEL_KEY)
    r["NMI"] = scib.me.nmi(ckpt, cluster_key="leiden_scib", label_key=LABEL_KEY)

    if has_emb:
        r["ASW_celltype"] = scib.me.silhouette(ckpt, label_key=LABEL_KEY, embed=embed)
        r["iso_label_F1"] = scib.me.isolated_labels_f1(
            ckpt, label_key=LABEL_KEY, batch_key=BATCH_KEY, embed=embed, verbose=False
        )
        r["iso_label_ASW"] = scib.me.isolated_labels_asw(
            ckpt, label_key=LABEL_KEY, batch_key=BATCH_KEY, embed=embed, verbose=False
        )
        r["cLISI"] = scib.me.clisi_graph(
            ckpt,
            label_key=LABEL_KEY,
            type_="embed",
            use_rep=embed,
            k0=90,
            subsample=100,
            scale=True,
            verbose=False,
        )
        r["ASW_batch"] = scib.me.silhouette_batch(
            ckpt, batch_key=BATCH_KEY, label_key=LABEL_KEY, embed=embed, verbose=False
        )
        r["iLISI"] = scib.me.ilisi_graph(
            ckpt,
            batch_key=BATCH_KEY,
            type_="embed",
            use_rep=embed,
            k0=90,
            subsample=100,
            scale=True,
            verbose=False,
        )
    else:
        for k in ["ASW_celltype", "iso_label_F1", "iso_label_ASW", "cLISI", "ASW_batch", "iLISI"]:
            r[k] = float("nan")

    r["graph_connectivity"] = scib.me.graph_connectivity(ckpt, label_key=LABEL_KEY)
    try:
        r["kBET"] = scib.me.kBET(
            ckpt, batch_key=BATCH_KEY, label_key=LABEL_KEY, type_="embed", embed=embed
        )
    except Exception:
        r["kBET"] = float("nan")
    try:
        r["PCR"] = scib.me.pcr_comparison(
            adata_pre, ckpt, covariate=BATCH_KEY, embed=embed, n_comps=8, scale=True, verbose=False
        )
    except Exception:
        r["PCR"] = float("nan")

    r["overall"] = overall_score(r)
    all_metrics.append(r)
    logger.info("%s done in %.2f min (overall=%.3f)", method, (time.time() - t0) / 60.0, r["overall"])
    del ckpt
    gc.collect()

if not all_metrics:
    raise RuntimeError(f"iter {ITER_ID}: no methods scored")

df = pd.DataFrame(all_metrics)
out_csv = checkpoint_path("scib_raw_metrics_full.csv")
df.to_csv(out_csv, index=False)
logger.info("Wrote %s (%d methods, label=%s)", out_csv, len(df), LABEL_KEY)
