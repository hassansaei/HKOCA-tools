# ==============================================================================
# HKOCA INTEGRATION BENCHMARKING — FAST METRICS EVALUATION
# Configuration: benchmark_config.yaml
# ==============================================================================

import os
import gc
import time
import logging
import psutil
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import scib
except ImportError:
    import subprocess, sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scib", "-q"])
    import scib

try:
    import yaml
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyyaml", "-q"])
    import yaml

from pathlib import Path


def load_config():
    argv = sys.argv[1:]
    if "--config" in argv:
        config_file = Path(argv[argv.index("--config") + 1]).resolve()
    else:
        config_file = Path(__file__).resolve().parent / "benchmark_config.yaml"
    with open(config_file, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def get_iter_id(cfg):
    return int(os.environ.get(cfg["iteration"]["slurm_env_var"], cfg["iteration"]["default_iter_id"]))


def get_checkpoint_dir(cfg, iter_id):
    pattern = cfg["iteration"]["checkpoint_subdir_pattern"]
    return Path(cfg["paths"]["checkpoints_dir"]) / pattern.format(iter_id=iter_id)


CFG = load_config()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

ITER_ID = get_iter_id(CFG)
CHECKPOINT_DIR = get_checkpoint_dir(CFG, ITER_ID)

METRICS_CFG = CFG["metrics"]
INT_CFG = CFG["integration"]
BATCH_KEY = CFG["metadata"]["batch_key"]
LABEL_KEY = CFG["metadata"]["label_key"]
BIO_KEYS = METRICS_CFG["bio_keys"]
BATCH_KEYS = METRICS_CFG["batch_keys"]
METHODS = METRICS_CFG["methods"]
EMBED_KEY = METRICS_CFG["embed_key"]
CLUSTER_KEY = METRICS_CFG["cluster_key"]
RANDOM_STATE = INT_CFG["random_state"]
N_NEIGHBORS = INT_CFG["n_neighbors"]

def checkpoint_path(filename):
    return os.path.join(CHECKPOINT_DIR, filename)

def overall_score(row):
    bio = np.nanmean([row.get(k, np.nan) for k in BIO_KEYS])
    batch = np.nanmean([row.get(k, np.nan) for k in BATCH_KEYS])
    return METRICS_CFG["bio_weight"] * bio + METRICS_CFG["batch_weight"] * batch

logger.info("PHASE 2.5: Importing R embeddings")
r_methods = CFG["r_integration"]["embedding_files"]

ref_path = checkpoint_path(METRICS_CFG["reference_checkpoint"])
master_adata = ad.read_h5ad(ref_path, backed='r')
master_obs = master_adata.obs.copy()
del master_adata

for method_name, csv_file in r_methods.items():
    csv_path = checkpoint_path(csv_file)
    if os.path.exists(csv_path):
        logger.info(f"[{method_name}] Loading coordinates from CSV...")
        emb_df = pd.read_csv(csv_path, index_col=0)
        ckpt = ad.AnnData(obs=master_obs.copy())
        ckpt.obsm[EMBED_KEY] = emb_df.reindex(ckpt.obs_names).values
        sc.pp.neighbors(ckpt, use_rep=EMBED_KEY, n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)
        ckpt.write_h5ad(checkpoint_path(f"checkpoint_{method_name}.h5ad"), compression=INT_CFG["checkpoint_compression"])
        del ckpt, emb_df; gc.collect()

logger.info("PHASE 3: Computing scIB metrics")
adata_pre = ad.read_h5ad(ref_path)
all_metrics = []

for method in METHODS:
    ckpt_file = checkpoint_path(f"checkpoint_{method}.h5ad")
    if not os.path.exists(ckpt_file):
        continue

    logger.info(f"\n  Computing metrics for: {method}")
    t_method = time.time()
    ckpt = ad.read_h5ad(ckpt_file)
    
    r = {"method": method, "seed": METRICS_CFG["seed_offset"] + ITER_ID}
    
    embed = EMBED_KEY
    has_emb = embed in ckpt.obsm
    
    scib.metrics.cluster_optimal_resolution(
        ckpt,
        label_key=LABEL_KEY,
        cluster_key=CLUSTER_KEY,
        cluster_function=sc.tl.leiden,
        metric=scib.metrics.nmi
    )
    
    r["ARI"] = scib.me.ari(ckpt, cluster_key=CLUSTER_KEY, label_key=LABEL_KEY)
    r["NMI"] = scib.me.nmi(ckpt, cluster_key=CLUSTER_KEY, label_key=LABEL_KEY)
    
    if has_emb:
        r["ASW_celltype"] = scib.me.silhouette(ckpt, label_key=LABEL_KEY, embed=embed)
        r["iso_label_F1"] = scib.me.isolated_labels_f1(ckpt, label_key=LABEL_KEY, batch_key=BATCH_KEY, embed=embed, verbose=False)
        r["iso_label_ASW"] = scib.me.isolated_labels_asw(ckpt, label_key=LABEL_KEY, batch_key=BATCH_KEY, embed=embed, verbose=False)
        r["cLISI"] = scib.me.clisi_graph(
            ckpt, label_key=LABEL_KEY, type_="embed", use_rep=embed,
            k0=METRICS_CFG["lisi_k0"], subsample=METRICS_CFG["lisi_subsample"],
            scale=METRICS_CFG["lisi_scale"], verbose=False
        )
        r["ASW_batch"] = scib.me.silhouette_batch(ckpt, batch_key=BATCH_KEY, label_key=LABEL_KEY, embed=embed, verbose=False)
        r["iLISI"] = scib.me.ilisi_graph(
            ckpt, batch_key=BATCH_KEY, type_="embed", use_rep=embed,
            k0=METRICS_CFG["lisi_k0"], subsample=METRICS_CFG["lisi_subsample"],
            scale=METRICS_CFG["lisi_scale"], verbose=False
        )
    else:
        for k in ["ASW_celltype", "iso_label_F1", "iso_label_ASW", "cLISI", "ASW_batch", "iLISI"]: r[k] = float("nan")

    r["graph_connectivity"] = scib.me.graph_connectivity(ckpt, label_key=LABEL_KEY)
    try: r["kBET"] = scib.me.kBET(ckpt, batch_key=BATCH_KEY, label_key=LABEL_KEY, type_="embed", embed=embed)
    except: r["kBET"] = float("nan")
    
    try: r["PCR"] = scib.me.pcr_comparison(
        adata_pre, ckpt, covariate=BATCH_KEY, embed=embed,
        n_comps=METRICS_CFG["pcr_n_comps"], scale=METRICS_CFG["pcr_scale"], verbose=False
    )
    except: r["PCR"] = float("nan")

    r["overall"] = overall_score(r)
    all_metrics.append(r)
    logger.info(f"  {method} complete in {(time.time() - t_method)/60:.2f} min (Score: {r['overall']:.3f})")
    
    del ckpt; gc.collect()

metrics_df = pd.DataFrame(all_metrics).set_index("method")
metrics_df.to_csv(checkpoint_path(METRICS_CFG["output_file"]))

logger.info("Metrics evaluation complete")
logger.info(f"Output: {checkpoint_path(METRICS_CFG['output_file'])}")
