# HKOCA Integration Benchmarking - Python Integration Methods
# Configuration: benchmark_config.yaml

import os
import gc
import random
import logging
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import psutil
import scipy.io
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

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

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def log_ram(stage=""):
    process = psutil.Process(os.getpid())
    ram_gb = process.memory_info().rss / 1024**3
    logger.info(f"RAM {stage}: {ram_gb:.2f} GB")

ITER_ID = 1
ITER_ID = get_iter_id(CFG)

SCRIPTS_DIR = CFG["paths"]["scripts_dir"]
CHECKPOINT_DIR = get_checkpoint_dir(CFG, ITER_ID)
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

MASTER_PATH = CFG["paths"]["master_atlas"]
BATCH_KEY = CFG["metadata"]["batch_key"]
LABEL_KEY = CFG["metadata"]["label_key"]
COUNTS_LAYER = CFG["metadata"]["counts_layer"]
INT_CFG = CFG["integration"]
RANDOM_STATE = INT_CFG["random_state"]
N_PCS = INT_CFG["n_pcs"]
N_NEIGHBORS = INT_CFG["n_neighbors"]
SCALE_MAX = INT_CFG["scale_max_value"]
PCA_SOLVER = INT_CFG["pca_svd_solver"]
CKPT_COMPRESSION = INT_CFG["checkpoint_compression"]
SAMPLES_PER_ITER = CFG["iteration"]["samples_per_iteration"]
N_ITERATIONS = CFG["iteration"]["n_iterations"]
SUBSAMPLE_SEED = CFG["iteration"]["subsample_seed"]
EXPORT_CFG = CFG["export"]

def checkpoint_path(filename):
    return os.path.join(CHECKPOINT_DIR, filename)

logger.info(f"Starting iteration {ITER_ID}")
logger.info(f"Checkpoint dir: {CHECKPOINT_DIR}")

def save_checkpoint(method_name, embedding, source_adata, n_components):
    """
    Build a lightweight AnnData containing only the embedding + kNN graph
    and write it to checkpoint_<method_name>.h5ad.

    Conditionally preserves the .X matrix ONLY for the 'pca' baseline so scIB
    can compute the Principal Component Regression (PCR) metric.
    """
    n_cells = source_adata.n_obs
    if embedding.shape[0] != n_cells:
        raise ValueError(
            f"[{method_name}] Embedding shape {embedding.shape} is wrong — "
            f"expected ({n_cells}, n_dims). Did you forget/double a .T?"
        )

    if method_name == "pca":
        logger.info(f"  [{method_name}] Preserving RAW counts for PCR metric...")
        ckpt = ad.AnnData(X=source_adata.layers[COUNTS_LAYER].copy(), obs=adata_meta.obs.copy())
        ckpt.var = source_adata.var.copy()
    else:
        ckpt = ad.AnnData(obs=adata_meta.obs.copy())

    ckpt.obsm["X_emb"]             = embedding.copy()
    ckpt.obsp["connectivities"]    = source_adata.obsp["connectivities"].copy()
    ckpt.obsp["distances"]         = source_adata.obsp["distances"].copy()
    ckpt.uns["neighbors"]          = source_adata.uns["neighbors"].copy()
    ckpt.uns["method"]             = method_name
    ckpt.uns["n_components"]       = n_components

    out_path = checkpoint_path(f"checkpoint_{method_name}.h5ad")
    ckpt.write_h5ad(out_path, compression=CKPT_COMPRESSION)
    logger.info(f"  [{method_name}] Checkpoint saved: {out_path}")
    
    del ckpt
    return out_path


logger.info("Loading master atlas...")
adata_full = sc.read_h5ad(MASTER_PATH)

all_samples = sorted(adata_full.obs[BATCH_KEY].unique().tolist())
total_available = len(all_samples)

random.seed(SUBSAMPLE_SEED)
random.shuffle(all_samples)

start_idx = (ITER_ID - 1) * SAMPLES_PER_ITER
end_idx = ITER_ID * SAMPLES_PER_ITER

required_samples = N_ITERATIONS * SAMPLES_PER_ITER
if end_idx > total_available:
    raise ValueError(
        f"STRATIFICATION ERROR: Your atlas only contains {total_available} unique samples. "
        f"To run {N_ITERATIONS} iterations of {SAMPLES_PER_ITER} samples with ABSOLUTELY ZERO OVERLAP, "
        f"you need at least {required_samples} unique {BATCH_KEY} values in your master dataset."
    )

selected_samples = all_samples[start_idx:end_idx]
logger.info(f"Iteration {ITER_ID} non-overlapping window: Indices [{start_idx}:{end_idx}]")
logger.info(f"Selected Samples for this task: {selected_samples}")

keep_mask = adata_full.obs[BATCH_KEY].isin(selected_samples)
adata = adata_full[keep_mask].copy()
adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].cat.remove_unused_categories()
logger.info(f"Subsampled shape: {adata.shape}")

del adata_full
gc.collect()

logger.info(f"Recalculating HVGs for {SAMPLES_PER_ITER}-sample subset...")
sc.pp.highly_variable_genes(
    adata,
    layer=COUNTS_LAYER,
    flavor=CFG["subsample"]["hvg_flavor"],
    n_top_genes=CFG["subsample"]["n_top_genes"],
    batch_key=BATCH_KEY,
    subset=True,
)
logger.info(f"Integration input shape: {adata.shape}")

adata_meta = ad.AnnData(obs=adata.obs.copy())

logger.info("Exporting MTX for R...")
scipy.io.mmwrite(
    checkpoint_path(EXPORT_CFG["mtx_counts"]),
    adata.layers[COUNTS_LAYER].T,
)
adata.obs.to_csv(checkpoint_path(EXPORT_CFG["mtx_metadata"]))
pd.DataFrame(adata.var_names).to_csv(
    checkpoint_path(EXPORT_CFG["mtx_genes"]), index=False, header=False
)
logger.info("MTX export complete")

METHOD = "pca"
logger.info(f"\n{'─'*60}\n[{METHOD}] Unintegrated PCA baseline")
log_ram(f"{METHOD} start")
t0 = time.time()

adata_pca = adata.copy()
sc.pp.scale(adata_pca, max_value=SCALE_MAX)
sc.tl.pca(adata_pca, n_comps=N_PCS, svd_solver=PCA_SOLVER, random_state=RANDOM_STATE)
logger.info(f"  [{METHOD}] PCA shape: {adata_pca.obsm['X_pca'].shape}")
sc.pp.neighbors(adata_pca, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_STATE)

save_checkpoint(METHOD, adata_pca.obsm["X_pca"], adata_pca, n_components=N_PCS)

del adata_pca
gc.collect()
log_ram(f"{METHOD} done")
logger.info(f"  [{METHOD}] elapsed: {(time.time()-t0)/60:.2f} min")

try:
    import harmonypy as hm
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "harmonypy", "-q"])
    import harmonypy as hm

HARMONY_CFG = INT_CFG["harmony"]
METHOD = "harmony"
logger.info(f"\n{'─'*60}\n[{METHOD}] Harmony integration")
log_ram(f"{METHOD} start")
t0 = time.time()

adata_h = adata.copy()
sc.pp.scale(adata_h, max_value=SCALE_MAX)
sc.tl.pca(adata_h, n_comps=N_PCS, svd_solver=PCA_SOLVER, random_state=RANDOM_STATE)

ho = hm.run_harmony(
    adata_h.obsm["X_pca"],
    adata_h.obs,
    vars_use=[BATCH_KEY],
    max_iter_harmony=HARMONY_CFG["max_iter_harmony"],
    max_iter_kmeans=HARMONY_CFG["max_iter_kmeans"],
    random_state=RANDOM_STATE,
    verbose=True,
)

harmony_emb = ho.Z_corr
logger.info(f"  [{METHOD}] Z_corr shape (no transpose): {harmony_emb.shape}")

adata_h.obsm["X_harmony"] = harmony_emb
sc.pp.neighbors(adata_h, use_rep="X_harmony", n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)

save_checkpoint(METHOD, harmony_emb, adata_h, n_components=N_PCS)

del adata_h, ho, harmony_emb
gc.collect()
log_ram(f"{METHOD} done")
logger.info(f"  [{METHOD}] elapsed: {(time.time()-t0)/60:.2f} min")

try:
    import scanorama
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scanorama", "-q"])
    import scanorama

SCANORAMA_DIMRED = INT_CFG["scanorama"]["dimred"]
METHOD = "scanorama"
logger.info(f"\n{'─'*60}\n[{METHOD}] Scanorama integration")
log_ram(f"{METHOD} start")
t0 = time.time()

adata_scn = adata.copy()
sample_ids = adata_scn.obs[BATCH_KEY].cat.categories.tolist()
logger.info(f"  [{METHOD}] Splitting into {len(sample_ids)} per-sample matrices")

adatas_list = [adata_scn[adata_scn.obs[BATCH_KEY] == sid].copy() for sid in sample_ids]
scanorama.integrate_scanpy(adatas_list, dimred=SCANORAMA_DIMRED, verbose=True)

emb_dict = {
    bc: sub.obsm["X_scanorama"][i]
    for sub in adatas_list
    for i, bc in enumerate(sub.obs_names)
}
scanorama_emb = np.vstack([emb_dict[bc] for bc in adata_scn.obs_names])
logger.info(f"  [{METHOD}] Concatenated embedding shape: {scanorama_emb.shape}")

adata_scn.obsm["X_scanorama"] = scanorama_emb
sc.pp.neighbors(adata_scn, use_rep="X_scanorama", n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)

save_checkpoint(METHOD, scanorama_emb, adata_scn, n_components=N_PCS)

del adata_scn, adatas_list, emb_dict, scanorama_emb
gc.collect()
log_ram(f"{METHOD} done")
logger.info(f"  [{METHOD}] elapsed: {(time.time()-t0)/60:.2f} min")

METHOD = "combat"
logger.info(f"\n{'─'*60}\n[{METHOD}] ComBat batch correction")
log_ram(f"{METHOD} start")
t0 = time.time()

adata_cb = adata.copy()
logger.info(f"  [{METHOD}] Densifying .X for ComBat...")
if sp.issparse(adata_cb.X):
    adata_cb.X = adata_cb.X.toarray().astype(np.float32)
else:
    adata_cb.X = np.asarray(adata_cb.X, dtype=np.float32)
log_ram(f"{METHOD} after densification")

sc.pp.combat(adata_cb, key=BATCH_KEY)
sc.pp.scale(adata_cb, max_value=SCALE_MAX)
sc.tl.pca(adata_cb, n_comps=N_PCS, svd_solver=PCA_SOLVER, random_state=RANDOM_STATE)
sc.pp.neighbors(adata_cb, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_STATE)

save_checkpoint(METHOD, adata_cb.obsm["X_pca"], adata_cb, n_components=N_PCS)

del adata_cb
gc.collect()
log_ram(f"{METHOD} done")
logger.info(f"  [{METHOD}] elapsed: {(time.time()-t0)/60:.2f} min")

try:
    import bbknn
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "bbknn", "-q"])
    import bbknn

METHOD = "bbknn"
logger.info(f"\n{'─'*60}\n[{METHOD}] BBKNN integration")
log_ram(f"{METHOD} start")
t0 = time.time()

adata_bbknn = adata.copy()
sc.pp.scale(adata_bbknn, max_value=SCALE_MAX)
sc.tl.pca(adata_bbknn, n_comps=N_PCS, svd_solver=PCA_SOLVER, random_state=RANDOM_STATE)

bbknn.bbknn(adata_bbknn, batch_key=BATCH_KEY, n_pcs=N_PCS)
logger.info(f"  [{METHOD}] BBKNN graph computed")

save_checkpoint(METHOD, adata_bbknn.obsm["X_pca"], adata_bbknn, n_components=N_PCS)

del adata_bbknn
gc.collect()
log_ram(f"{METHOD} done")
logger.info(f"  [{METHOD}] elapsed: {(time.time()-t0)/60:.2f} min")

import torch

try:
    import scvi
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scvi-tools", "-q"])
    import scvi

SCVI_CFG = INT_CFG["scvi"]
METHOD = "scvi"
logger.info(f"\n{'─'*60}\n[{METHOD}] scVI VAE")
logger.info(f"  [{METHOD}] CUDA: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    logger.info(f"  [{METHOD}] GPU: {torch.cuda.get_device_name(0)}")
log_ram(f"{METHOD} start")
t0 = time.time()

adata_scvi = adata.copy()
adata_scvi.X = adata_scvi.layers[COUNTS_LAYER].copy()

scvi.model.SCVI.setup_anndata(adata_scvi, layer=None, batch_key=BATCH_KEY)

model_scvi = scvi.model.SCVI(
    adata_scvi,
    n_latent=SCVI_CFG["n_latent"],
    n_layers=SCVI_CFG["n_layers"],
    n_hidden=SCVI_CFG["n_hidden"],
    dropout_rate=SCVI_CFG["dropout_rate"],
    dispersion=SCVI_CFG["dispersion"],
    gene_likelihood=SCVI_CFG["gene_likelihood"],
    latent_distribution=SCVI_CFG["latent_distribution"],
)

USE_GPU = torch.cuda.is_available()
model_scvi.train(
    max_epochs=SCVI_CFG["max_epochs"],
    batch_size=SCVI_CFG["batch_size"],
    early_stopping=True,
    early_stopping_patience=SCVI_CFG["early_stopping_patience"],
    accelerator="gpu" if USE_GPU else "cpu",
    devices=1 if USE_GPU else "auto",
)

scvi_latent = model_scvi.get_latent_representation(adata_scvi)
logger.info(f"  [{METHOD}] Latent shape: {scvi_latent.shape}")

scvi_model_path = os.path.join(CHECKPOINT_DIR, SCVI_CFG["model_dir"])
model_scvi.save(scvi_model_path, overwrite=True)
logger.info(f"  [{METHOD}] Model saved to: {scvi_model_path}")

adata_scvi.obsm["X_scVI"] = scvi_latent
sc.pp.neighbors(adata_scvi, use_rep="X_scVI", n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)

save_checkpoint(METHOD, scvi_latent, adata_scvi, n_components=SCVI_CFG["n_latent"])

del adata_scvi, scvi_latent
gc.collect()
if USE_GPU:
    torch.cuda.empty_cache()
log_ram(f"{METHOD} done")
logger.info(f"  [{METHOD}] elapsed: {(time.time()-t0)/60:.2f} min")

SCANVI_CFG = INT_CFG["scanvi"]
METHOD = "scanvi"
logger.info(f"\n{'─'*60}\n[{METHOD}] scANVI semi-supervised VAE (warm-start from scVI)")
log_ram(f"{METHOD} start")
t0 = time.time()

adata_scanvi = adata.copy()
adata_scanvi.X = adata_scanvi.layers[COUNTS_LAYER].copy()
scvi.model.SCVI.setup_anndata(adata_scanvi, layer=None, batch_key=BATCH_KEY)

model_scvi_parent = scvi.model.SCVI.load(scvi_model_path, adata=adata_scanvi)

model_scanvi = scvi.model.SCANVI.from_scvi_model(
    model_scvi_parent,
    unlabeled_category=SCANVI_CFG["unlabeled_category"],
    labels_key=LABEL_KEY,
)

model_scanvi.train(
    max_epochs=SCANVI_CFG["max_epochs"],
    batch_size=SCANVI_CFG["batch_size"],
    early_stopping=True,
    early_stopping_patience=SCANVI_CFG["early_stopping_patience"],
    accelerator="gpu" if USE_GPU else "cpu",
    devices=1 if USE_GPU else "auto",
)

scanvi_latent = model_scanvi.get_latent_representation(adata_scanvi)
logger.info(f"  [{METHOD}] Latent shape: {scanvi_latent.shape}")

adata_scanvi.obsm["X_scANVI"] = scanvi_latent
sc.pp.neighbors(adata_scanvi, use_rep="X_scANVI", n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)

save_checkpoint(METHOD, scanvi_latent, adata_scanvi, n_components=SCVI_CFG["n_latent"])

del adata_scanvi, model_scvi_parent, model_scanvi, scanvi_latent
gc.collect()
if USE_GPU:
    torch.cuda.empty_cache()
log_ram(f"{METHOD} done")
logger.info(f"  [{METHOD}] elapsed: {(time.time()-t0)/60:.2f} min")

SCPOLI_CFG = INT_CFG["scpoli"]
METHOD = "scpoli"
logger.info(f"\n{'─'*60}\n[{METHOD}] scPoli conditional VAE")
log_ram(f"{METHOD} start")
t0 = time.time()

if not hasattr(ad, "read"):
    ad.read = ad.read_h5ad

try:
    from scarches.models import scPoli
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scarches", "-q"])
    from scarches.models import scPoli

try:
    adata_scpoli = adata.copy()

    hidden_size = int(np.sqrt(adata_scpoli.n_obs))
    logger.info(
        f"  [{METHOD}] n_obs={adata_scpoli.n_obs}, hidden_size={hidden_size}, "
        f"latent_dim={SCPOLI_CFG['latent_dim']}"
    )

    model = scPoli(
        adata=adata_scpoli,
        condition_keys=BATCH_KEY,
        cell_type_keys=LABEL_KEY,
        hidden_layer_sizes=[hidden_size],
        latent_dim=SCPOLI_CFG["latent_dim"],
        embedding_dims=SCPOLI_CFG["embedding_dims"],
        recon_loss=SCPOLI_CFG["recon_loss"],
    )

    model.train(
        n_epochs=SCPOLI_CFG["n_epochs"],
        pretraining_epochs=SCPOLI_CFG["pretraining_epochs"],
        eta=SCPOLI_CFG["eta"],
        alpha_epoch_anneal=SCPOLI_CFG["alpha_epoch_anneal"],
        early_stopping_kwargs={
            "early_stopping_metric": SCPOLI_CFG["early_stopping_metric"],
            "mode": SCPOLI_CFG["early_stopping_mode"],
            "threshold": SCPOLI_CFG["early_stopping_threshold"],
            "patience": SCPOLI_CFG["early_stopping_patience"],
            "reduce_lr": SCPOLI_CFG["reduce_lr"],
            "lr_patience": SCPOLI_CFG["lr_patience"],
            "lr_factor": SCPOLI_CFG["lr_factor"],
        },
        use_early_stopping=True,
    )

    scpoli_emb = model.get_latent(adata_scpoli)
    logger.info(f"  [{METHOD}] Raw embedding shape: {scpoli_emb.shape}")

    n_cells = adata_scpoli.n_obs
    if scpoli_emb.shape[0] == n_cells:
        final_emb = scpoli_emb
    elif scpoli_emb.shape[1] == n_cells:
        logger.warning(f"  [{METHOD}] Transposing embedding from {scpoli_emb.shape}")
        final_emb = scpoli_emb.T
    else:
        raise ValueError(
            f"[{METHOD}] Shape {scpoli_emb.shape} matches neither dim for n_cells={n_cells}"
        )

    adata_scpoli.obsm["X_scpoli"] = final_emb
    sc.pp.neighbors(adata_scpoli, use_rep="X_scpoli", n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE)

    save_checkpoint(METHOD, final_emb, adata_scpoli, n_components=final_emb.shape[1])

    elapsed = time.time() - t0
    logger.info(f"  [{METHOD}] elapsed: {elapsed/60:.2f} min")

except Exception as e:
    logger.error(f"[{METHOD}] FAILED: {e}", exc_info=True)

finally:
    for var in ["adata_scpoli", "model", "scpoli_emb", "final_emb"]:
        if var in locals():
            del locals()[var]
    gc.collect()
    log_ram(f"{METHOD} done")

logger.info("=" * 60)
logger.info(f"All integration methods complete")
logger.info(f"Checkpoints: {CHECKPOINT_DIR}")

import glob
h5ad_files = sorted(glob.glob(os.path.join(CHECKPOINT_DIR, "checkpoint_*.h5ad")))
logger.info(f"  {len(h5ad_files)} checkpoint(s) written:")
for f in h5ad_files:
    size_mb = os.path.getsize(f) / 1024**2
    logger.info(f"    {os.path.basename(f)}  ({size_mb:.1f} MB)")
