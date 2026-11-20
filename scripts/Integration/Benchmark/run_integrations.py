#!/usr/bin/env python3
"""Run Python integration methods for one subsample iteration.

Loads the master atlas, takes a non-overlapping sample_id window, recomputes HVGs,
exports MTX for R, then runs label-free methods plus scVI / scANVI / scPoli (Level_2).
Writes checkpoint_<method>.h5ad under checkpoints/subsample_iter_N/.
"""
from __future__ import annotations

import gc
import glob
import logging
import os
import random
import shutil
import subprocess
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import psutil
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCH_DIR = os.path.dirname(SCRIPT_DIR)
with open(os.path.join(BENCH_DIR, "config.yaml"), encoding="utf-8") as _fh:
    _CFG = yaml.safe_load(_fh)

MASTER_H5AD = str(_CFG["paths"]["master_h5ad"])
CHECKPOINTS = os.path.join(BENCH_DIR, "checkpoints")
N_ITERS = int(_CFG["benchmark"]["n_iters"])
SUBSAMPLE_SEED = int(_CFG["benchmark"]["subsample_seed"])
RANDOM_STATE = int(_CFG["benchmark"]["random_state"])
N_TOP_HVG = int(_CFG["benchmark"]["n_top_hvg"])
BATCH_KEY = str(_CFG["benchmark"]["batch_key"])
LABEL_KEY = str(_CFG["benchmark"]["label_key"])
_METHOD_FLAGS = {str(k): bool(v) for k, v in _CFG["methods"].items()}
_R_METHODS = {"seurat_v3", "seurat_v4", "css_pearson", "css_spearman"}


def method_enabled(name: str) -> bool:
    return bool(_METHOD_FLAGS.get(name, False))


def enabled_methods() -> list[str]:
    return [m for m, on in _METHOD_FLAGS.items() if on]


def any_r_method_enabled() -> bool:
    return any(method_enabled(m) for m in _R_METHODS)


def iter_id() -> int:
    return int(os.environ.get("BENCHMARK_ITER_ID", "1"))


def iter_dir(i: int | None = None) -> str:
    i = iter_id() if i is None else i
    return os.path.join(CHECKPOINTS, f"subsample_iter_{i}")

if not hasattr(np, "infty"):
    np.infty = np.inf  # type: ignore[attr-defined]
if not hasattr(ad, "read"):
    ad.read = ad.read_h5ad

def _patch_torch_numpy2_early() -> None:
    import torch

    def _wrap(fn, flag_name: str):
        if getattr(fn, flag_name, False):
            return fn
        _orig = fn

        def _compat(data, *args, **kwargs):  # noqa: ANN001
            if isinstance(data, np.ndarray) and "dtype" not in kwargs:
                try:
                    return _orig(data, *args, **kwargs)
                except RuntimeError as exc:
                    if "Could not infer dtype" not in str(exc):
                        raise
                    if np.issubdtype(data.dtype, np.floating):
                        return _orig(np.asarray(data, dtype=np.float64), *args, **kwargs)
                    if np.issubdtype(data.dtype, np.integer):
                        return _orig(data.tolist(), *args, **kwargs)
                    return _orig(np.asarray(data).tolist(), *args, **kwargs)
            return _orig(data, *args, **kwargs)

        setattr(_compat, flag_name, True)
        return _compat

    if not getattr(torch.tensor, "_atlas_numpy2_patched", False):
        torch.tensor = _wrap(torch.tensor, "_atlas_numpy2_patched")
    if not getattr(torch.as_tensor, "_atlas_numpy2_as_patched", False):
        torch.as_tensor = _wrap(torch.as_tensor, "_atlas_numpy2_as_patched")

_patch_torch_numpy2_early()

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("benchmark_integrations")

def log_ram(stage: str = "") -> None:
    ram_gb = psutil.Process(os.getpid()).memory_info().rss / 1024**3
    logger.info("RAM %s: %.2f GB", stage, ram_gb)

ITER_ID = iter_id()
CHECKPOINT_DIR = iter_dir(ITER_ID)
os.makedirs(CHECKPOINT_DIR, exist_ok=True)
os.makedirs(CHECKPOINTS, exist_ok=True)

def checkpoint_path(filename: str) -> str:
    return os.path.join(CHECKPOINT_DIR, filename)

def save_checkpoint(method_name: str, embedding: np.ndarray, source_adata: ad.AnnData, n_components: int) -> str:
    n_cells = source_adata.n_obs
    if embedding.shape[0] != n_cells:
        raise ValueError(
            f"[{method_name}] embedding shape {embedding.shape} != n_cells={n_cells}"
        )

    if method_name == "pca":
        ckpt = ad.AnnData(X=source_adata.layers["counts"].copy(), obs=adata_meta.obs.copy())
        ckpt.var = source_adata.var.copy()
    else:
        ckpt = ad.AnnData(obs=adata_meta.obs.copy())

    ckpt.obsm["X_emb"] = embedding.copy()
    ckpt.obsp["connectivities"] = source_adata.obsp["connectivities"].copy()
    ckpt.obsp["distances"] = source_adata.obsp["distances"].copy()
    ckpt.uns["neighbors"] = source_adata.uns["neighbors"].copy()
    ckpt.uns["method"] = method_name
    ckpt.uns["n_components"] = n_components
    ckpt.uns["label_key_used"] = LABEL_KEY if method_name in {"scanvi", "scpoli"} else "none"
    ckpt.uns["benchmark_variant"] = "level2_org"

    out_path = checkpoint_path(f"checkpoint_{method_name}.h5ad")
    ckpt.write_h5ad(out_path, compression="gzip")
    logger.info("[%s] wrote %s", method_name, out_path)
    del ckpt
    return out_path

def has_valid_checkpoint(method_name: str) -> bool:
    path = checkpoint_path(f"checkpoint_{method_name}.h5ad")
    if not os.path.isfile(path):
        return False
    try:
        ckpt = ad.read_h5ad(path, backed="r")
        ok = ckpt.n_obs > 0 and "X_emb" in ckpt.obsm and ckpt.obsm["X_emb"].shape[0] == ckpt.n_obs
        if method_name == "pca":
            ok = ok and (ckpt.X is not None)
        ckpt.file.close()
        if not ok:
            raise ValueError("missing X_emb / X or shape mismatch")
        logger.info("SKIP [%s] valid checkpoint already present", method_name)
        return True
    except Exception as exc:
        broken = path + ".broken"
        logger.warning("Broken checkpoint %s (%s) → quarantine as %s", path, exc, broken)
        try:
            if os.path.isfile(broken):
                os.remove(broken)
            os.rename(path, broken)
        except OSError as rename_exc:
            logger.warning("Could not quarantine %s: %s", path, rename_exc)
        return False

logger.info("=" * 72)
logger.info("Level_2 integrations | iter=%s/%s | label=%s", ITER_ID, N_ITERS, LABEL_KEY)
logger.info("Enabled methods: %s", ", ".join(enabled_methods()) or "(none)")
logger.info("Master: %s", MASTER_H5AD)
logger.info("Output: %s", CHECKPOINT_DIR)
logger.info("=" * 72)

if not enabled_methods():
    raise RuntimeError("No methods enabled in config.yaml")

if ITER_ID < 1 or ITER_ID > N_ITERS:
    raise ValueError(f"ITER_ID must be in 1..{N_ITERS}, got {ITER_ID}")

adata_full = sc.read_h5ad(MASTER_H5AD)
if LABEL_KEY not in adata_full.obs.columns:
    raise KeyError(f"{LABEL_KEY} missing from master obs")
if BATCH_KEY not in adata_full.obs.columns:
    raise KeyError(f"{BATCH_KEY} missing from master obs")
if "counts" not in adata_full.layers:
    raise KeyError("layers['counts'] required on master")

all_samples = sorted(adata_full.obs[BATCH_KEY].unique().tolist())
total_available = len(all_samples)
random.seed(SUBSAMPLE_SEED)
random.shuffle(all_samples)

if total_available < N_ITERS:
    raise ValueError(f"Need >= {N_ITERS} sample_ids, found {total_available}")

base = total_available // N_ITERS
rem = total_available % N_ITERS
sizes = [base + (1 if i < rem else 0) for i in range(N_ITERS)]
starts = [0]
for s in sizes[:-1]:
    starts.append(starts[-1] + s)
start_idx = starts[ITER_ID - 1]
end_idx = start_idx + sizes[ITER_ID - 1]
selected_samples = all_samples[start_idx:end_idx]

logger.info(
    "Iteration %s/%s non-overlapping window [%s:%s] (n=%s of %s)",
    ITER_ID, N_ITERS, start_idx, end_idx, len(selected_samples), total_available,
)
logger.info("Selected samples: %s", selected_samples)

keep_mask = adata_full.obs[BATCH_KEY].isin(selected_samples)
adata = adata_full[keep_mask].copy()
adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].cat.remove_unused_categories()
adata.obs[LABEL_KEY] = adata.obs[LABEL_KEY].astype(str).astype("category")
logger.info("Subsampled shape: %s", adata.shape)

del adata_full
gc.collect()

logger.info("Recalculating HVGs (n_top=%s, batch=%s)...", N_TOP_HVG, BATCH_KEY)
sc.pp.highly_variable_genes(
    adata,
    layer="counts",
    flavor="seurat_v3",
    n_top_genes=N_TOP_HVG,
    batch_key=BATCH_KEY,
    subset=True,
)
logger.info("Integration input shape: %s", adata.shape)

adata_meta = ad.AnnData(obs=adata.obs.copy())

# Shared reference for metrics (obs + counts) even if pca is disabled for scoring
ref = ad.AnnData(X=adata.layers["counts"].copy(), obs=adata.obs.copy(), var=adata.var.copy())
ref.layers["counts"] = ref.X.copy()
ref.write_h5ad(checkpoint_path("reference.h5ad"), compression="gzip")
logger.info("Wrote reference.h5ad")

if any_r_method_enabled():
    logger.info("Exporting MTX for R methods...")
    scipy.io.mmwrite(checkpoint_path("heoca_counts.mtx"), adata.layers["counts"].T)
    adata.obs.to_csv(checkpoint_path("heoca_metadata.csv"))
    pd.DataFrame(adata.var_names).to_csv(checkpoint_path("heoca_genes.csv"), index=False, header=False)
    logger.info("MTX export complete")
else:
    logger.info("SKIP MTX export — no R methods enabled")

METHOD = "pca"
if method_enabled(METHOD) and not has_valid_checkpoint(METHOD):
    logger.info("[%s] unintegrated PCA baseline", METHOD)
    t0 = time.time()
    adata_pca = adata.copy()
    sc.pp.scale(adata_pca, max_value=10)
    sc.tl.pca(adata_pca, n_comps=50, svd_solver="arpack", random_state=RANDOM_STATE)
    sc.pp.neighbors(adata_pca, n_neighbors=15, n_pcs=50, random_state=RANDOM_STATE)
    save_checkpoint(METHOD, adata_pca.obsm["X_pca"], adata_pca, n_components=50)
    del adata_pca
    gc.collect()
    logger.info("[%s] done in %.2f min", METHOD, (time.time() - t0) / 60.0)
elif not method_enabled(METHOD):
    logger.info("SKIP [%s] disabled in config.yaml", METHOD)

METHOD = "harmony"
if method_enabled(METHOD) and not has_valid_checkpoint(METHOD):
    try:
        import harmonypy as hm
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "harmonypy", "-q"])
        import harmonypy as hm
    logger.info("[%s] Harmony", METHOD)
    t0 = time.time()
    adata_h = adata.copy()
    sc.pp.scale(adata_h, max_value=10)
    sc.tl.pca(adata_h, n_comps=50, svd_solver="arpack", random_state=RANDOM_STATE)
    ho = hm.run_harmony(
        adata_h.obsm["X_pca"],
        adata_h.obs,
        vars_use=[BATCH_KEY],
        max_iter_harmony=30,
        max_iter_kmeans=20,
        random_state=RANDOM_STATE,
        verbose=True,
    )
    harmony_emb = ho.Z_corr
    adata_h.obsm["X_harmony"] = harmony_emb
    sc.pp.neighbors(adata_h, use_rep="X_harmony", n_neighbors=15, random_state=RANDOM_STATE)
    save_checkpoint(METHOD, harmony_emb, adata_h, n_components=50)
    del adata_h, ho, harmony_emb
    gc.collect()
    logger.info("[%s] done in %.2f min", METHOD, (time.time() - t0) / 60.0)
elif not method_enabled(METHOD):
    logger.info("SKIP [%s] disabled in config.yaml", METHOD)

METHOD = "scanorama"
if method_enabled(METHOD) and not has_valid_checkpoint(METHOD):
    try:
        import scanorama
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scanorama", "-q"])
        import scanorama
    logger.info("[%s] Scanorama", METHOD)
    t0 = time.time()
    adata_scn = adata.copy()
    sample_ids = adata_scn.obs[BATCH_KEY].cat.categories.tolist()
    adatas_list = [adata_scn[adata_scn.obs[BATCH_KEY] == sid].copy() for sid in sample_ids]
    scanorama.integrate_scanpy(adatas_list, dimred=50, verbose=True)
    emb_dict = {
        bc: sub.obsm["X_scanorama"][i]
        for sub in adatas_list
        for i, bc in enumerate(sub.obs_names)
    }
    scanorama_emb = np.vstack([emb_dict[bc] for bc in adata_scn.obs_names])
    adata_scn.obsm["X_scanorama"] = scanorama_emb
    sc.pp.neighbors(adata_scn, use_rep="X_scanorama", n_neighbors=15, random_state=RANDOM_STATE)
    save_checkpoint(METHOD, scanorama_emb, adata_scn, n_components=50)
    del adata_scn, adatas_list, emb_dict, scanorama_emb
    gc.collect()
    logger.info("[%s] done in %.2f min", METHOD, (time.time() - t0) / 60.0)
elif not method_enabled(METHOD):
    logger.info("SKIP [%s] disabled in config.yaml", METHOD)

METHOD = "combat"
if method_enabled(METHOD) and not has_valid_checkpoint(METHOD):
    logger.info("[%s] ComBat", METHOD)
    t0 = time.time()
    adata_cb = adata.copy()
    if sp.issparse(adata_cb.X):
        adata_cb.X = adata_cb.X.toarray().astype(np.float32)
    else:
        adata_cb.X = np.asarray(adata_cb.X, dtype=np.float32)
    sc.pp.combat(adata_cb, key=BATCH_KEY)
    sc.pp.scale(adata_cb, max_value=10)
    sc.tl.pca(adata_cb, n_comps=50, svd_solver="arpack", random_state=RANDOM_STATE)
    sc.pp.neighbors(adata_cb, n_neighbors=15, n_pcs=50, random_state=RANDOM_STATE)
    save_checkpoint(METHOD, adata_cb.obsm["X_pca"], adata_cb, n_components=50)
    del adata_cb
    gc.collect()
    logger.info("[%s] done in %.2f min", METHOD, (time.time() - t0) / 60.0)
elif not method_enabled(METHOD):
    logger.info("SKIP [%s] disabled in config.yaml", METHOD)

METHOD = "bbknn"
if method_enabled(METHOD) and not has_valid_checkpoint(METHOD):
    try:
        import bbknn
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "bbknn", "-q"])
        import bbknn
    logger.info("[%s] BBKNN", METHOD)
    t0 = time.time()
    adata_bbknn = adata.copy()
    sc.pp.scale(adata_bbknn, max_value=10)
    sc.tl.pca(adata_bbknn, n_comps=50, svd_solver="arpack", random_state=RANDOM_STATE)
    bbknn.bbknn(adata_bbknn, batch_key=BATCH_KEY, n_pcs=50)
    save_checkpoint(METHOD, adata_bbknn.obsm["X_pca"], adata_bbknn, n_components=50)
    del adata_bbknn
    gc.collect()
    logger.info("[%s] done in %.2f min", METHOD, (time.time() - t0) / 60.0)
elif not method_enabled(METHOD):
    logger.info("SKIP [%s] disabled in config.yaml", METHOD)

need_scvi = method_enabled("scvi") or method_enabled("scanvi")
need_scpoli = method_enabled("scpoli")

if need_scvi or need_scpoli:
    import torch

    USE_GPU = torch.cuda.is_available()
    logger.info("CUDA available: %s", USE_GPU)
else:
    USE_GPU = False

if need_scvi:
    try:
        import scvi
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scvi-tools", "-q"])
        import scvi

    scvi_model_path = os.path.join(CHECKPOINT_DIR, "scvi_model")
    train_scvi = method_enabled("scvi") and not has_valid_checkpoint("scvi")
    # scANVI warm-start needs a saved scVI model even if scvi scoring is off
    if method_enabled("scanvi") and not has_valid_checkpoint("scanvi") and not os.path.isdir(scvi_model_path):
        train_scvi = True

    if train_scvi:
        logger.info("[scvi] training%s", "" if method_enabled("scvi") else " (dependency for scANVI)")
        t0 = time.time()
        adata_scvi = adata.copy()
        x_counts = adata_scvi.layers["counts"].copy()
        adata_scvi.X = (
            x_counts.tocsr().astype(np.float32) if sp.issparse(x_counts) else np.asarray(x_counts, dtype=np.float32)
        )
        scvi.model.SCVI.setup_anndata(adata_scvi, layer=None, batch_key=BATCH_KEY)
        model_scvi = scvi.model.SCVI(
            adata_scvi,
            n_latent=30,
            n_layers=2,
            n_hidden=256,
            dropout_rate=0.1,
            dispersion="gene",
            gene_likelihood="nb",
            latent_distribution="normal",
        )
        model_scvi.train(
            max_epochs=200,
            batch_size=2048,
            early_stopping=True,
            early_stopping_patience=15,
            accelerator="gpu" if USE_GPU else "cpu",
            devices=1 if USE_GPU else "auto",
        )
        scvi_latent = model_scvi.get_latent_representation(adata_scvi)
        model_scvi.save(scvi_model_path, overwrite=True)
        if method_enabled("scvi"):
            adata_scvi.obsm["X_scVI"] = scvi_latent
            sc.pp.neighbors(adata_scvi, use_rep="X_scVI", n_neighbors=15, random_state=RANDOM_STATE)
            save_checkpoint("scvi", scvi_latent, adata_scvi, n_components=30)
        del adata_scvi, scvi_latent, model_scvi
        gc.collect()
        if USE_GPU:
            torch.cuda.empty_cache()
        logger.info("[scvi] done in %.2f min", (time.time() - t0) / 60.0)
    elif not method_enabled("scvi"):
        logger.info("SKIP [scvi] disabled in config.yaml")

    METHOD = "scanvi"
    if method_enabled(METHOD) and not has_valid_checkpoint(METHOD):
        if not os.path.isdir(scvi_model_path):
            raise FileNotFoundError(f"scVI model required for scANVI warm-start: {scvi_model_path}")
        logger.info("[%s] scANVI warm-start from scVI | labels=%s", METHOD, LABEL_KEY)
        t0 = time.time()
        adata_scanvi = adata.copy()
        x_counts = adata_scanvi.layers["counts"].copy()
        adata_scanvi.X = (
            x_counts.tocsr().astype(np.float32) if sp.issparse(x_counts) else np.asarray(x_counts, dtype=np.float32)
        )
        scvi.model.SCVI.setup_anndata(adata_scanvi, layer=None, batch_key=BATCH_KEY)
        model_scvi_parent = scvi.model.SCVI.load(scvi_model_path, adata=adata_scanvi)
        model_scanvi = scvi.model.SCANVI.from_scvi_model(
            model_scvi_parent,
            unlabeled_category="Unknown",
            labels_key=LABEL_KEY,
        )
        model_scanvi.train(
            max_epochs=200,
            batch_size=2048,
            early_stopping=True,
            early_stopping_patience=15,
            accelerator="gpu" if USE_GPU else "cpu",
            devices=1 if USE_GPU else "auto",
        )
        scanvi_latent = model_scanvi.get_latent_representation(adata_scanvi)
        adata_scanvi.obsm["X_scANVI"] = scanvi_latent
        sc.pp.neighbors(adata_scanvi, use_rep="X_scANVI", n_neighbors=15, random_state=RANDOM_STATE)
        save_checkpoint(METHOD, scanvi_latent, adata_scanvi, n_components=30)
        scanvi_model_dir = checkpoint_path("scanvi_model")
        try:
            if os.path.isdir(scanvi_model_dir):
                shutil.rmtree(scanvi_model_dir)
            model_scanvi.save(scanvi_model_dir, overwrite=True)
        except Exception as exc:
            logger.warning("[%s] model.save failed (non-fatal): %s", METHOD, exc)
        del adata_scanvi, model_scvi_parent, model_scanvi, scanvi_latent
        gc.collect()
        if USE_GPU:
            torch.cuda.empty_cache()
        logger.info("[%s] done in %.2f min", METHOD, (time.time() - t0) / 60.0)
    elif not method_enabled(METHOD):
        logger.info("SKIP [%s] disabled in config.yaml", METHOD)
else:
    logger.info("SKIP [scvi/scanvi] disabled in config.yaml")

METHOD = "scpoli"
if method_enabled(METHOD) and not has_valid_checkpoint(METHOD):
    import torch

    logger.info("[%s] scPoli | cell_type_keys=%s", METHOD, LABEL_KEY)
    t0 = time.time()
    try:
        from scarches.models import scPoli
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "scarches", "-q"])
        from scarches.models import scPoli

    def _patch_scpoli_get_latent_sparse_indexing():
        import torch as _torch
        from scipy import sparse

        if getattr(scPoli.get_latent, "_atlas_sparse_index_patched", False):
            return

        def get_latent(self, adata_in, mean=False):  # noqa: ANN001
            device = next(self.model.parameters()).device
            x = adata_in.X
            c = {k: adata_in.obs[k].values for k in self.condition_keys_}
            label_tensor = []
            for cond in c.keys():
                query_conditions = c[cond]
                if not set(query_conditions).issubset(self.conditions_[cond]):
                    raise ValueError("Incorrect conditions")
                labels = np.zeros(query_conditions.shape[0])
                for condition, label in self.model.condition_encoders[cond].items():
                    labels[query_conditions == condition] = label
                label_tensor.append(labels)
            c = _torch.tensor(np.asarray(label_tensor).T, device=device)
            latents = []
            n_obs = x.shape[0]
            for start in range(0, n_obs, 512):
                batch_np = np.arange(start, min(start + 512, n_obs))
                x_batch = x[batch_np, :]
                if sparse.issparse(x_batch):
                    x_batch = x_batch.toarray()
                x_batch = _torch.tensor(np.asarray(x_batch, dtype=np.float32), device=device)
                batch_t = _torch.as_tensor(batch_np, device=device)
                latent = self.model.get_latent(x_batch, c[batch_t, :], mean)
                latents.append(latent.cpu().detach())
            return np.array(_torch.cat(latents))

        get_latent._atlas_sparse_index_patched = True  # type: ignore[attr-defined]
        scPoli.get_latent = get_latent  # type: ignore[method-assign]

    _patch_scpoli_get_latent_sparse_indexing()

    adata_scpoli = adata.copy()
    x_counts = adata_scpoli.layers["counts"].copy()
    adata_scpoli.X = (
        x_counts.tocsr().astype(np.float32) if sp.issparse(x_counts) else np.asarray(x_counts, dtype=np.float32)
    )
    hidden_size = int(np.sqrt(adata_scpoli.n_obs))
    logger.info("[%s] n_obs=%s hidden=%s", METHOD, adata_scpoli.n_obs, hidden_size)

    model = scPoli(
        adata=adata_scpoli,
        condition_keys=BATCH_KEY,
        cell_type_keys=LABEL_KEY,
        hidden_layer_sizes=[hidden_size],
        latent_dim=30,
        embedding_dims=10,
        recon_loss="nb",
    )
    model.train(
        n_epochs=100,
        pretraining_epochs=80,
        eta=10,
        alpha_epoch_anneal=100,
        early_stopping_kwargs={
            "early_stopping_metric": "val_loss",
            "mode": "min",
            "threshold": 0,
            "patience": 20,
            "reduce_lr": True,
            "lr_patience": 13,
            "lr_factor": 0.1,
        },
        use_early_stopping=True,
    )

    scpoli_model_dir = checkpoint_path("scpoli_model")
    try:
        if os.path.isdir(scpoli_model_dir):
            shutil.rmtree(scpoli_model_dir)
        model.save(scpoli_model_dir, overwrite=True)
    except Exception as save_exc:
        logger.warning("[%s] model.save failed (non-fatal): %s", METHOD, save_exc)

    scpoli_emb = model.get_latent(adata_scpoli)
    n_cells = adata_scpoli.n_obs
    if scpoli_emb.shape[0] == n_cells:
        final_emb = scpoli_emb
    elif scpoli_emb.shape[1] == n_cells:
        final_emb = scpoli_emb.T
    else:
        raise ValueError(f"[{METHOD}] bad emb shape {scpoli_emb.shape} for n={n_cells}")

    adata_scpoli.obsm["X_scpoli"] = final_emb
    sc.pp.neighbors(adata_scpoli, use_rep="X_scpoli", n_neighbors=15, random_state=RANDOM_STATE)
    save_checkpoint(METHOD, final_emb, adata_scpoli, n_components=final_emb.shape[1])
    del adata_scpoli, model, scpoli_emb, final_emb
    gc.collect()
    if USE_GPU:
        torch.cuda.empty_cache()
    logger.info("[%s] done in %.2f min", METHOD, (time.time() - t0) / 60.0)
elif not method_enabled(METHOD):
    logger.info("SKIP [%s] disabled in config.yaml", METHOD)

logger.info("=" * 72)
logger.info("Python integrations complete for iter %s", ITER_ID)
for path in sorted(glob.glob(os.path.join(CHECKPOINT_DIR, "checkpoint_*.h5ad"))):
    logger.info("  %s", os.path.basename(path))
logger.info("=" * 72)
