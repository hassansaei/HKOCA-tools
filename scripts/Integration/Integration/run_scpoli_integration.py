#!/usr/bin/env python3
"""Train scPoli on the annotated organoid atlas and write the integrated object."""
from __future__ import annotations

import gc
import logging
import os
import socket
import sys
import types
import warnings
from pathlib import Path
from typing import Mapping

import anndata as ad
import numpy as np
import scipy.sparse as sp

if not hasattr(np, "infty"):
    np.infty = np.inf  # type: ignore[attr-defined]

_boot_logger = logging.getLogger(__name__)


def _patch_torch_numpy2_dtype_inference() -> None:
    import torch
    from scarches.dataset.scpoli.anndata import MultiConditionAnnotatedDataset

    if not getattr(torch.tensor, "_atlas_numpy2_patched", False):
        _orig_tensor = torch.tensor

        def _tensor_compat(data, *args, **kwargs):  # noqa: ANN001
            if isinstance(data, np.ndarray) and "dtype" not in kwargs:
                try:
                    return _orig_tensor(data, *args, **kwargs)
                except RuntimeError as exc:
                    if "Could not infer dtype" not in str(exc):
                        raise
                    if np.issubdtype(data.dtype, np.floating):
                        return _orig_tensor(
                            np.asarray(data, dtype=np.float64), *args, **kwargs
                        )
                    if np.issubdtype(data.dtype, np.integer):
                        return _orig_tensor(data.tolist(), *args, **kwargs)
                    return _orig_tensor(np.asarray(data).tolist(), *args, **kwargs)
            return _orig_tensor(data, *args, **kwargs)

        _tensor_compat._atlas_numpy2_patched = True  # type: ignore[attr-defined]
        torch.tensor = _tensor_compat  # type: ignore[assignment]

    if getattr(MultiConditionAnnotatedDataset, "_atlas_numpy2_index_patched", False):
        return

    def _coerce_index(index):  # noqa: ANN001
        if isinstance(index, (int, slice)):
            return index
        if isinstance(index, np.integer):
            return int(index)
        if isinstance(index, np.ndarray):
            if index.ndim == 0:
                return int(index)
            return [int(i) for i in index.tolist()]
        if isinstance(index, torch.Tensor):
            if index.dim() == 0:
                return int(index.item())
            return [int(i) for i in index.tolist()]
        return index

    _orig_getitem = MultiConditionAnnotatedDataset.__getitem__

    def _getitem_compat(self, index):  # noqa: ANN001
        return _orig_getitem(self, _coerce_index(index))

    MultiConditionAnnotatedDataset.__getitem__ = _getitem_compat  # type: ignore[method-assign]
    MultiConditionAnnotatedDataset._atlas_numpy2_index_patched = True  # type: ignore[attr-defined]


if not hasattr(ad, "read"):
    ad.read = ad.read_h5ad

if not hasattr(ad, "io"):
    io_module = types.ModuleType("anndata.io")
    _names = [
        "read_csv",
        "read_loom",
        "read_text",
        "read_excel",
        "read_hdf",
        "read_mtx",
        "read_zarr",
        "read_h5ad",
        "read_umi_tools",
        "read_elem",
        "write_elem",
    ]
    _experimental = None
    try:
        import anndata.experimental as _experimental
    except ImportError:
        pass
    for _name in _names:
        if hasattr(ad, _name):
            setattr(io_module, _name, getattr(ad, _name))
        elif _experimental is not None and hasattr(_experimental, _name):
            setattr(io_module, _name, getattr(_experimental, _name))
        else:
            _boot_logger.warning("anndata.%s not found; may cause issues.", _name)
    ad.io = io_module
    sys.modules["anndata.io"] = io_module

if not hasattr(ad, "abc"):
    from anndata._core.sparse_dataset import CSRDataset, CSCDataset

    abc_module = types.ModuleType("anndata.abc")
    abc_module.CSRDataset = CSRDataset
    abc_module.CSCDataset = CSCDataset
    ad.abc = abc_module
    sys.modules["anndata.abc"] = abc_module

sys.modules["anndata"] = ad

import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc
import yaml
from scarches.models import scPoli

_patch_torch_numpy2_dtype_inference()


def _patch_scpoli_get_latent_sparse_indexing() -> None:
    import torch
    from scipy import sparse

    if getattr(scPoli.get_latent, "_atlas_sparse_index_patched", False):
        return

    def get_latent(self, adata, mean=False):  # noqa: ANN001
        device = next(self.model.parameters()).device
        x = adata.X
        c = {k: adata.obs[k].values for k in self.condition_keys_}

        label_tensor = []
        for cond in c.keys():
            query_conditions = c[cond]
            if not set(query_conditions).issubset(self.conditions_[cond]):
                raise ValueError("Incorrect conditions")
            labels = np.zeros(query_conditions.shape[0])
            for condition, label in self.model.condition_encoders[cond].items():
                labels[query_conditions == condition] = label
            label_tensor.append(labels)
        c = torch.tensor(np.asarray(label_tensor).T, device=device)

        latents = []
        n_obs = x.shape[0]
        for start in range(0, n_obs, 512):
            batch_np = np.arange(start, min(start + 512, n_obs))
            x_batch = x[batch_np, :]
            if sparse.issparse(x_batch):
                x_batch = x_batch.toarray()
            x_batch = torch.tensor(
                np.asarray(x_batch, dtype=np.float32), device=device
            )
            batch_t = torch.as_tensor(batch_np, device=device)
            latent = self.model.get_latent(x_batch, c[batch_t, :], mean)
            latents.append(latent.cpu().detach())
        return np.array(torch.cat(latents))

    get_latent._atlas_sparse_index_patched = True  # type: ignore[attr-defined]
    scPoli.get_latent = get_latent  # type: ignore[method-assign]


_patch_scpoli_get_latent_sparse_indexing()
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Color palettes (from Color_palette.md)
# ---------------------------------------------------------------------------
LEVEL_1_COLORS: dict[str, str] = {
    "Nephron": "#0072B2",
    "Ureteric_Epithelium": "#6A51A3",
    "Stromal_cells": "#D55E00",
    "Endothelial_cells": "#009E73",
    "Off_target_cells": "#737373",
    "undifferentiated_pluripotent_cells": "#F0E442",
}

LEVEL_2_COLORS: dict[str, str] = {
    "Nephron_progenitor_cells": "#081D58",
    "Podocytes": "#225EA8",
    "Proximal_tubule": "#1D91C0",
    "LOH": "#008080",
    "Distal_tubule": "#00BFFF",
    "Proliferative_tubule": "#41B6C4",
    "Ureteric_Bud": "#3F007D",
    "Collecting_duct": "#6A51A3",
    "Stromal_progenitors": "#7A0177",
    "Proliferative_stroma": "#C51B8A",
    "Mesangial_cells": "#F768A1",
    "Pericytes": "#E31A1C",
    "Interstitial_Fibroblasts": "#FF7F00",
    "Activated Fibroblasts": "#FDBF6F",
    "Early_endothelium": "#00441B",
    "Late_endothelium": "#238B45",
    "Off_target_cells": "#525252",
    "undifferentiated_pluripotent_cells": "#FFD700",
    "S1": "#E31A1C",
    "S2": "#FF7F00",
    "S3": "#FDBF6F",
    "PEC-like": "#146A8F",
    "Unassigned": "#999999",
    "Endodermal_cells": "#B49C1B",
    "Neural_progenitor_cells": "#F05270",
    "Muscle_progenitor_cells": "#A59430",
    "Glial_cells": "#CFED24",
    "Endothelial_cells": "#238B45",
    "Precursor_podocytes": "#1F4E99",
    "Mature_podocytes": "#438CE3",
}

LEVEL_3_COLORS: dict[str, str] = {
    "Nephron_progenitor_cells": "#081D58",
    "Precursor_podocytes": "#1F4E99",
    "Mature_podocytes": "#438CE3",
    "PT (PEC-like)": "#146A8F",
    "Dev_PT": "#1D91C0",
    "PT": "#29AADB",
    "LOH": "#008080",
    "Dev_DT": "#0086B3",
    "DT": "#33CCFF",
    "Proliferative_tubule": "#41B6C4",
    "Ureteric_Bud": "#3F007D",
    "Collecting_duct": "#544082",
    "Collecting_duct_Principal_cells": "#544082",
    "CD_Principal_cells": "#544082",
    "Collecting_duct_Intercalated_cells": "#8066C2",
    "CD_Intercalated_cells": "#8066C2",
    "Proliferative_collecting_duct_cells": "#9E8BD4",
    "Stromal_progenitors": "#7A0177",
    "Proliferative_stroma": "#C51B8A",
    "Mesangial_cells": "#F768A1",
    "Pericytes": "#E31A1C",
    "Interstitial_Fibroblasts": "#FF7F00",
    "Activated Fibroblasts": "#FDBF6F",
    "Early_endothelium": "#00441B",
    "Late_endothelium": "#238B45",
    "Endodermal_cells": "#B49C1B",
    "Muscle_progenitor_cells": "#A59430",
    "Neural_progenitor_cells": "#F05270",
    "Glial_cells": "#CFED24",
    "Off_target_cells": "#252525",
    "undifferentiated_pluripotent_cells": "#F0E442",
    "Unassigned": "#999999",
    "S1": "#E31A1C",
    "S2": "#FF7F00",
    "S3": "#FDBF6F",
    "PEC-like": "#146A8F",
    "Distal_tubule": "#00BFFF",
    "Endothelial_cells": "#238B45",
}

DIFF_PROTOCOL_COLORS: dict[str, str] = {
    "Morizane_modified": "#1F77B4",
    "Morizane": "#1F77B4",
    "Takasato": "#D62728",
    "Takasato_modified": "#FF7F0E",
    "Freedman": "#2CA02C",
    "Garreta_suspension": "#8C510A",
    "Low": "#FF69B4",
    "Hybrid_Takasato_Morizane": "#4B0082",
    "Hybrid_Takasato_Uchimura": "#F7C948",
    "Uchimura": "#7F7F7F",
    "Vanslambrouck_2022": "#9ECAE1",
}

SC_PROTOCOL_COLORS: dict[str, str] = {
    "10X_5": "#1f77b4",
    "10x_5": "#1f77b4",
    "10x_3_v1": "#ff7f0e",
    "10x_3_v2": "#2ca02c",
    "10x_3_v3": "#d62728",
    "10x_3_v3.1": "#d62728",
    "10x_3_v4": "#8c564b",
    "Drop-seq": "#9467bd",
    "drop seq": "#9467bd",
    "Smart-seq2": "#7f7f7f",
    "PIPseq_V_T10_3": "#e377c2",
}

FALLBACK_CYCLE = [
    "#9C755F",
    "#BAB0AC",
    "#86BCB6",
    "#F1CE63",
    "#C0A0C3",
    "#F4A6B5",
    "#aec7e8",
    "#ffbb78",
    "#98df8a",
    "#ff9896",
    "#c5b0d5",
    "#bcbd22",
    "#17becf",
]


def palette_for_categories(
    categories,
    base: Mapping[str, str],
    fallback: list[str] | None = None,
) -> dict[str, str]:
    fb = list(fallback or FALLBACK_CYCLE)
    out: dict[str, str] = {}
    fi = 0
    for cat in categories:
        key = str(cat)
        if key in base:
            out[key] = base[key]
        else:
            out[key] = fb[fi % len(fb)]
            fi += 1
    return out


HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = HERE / "config.yaml"


def refuse_login() -> None:
    host = socket.gethostname()
    if host == "master" or host.startswith("master"):
        raise SystemExit(
            f"Refusing to run on login node ({host}). Submit on a GPU compute node."
        )


def ensure_new(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")


def ensure_new_dir(path: Path) -> None:
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty dir: {path}")


def save_model_checkpoint(
    model: scPoli,
    path: Path,
    logger: logging.Logger,
    label: str,
) -> None:
    ensure_new_dir(path)
    if path.exists():
        path.rmdir()
    logger.info("Saving %s checkpoint to %s", label, path)
    model.save(str(path), overwrite=False)
    logger.info("%s checkpoint written: %s", label, path)


def rss_gb() -> float:
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1e6
    except Exception:
        pass
    return float("nan")


def log_mem(logger: logging.Logger, msg: str) -> None:
    logger.info("[mem] %s: RSS=%.1f GB", msg, rss_gb())


def load_config() -> dict:
    argv = sys.argv[1:]
    if "--config" in argv:
        config_file = Path(argv[argv.index("--config") + 1]).resolve()
    else:
        config_file = DEFAULT_CONFIG
    with open(config_file, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _materialize_sparse(X):
    if hasattr(X, "to_memory"):
        X = X.to_memory()
    if not sp.issparse(X):
        raise TypeError(f"Expected sparse matrix, got {type(X)}")
    return X.tocsr()


def build_counts_adata_from_annotated(
    annotated_path: Path,
    label_keys: list[str],
    logger: logging.Logger,
) -> ad.AnnData:
    logger.info("Loading annotated object: %s", annotated_path)
    src = ad.read_h5ad(annotated_path, backed="r")
    names = src.obs_names.to_numpy()
    obs = src.obs.copy()
    var = src.var.copy()

    if "counts" not in src.layers:
        src.file.close()
        raise KeyError(f"layers['counts'] missing in {annotated_path}")

    X = _materialize_sparse(src.layers["counts"])
    logger.info("Using layers['counts'] from annotated object")
    src.file.close()
    del src
    gc.collect()
    log_mem(logger, "after counts materialize")

    for key in label_keys:
        if key not in obs.columns:
            logger.warning("Label key missing in annotated object: %s", key)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = names
    adata.var_names_make_unique()
    return adata


def main() -> None:
    refuse_login()
    cfg = load_config()

    out_root = Path(cfg["paths"]["output_dir"])
    obj_dir = out_root / "objects"
    fig_dir = out_root / "figures"
    log_dir = out_root / "logs"
    for d in (
        obj_dir,
        fig_dir,
        log_dir,
        out_root / "tmp" / "mpl",
        out_root / "tmp" / "numba",
        out_root / "tmp" / "xdg",
    ):
        d.mkdir(parents=True, exist_ok=True)

    os.environ.setdefault("TMPDIR", str(out_root / "tmp"))
    os.environ.setdefault("MPLCONFIGDIR", str(out_root / "tmp" / "mpl"))
    os.environ.setdefault("NUMBA_CACHE_DIR", str(out_root / "tmp" / "numba"))
    os.environ.setdefault("XDG_CACHE_HOME", str(out_root / "tmp" / "xdg"))
    os.environ.setdefault("JAX_PLATFORMS", "cpu")

    final_h5ad = obj_dir / cfg["paths"]["output_h5ad"]
    model_dir = out_root / cfg["paths"]["model_dir"]
    pretrain_ckpt_dir = out_root / cfg["paths"]["checkpoint_pretrain_dir"]
    resume_from_model = (
        (model_dir / "model_params.pt").is_file()
        and (model_dir / "attr.pkl").is_file()
        and not final_h5ad.exists()
    )
    ensure_new(final_h5ad)
    if not resume_from_model:
        ensure_new_dir(model_dir)
        ensure_new_dir(pretrain_ckpt_dir)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_dir / "scpoli_integration.log"),
        ],
    )
    logger = logging.getLogger("scpoli")
    logger.info("host=%s", socket.gethostname())
    if resume_from_model:
        logger.info(
            "RESUME: found trained model at %s — skip training, extract latent + write outputs",
            model_dir,
        )

    try:
        import torch

        logger.info(
            "torch %s | cuda_available=%s",
            torch.__version__,
            torch.cuda.is_available(),
        )
        if not torch.cuda.is_available():
            logger.warning(
                "CUDA not available. scPoli will run on CPU and will be very slow."
            )
    except Exception as exc:
        logger.warning("Could not import torch for CUDA check: %s", exc)

    feat_cfg = cfg["feature_selection"]
    pca_cfg = cfg["pca"]
    scpoli_cfg = cfg["scpoli"]
    train_cfg = cfg["training"]
    plot_cfg = cfg["plotting"]

    sc.settings.figdir = str(fig_dir)
    mpl.rcParams["figure.figsize"] = tuple(plot_cfg["figsize"])
    mpl.rcParams["figure.dpi"] = plot_cfg["dpi"]
    mpl.rcParams["savefig.dpi"] = plot_cfg["dpi"]
    mpl.rcParams["savefig.format"] = "png"
    sc.set_figure_params(dpi=plot_cfg["dpi"], facecolor="white", frameon=False)

    adata = build_counts_adata_from_annotated(
        Path(cfg["paths"]["annotated_h5ad"]),
        list(cfg["metadata"]["label_keys"]),
        logger,
    )
    logger.info(
        "Working matrix: %s cells x %s genes",
        f"{adata.n_obs:,}",
        f"{adata.n_vars:,}",
    )
    log_mem(logger, "after AnnData build")

    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    gc.collect()

    for key in cfg["metadata"]["sanitize_keys"]:
        if key not in adata.obs.columns:
            raise KeyError(f"Required obs key missing: {key}")
        adata.obs[key] = adata.obs[key].astype(str).replace(
            ["nan", "NaN", "None", ""], np.nan
        )
        nan_count = int(adata.obs[key].isna().sum())
        if nan_count > 0:
            logger.warning("Found %s missing values in '%s'. Purging...", nan_count, key)
            adata = adata[adata.obs[key].notna()].copy()
            gc.collect()

    for key in ("sample_id", "study", "Level_2", "Level_1"):
        if key in adata.obs.columns:
            adata.obs[key] = adata.obs[key].astype("category")

    logger.info("Sanitization complete. %s cells remaining.", f"{adata.n_obs:,}")

    if sp.issparse(adata.X):
        xmax = float(adata.X.data.max()) if adata.X.nnz else 0.0
    else:
        xmax = float(np.max(adata.X))
    logger.info("Counts matrix max=%.1f (expect raw-like >> 20)", xmax)
    if xmax < 20:
        raise ValueError(
            f"Expression max={xmax:.2f} looks log-normalized. "
            "scPoli NB + seurat_v3 require raw counts."
        )

    logger.info(
        "Selecting HVGs (%s, flavor=%s, batch_key=%s)...",
        feat_cfg["n_top_genes"],
        feat_cfg["hvg_flavor"],
        feat_cfg["batch_key"],
    )
    model_var_csv = model_dir / "var_names.csv"
    if resume_from_model and model_var_csv.is_file():
        import pandas as pd

        model_genes = pd.read_csv(model_var_csv, header=None).iloc[:, 0].astype(str).tolist()
        missing = [g for g in model_genes if g not in adata.var_names]
        if missing:
            raise KeyError(
                f"Resume: {len(missing)} model genes missing from counts "
                f"(e.g. {missing[:5]})"
            )
        adata = adata[:, model_genes].copy()
        logger.info(
            "RESUME: subset to %s genes from %s",
            adata.n_vars,
            model_var_csv,
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=int(feat_cfg["n_top_genes"]),
            flavor=feat_cfg["hvg_flavor"],
            batch_key=feat_cfg["batch_key"],
        )
        adata = adata[:, adata.var["highly_variable"]].copy()
        logger.info("HVG subset: %s genes", adata.n_vars)
    gc.collect()
    log_mem(logger, "after HVG subset")

    logger.info("Baseline PCA representation...")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=float(pca_cfg["normalize_target_sum"]))
    sc.pp.log1p(adata)

    sc.pp.pca(adata, n_comps=int(pca_cfg["n_comps"]))
    sc.pp.neighbors(
        adata,
        use_rep="X_pca",
        key_added="neighbors_pca",
        n_neighbors=int(pca_cfg["n_neighbors"]),
        random_state=int(pca_cfg["random_state"]),
    )
    sc.tl.umap(
        adata,
        neighbors_key="neighbors_pca",
        random_state=int(pca_cfg["random_state"]),
    )
    adata.obsm["X_umap_pca"] = adata.obsm["X_umap"].copy()
    log_mem(logger, "after baseline UMAP")

    logger.info("Restoring counts for scPoli NB loss...")
    adata.X = adata.layers["counts"].copy()
    del adata.layers["counts"]
    gc.collect()

    n_epochs = int(train_cfg["n_epochs"])
    pretraining_epochs = int(train_cfg["pretraining_epochs"])

    if resume_from_model:
        logger.info("Loading trained scPoli from %s", model_dir)
        model = scPoli.load(str(model_dir), adata)
    else:
        logger.info(
            "Initializing scPoli (conditions=%s, cell_type=%s)...",
            scpoli_cfg["condition_keys"],
            scpoli_cfg["cell_type_key"],
        )
        model = scPoli(
            adata=adata,
            condition_keys=list(scpoli_cfg["condition_keys"]),
            cell_type_keys=scpoli_cfg["cell_type_key"],
            hidden_layer_sizes=list(scpoli_cfg["hidden_layer_sizes"]),
            latent_dim=int(scpoli_cfg["latent_dim"]),
            embedding_dims=int(scpoli_cfg["embedding_dims"]),
            recon_loss=scpoli_cfg["recon_loss"],
        )

        if not (0 < pretraining_epochs < n_epochs):
            raise ValueError(
                f"Need 0 < pretraining_epochs < n_epochs; got "
                f"pretraining_epochs={pretraining_epochs}, n_epochs={n_epochs}"
            )
        prototype_epochs = n_epochs - pretraining_epochs
        train_kwargs = dict(
            eta=float(train_cfg["eta"]),
            alpha_epoch_anneal=float(train_cfg["alpha_epoch_anneal"]),
            batch_size=int(train_cfg["batch_size"]),
            early_stopping_kwargs={
                "early_stopping_metric": train_cfg["early_stopping_metric"],
                "patience": int(train_cfg["patience"]),
                "threshold": float(train_cfg["threshold"]),
                "reduce_lr": bool(train_cfg["reduce_lr"]),
                "lr_patience": int(train_cfg["lr_patience"]),
                "lr_factor": float(train_cfg["lr_factor"]),
            },
        )

        logger.info(
            "Training scPoli phase 1/2: pretraining (%s epochs)...",
            pretraining_epochs,
        )
        model.train(
            n_epochs=pretraining_epochs,
            pretraining_epochs=pretraining_epochs,
            **train_kwargs,
        )
        save_model_checkpoint(
            model, pretrain_ckpt_dir, logger, "after-pretraining"
        )

        logger.info(
            "Training scPoli phase 2/2: prototypes + fine-tune (%s epochs)...",
            prototype_epochs,
        )
        model.train(
            n_epochs=prototype_epochs,
            pretraining_epochs=0,
            **train_kwargs,
        )
        save_model_checkpoint(model, model_dir, logger, "after-training")

    logger.info("Extracting latent embedding...")
    adata.obsm["X_emb"] = model.get_latent(adata, mean=True)

    logger.info("Neighbors + UMAP on X_emb...")
    sc.pp.neighbors(
        adata,
        use_rep="X_emb",
        key_added="neighbors_scpoli",
        n_neighbors=int(pca_cfg["n_neighbors"]),
        random_state=int(pca_cfg["random_state"]),
    )
    umap_kwargs = {
        "neighbors_key": "neighbors_scpoli",
        "random_state": int(pca_cfg["random_state"]),
    }
    if "min_dist" in pca_cfg:
        umap_kwargs["min_dist"] = float(pca_cfg["min_dist"])
    if "spread" in pca_cfg:
        umap_kwargs["spread"] = float(pca_cfg["spread"])
    sc.tl.umap(adata, **umap_kwargs)
    adata.obsm["X_umap_scpoli"] = adata.obsm["X_umap"].copy()

    adata.uns["scpoli"] = {
        "output_name": final_h5ad.name,
        "annotated_source": str(cfg["paths"]["annotated_h5ad"]),
        "params": {
            "condition_keys": list(scpoli_cfg["condition_keys"]),
            "cell_type_keys": scpoli_cfg["cell_type_key"],
            "hidden_layer_sizes": list(scpoli_cfg["hidden_layer_sizes"]),
            "latent_dim": int(scpoli_cfg["latent_dim"]),
            "embedding_dims": int(scpoli_cfg["embedding_dims"]),
            "n_epochs": n_epochs,
            "pretraining_epochs": pretraining_epochs,
            "eta": float(train_cfg["eta"]),
            "batch_size": int(train_cfg["batch_size"]),
            "n_top_genes": int(feat_cfg["n_top_genes"]),
            "checkpoint_pretrain_dir": str(pretrain_ckpt_dir),
            "model_dir": str(model_dir),
        },
    }

    logger.info("Saving integrated object to %s", final_h5ad)
    ensure_new(final_h5ad)
    adata.write_h5ad(final_h5ad, compression="gzip")
    logger.info("H5AD written: %s", final_h5ad)

    logger.info("Rendering comparison PNGs...")
    point_size = float(plot_cfg.get("point_size", 1.5))
    point_alpha = float(plot_cfg.get("point_alpha", 1.0))
    title_fs = int(plot_cfg.get("title_fontsize", 18))
    label_fs = int(plot_cfg.get("label_fontsize", 14))
    tick_fs = int(plot_cfg.get("tick_fontsize", 12))
    legend_fs = int(plot_cfg.get("legend_fontsize", 11))

    def _palette_for(key: str, base: dict | None) -> dict[str, str]:
        cats = list(adata.obs[key].astype("category").cat.categories)
        return palette_for_categories(cats, base or {})

    def generate_clean_comparison_plot(
        color_key, filename, title_suffix, palette=None, legend_loc="right margin"
    ):
        out = fig_dir / filename
        ensure_new(out)
        if palette is None:
            palette = _palette_for(color_key, None)
        elif not isinstance(palette, dict) or set(
            map(str, adata.obs[color_key].astype("category").cat.categories)
        ) - set(map(str, palette)):
            base = palette if isinstance(palette, dict) else {}
            palette = _palette_for(color_key, base)
        adata.uns[f"{color_key}_colors"] = [
            palette[str(c)]
            for c in adata.obs[color_key].astype("category").cat.categories
        ]
        fig, axes = plt.subplots(1, 2, figsize=tuple(plot_cfg["comparison_figsize"]))
        for ax, basis, prefix in (
            (axes[0], "X_umap_pca", "BEFORE: PCA"),
            (axes[1], "X_umap_scpoli", "AFTER: scPoli"),
        ):
            sc.pl.embedding(
                adata,
                basis=basis,
                color=color_key,
                ax=ax,
                legend_loc=legend_loc if ax is axes[1] else None,
                show=False,
                frameon=False,
                palette=palette,
                size=point_size,
                alpha=point_alpha,
                title=f"{prefix} ({title_suffix})",
            )
            ax.set_title(ax.get_title(), fontsize=title_fs, pad=12)
            ax.set_xlabel(ax.get_xlabel(), fontsize=label_fs)
            ax.set_ylabel(ax.get_ylabel(), fontsize=label_fs)
            ax.tick_params(labelsize=tick_fs)
        leg = axes[1].get_legend()
        if leg is not None:
            for text in leg.get_texts():
                text.set_fontsize(legend_fs)
            leg.set_title(color_key, prop={"size": label_fs})
        fig.tight_layout()
        fig.savefig(out, dpi=int(plot_cfg["dpi"]), bbox_inches="tight", facecolor="white")
        plt.close(fig)
        logger.info("saved %s", out)

    for key in (
        "sample_id",
        "study",
        "Level_1",
        "Level_2",
        "Level_3",
        "Age",
        "diff_protocol",
        "sc_protocol",
    ):
        if key in adata.obs.columns:
            adata.obs[key] = adata.obs[key].astype(str).astype("category")

    generate_clean_comparison_plot(
        "sample_id", "Comparison_Sample_Alignment.png", "Sample ID", legend_loc=None
    )
    generate_clean_comparison_plot("study", "Comparison_Study_Alignment.png", "Study")
    generate_clean_comparison_plot(
        "Level_1",
        "Comparison_Level1_Biology.png",
        "Level 1 Broad Categories",
        palette=LEVEL_1_COLORS,
    )
    generate_clean_comparison_plot(
        "Level_2",
        "Comparison_Level2_Biology.png",
        "Level 2 Cell Types",
        palette=LEVEL_2_COLORS,
    )
    generate_clean_comparison_plot(
        "Level_3",
        "Comparison_Level3_Biology.png",
        "Level 3 Cell Types",
        palette=LEVEL_3_COLORS,
    )
    for key, filename, title, base in (
        ("Age", "Comparison_Age.png", "Age", None),
        (
            "diff_protocol",
            "Comparison_DiffProtocol.png",
            "Differentiation Protocol",
            DIFF_PROTOCOL_COLORS,
        ),
        (
            "sc_protocol",
            "Comparison_ScTechnology.png",
            "Single-cell Technology",
            SC_PROTOCOL_COLORS,
        ),
    ):
        if key not in adata.obs.columns:
            logger.warning("Skipping plot; obs key missing: %s", key)
            continue
        generate_clean_comparison_plot(key, filename, title, palette=base)

    logger.info("DONE")


if __name__ == "__main__":
    main()
