"""scPoli surgery helpers (scArches) for HKOCA atlas projection."""

from __future__ import annotations

import io
import logging
import pickle
import sys
import types
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger("hkoca.projection")

UNKNOWN_LABEL = "Unknown"
EARLY_STOPPING_KWARGS = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}


def patch_anndata_for_scarches() -> None:
    """Restore anndata APIs that scarches 0.6 still imports."""
    import anndata as ad

    if not hasattr(ad, "read"):
        ad.read = ad.read_h5ad
    if not hasattr(ad, "io"):
        io_module = types.ModuleType("anndata.io")
        for name in (
            "read_csv", "read_loom", "read_text", "read_excel", "read_hdf",
            "read_mtx", "read_zarr", "read_h5ad", "read_umi_tools", "read_elem",
            "write_elem",
        ):
            if hasattr(ad, name):
                setattr(io_module, name, getattr(ad, name))
        ad.io = io_module
        sys.modules["anndata.io"] = io_module
    if not hasattr(ad, "abc"):
        from anndata._core.sparse_dataset import CSCDataset, CSRDataset

        abc_module = types.ModuleType("anndata.abc")
        abc_module.CSRDataset = CSRDataset
        abc_module.CSCDataset = CSCDataset
        ad.abc = abc_module
        sys.modules["anndata.abc"] = abc_module


def is_counts_like(x, max_check: int = 200_000) -> bool:
    if sp.issparse(x):
        data = x.data
        if data.size == 0:
            return True
        if data.size > max_check:
            rng = np.random.default_rng(0)
            data = data[rng.choice(data.size, size=max_check, replace=False)]
    else:
        flat = np.asarray(x).ravel()
        if flat.size > max_check:
            rng = np.random.default_rng(0)
            flat = flat[rng.choice(flat.size, size=max_check, replace=False)]
        data = flat
    data = np.asarray(data, dtype=np.float64)
    if np.any(data < -1e-6):
        return False
    return bool(np.nanmax(np.abs(data - np.round(data))) < 1e-3)


def ensure_raw_counts(adata, prefer_layer: str = "counts"):
    import anndata as ad

    out = adata.copy()
    if prefer_layer in out.layers and is_counts_like(out.layers[prefer_layer]):
        out.X = out.layers[prefer_layer].copy()
        src = f"layers[{prefer_layer!r}]"
    elif is_counts_like(out.X):
        src = ".X"
    elif prefer_layer in out.layers:
        out.X = out.layers[prefer_layer].copy()
        src = f"layers[{prefer_layer!r}] (heuristic soft-fail)"
    else:
        raise ValueError(
            "Could not find raw counts. Provide integer UMIs in layers['counts'] "
            "or raw .X (SCT residuals are not valid scPoli input)."
        )
    if not sp.issparse(out.X):
        out.X = sp.csr_matrix(out.X)
    else:
        out.X = out.X.tocsr()
    logger.info("Raw counts source: %s | shape=%s nnz=%s", src, out.X.shape, out.X.nnz)
    return out


def align_to_reference_genes(adata, ref_genes: list[str]):
    import anndata as ad

    adata = adata.copy()
    adata.var_names_make_unique()
    q_genes = adata.var_names.astype(str)
    q_set = set(q_genes)
    shared = [g for g in ref_genes if g in q_set]
    missing = [g for g in ref_genes if g not in q_set]
    extra = [g for g in q_genes if g not in set(ref_genes)]
    logger.info(
        "Gene overlap with reference model: shared=%d / %d (%.1f%%) | missing=%d | dropped=%d",
        len(shared), len(ref_genes), 100.0 * len(shared) / max(len(ref_genes), 1),
        len(missing), len(extra),
    )
    if len(shared) < 0.5 * len(ref_genes):
        raise ValueError(
            f"Too few shared genes ({len(shared)}/{len(ref_genes)}). "
            "Check that query var_names are human gene symbols."
        )
    ad_shared = adata[:, shared].copy()
    if missing:
        pad = ad.AnnData(
            X=sp.csr_matrix((ad_shared.n_obs, len(missing)), dtype=np.float32),
            obs=ad_shared.obs.copy(),
            var=pd.DataFrame(index=pd.Index(missing, name=ad_shared.var.index.name)),
        )
        aligned = ad.concat([ad_shared, pad], axis=1, merge="same", join="outer")
    else:
        aligned = ad_shared
    aligned = aligned[:, ref_genes].copy()
    aligned.var_names = pd.Index(ref_genes)
    if not sp.issparse(aligned.X):
        aligned.X = sp.csr_matrix(aligned.X)
    return aligned


def scpoli_placeholder_col(model_dir: Path) -> str:
    import torch

    with open(Path(model_dir) / "attr.pkl", "rb") as fh:
        blob = fh.read()

    class _CPUUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            if module == "torch.storage" and name == "_load_from_bytes":
                def _load_from_bytes(b):
                    try:
                        return torch.load(io.BytesIO(b), map_location="cpu", weights_only=False)
                    except TypeError:
                        return torch.load(io.BytesIO(b), map_location="cpu")

                return _load_from_bytes
            return super().find_class(module, name)

    attrs = _CPUUnpickler(io.BytesIO(blob)).load()
    return attrs["cell_type_keys_"][0]


def load_reference_genes(model_dir: Path) -> list[str]:
    csv_path = Path(model_dir) / "var_names.csv"
    if not csv_path.is_file():
        raise FileNotFoundError(f"Missing var_names.csv in model dir: {model_dir}")
    return pd.read_csv(csv_path, header=None)[0].astype(str).tolist()


def validate_model_dir(model_dir: Path) -> None:
    model_dir = Path(model_dir)
    missing = [n for n in ("model_params.pt", "attr.pkl", "var_names.csv") if not (model_dir / n).is_file()]
    if missing:
        raise FileNotFoundError(
            f"scPoli model dir {model_dir} is missing: {', '.join(missing)}"
        )


def majority_map(obs_df: pd.DataFrame, from_col: str, to_col: str) -> dict[str, str]:
    ct = (
        obs_df.groupby([from_col, to_col], observed=True)
        .size()
        .rename("n")
        .reset_index()
        .sort_values([from_col, "n"], ascending=[True, False])
    )
    return ct.drop_duplicates(from_col).set_index(from_col)[to_col].astype(str).to_dict()


def load_prototype_bridge(tsv_path: Path, celltype_key: str) -> dict[str, str]:
    df = pd.read_csv(tsv_path, sep="\t")
    if "prototype_label" not in df.columns or celltype_key not in df.columns:
        raise ValueError(f"Unexpected columns in {tsv_path}: {list(df.columns)}")
    return df.set_index("prototype_label")[celltype_key].astype(str).to_dict()


def _condition_tensor(model, adata, device):
    import torch

    label_tensor = []
    for cond in model.condition_keys_:
        query_conditions = adata.obs[cond].values
        if not set(query_conditions).issubset(model.conditions_[cond]):
            missing = set(query_conditions) - set(model.conditions_[cond])
            raise ValueError(f"Unknown {cond} values for surgery model: {sorted(missing)[:10]}")
        labels = np.zeros(query_conditions.shape[0])
        for condition, label in model.model.condition_encoders[cond].items():
            labels[query_conditions == condition] = label
        label_tensor.append(labels)
    return torch.tensor(np.stack(label_tensor), device=device).T


def get_latent_safe(model, adata, mean: bool = True) -> np.ndarray:
    import torch

    device = next(model.model.parameters()).device
    x = adata.X
    c = _condition_tensor(model, adata, device)
    latents = []
    n = x.shape[0]
    for start in range(0, n, 512):
        batch = np.arange(start, min(start + 512, n))
        x_batch = x[batch, :]
        if sp.issparse(x_batch):
            x_batch = x_batch.toarray()
        x_batch = torch.tensor(np.asarray(x_batch), device=device).float()
        latent = model.model.get_latent(x_batch, c[batch, :], mean)
        latents.append(latent.cpu().detach())
    return torch.cat(latents).numpy()


def classify_safe(model, adata, scale_uncertainties: bool = True) -> dict[str, Any]:
    import torch
    from sklearn.preprocessing import MinMaxScaler, RobustScaler

    if model.prototypes_labeled_["mean"] is None:
        raise ValueError("Model has no labeled prototypes.")
    device = next(model.model.parameters()).device
    model.model.eval()
    x = adata.X
    c = _condition_tensor(model, adata, device)
    results: dict[str, Any] = {}
    n = x.shape[0]
    for cell_type_key in model.cell_type_keys_:
        prototypes_idx = [
            i for i, key in enumerate(model.cell_types_.keys())
            if cell_type_key in model.cell_types_[key]
        ]
        prototypes_idx = torch.tensor(prototypes_idx, device=device)
        preds, uncert = [], []
        for start in range(0, n, 512):
            batch = np.arange(start, min(start + 512, n))
            x_batch = x[batch, :]
            if sp.issparse(x_batch):
                x_batch = x_batch.toarray()
            x_batch = torch.tensor(np.asarray(x_batch), device=device).float()
            pred, prob, _ = model.model.classify(
                x_batch,
                c[batch].to(device),
                prototype=False,
                classes_list=prototypes_idx,
                p=2,
                get_prob=False,
                log_distance=True,
            )
            preds.append(pred.cpu().detach())
            uncert.append(prob.cpu().detach())
        full_pred = np.array(torch.cat(preds))
        full_uncert = np.array(torch.cat(uncert))
        inv_ct_encoder = {v: k for k, v in model.model.cell_type_encoder.items()}
        full_pred_names = [inv_ct_encoder[int(p)] for p in full_pred]
        if scale_uncertainties:
            full_uncert = RobustScaler().fit_transform(full_uncert.reshape(-1, 1))
            full_uncert = MinMaxScaler(feature_range=(0, 1)).fit_transform(full_uncert).reshape(-1)
        results[cell_type_key] = {"preds": np.array(full_pred_names), "uncert": full_uncert}
    return results


def load_query_model(query, model_dir: Path, unknown_label: str = UNKNOWN_LABEL):
    """Load scPoli query surgery model (scarches DataFrame setitem workaround)."""
    import pandas as pd
    import torch
    from scarches.models.scpoli import scPoli

    patch_anndata_for_scarches()
    map_location = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info("Loading scPoli reference weights from %s (%s)", model_dir, map_location)
    orig_setitem = pd.DataFrame.__setitem__

    def _setitem(self, key, value):
        if key == "conditions_combined" and isinstance(value, pd.DataFrame) and value.shape[1] == 1:
            value = value.iloc[:, 0]
        return orig_setitem(self, key, value)

    pd.DataFrame.__setitem__ = _setitem
    try:
        return scPoli.load_query_data(
            adata=query,
            reference_model=str(model_dir),
            labeled_indices=[],
            unknown_ct_names=[unknown_label],
            map_location=map_location,
        )
    finally:
        pd.DataFrame.__setitem__ = orig_setitem


def train_query_model(
    model,
    *,
    n_epochs: int = 50,
    pretrain_epochs: int = 40,
    eta: float = 10,
) -> None:
    logger.info(
        "Training scPoli query surgery (n_epochs=%s, pretraining_epochs=%s, eta=%s)",
        n_epochs,
        pretrain_epochs,
        eta,
    )
    model.train(
        n_epochs=n_epochs,
        pretraining_epochs=pretrain_epochs,
        eta=eta,
        early_stopping_kwargs=EARLY_STOPPING_KWARGS,
        use_early_stopping=True,
        unlabeled_prototype_training=False,
    )


def load_atlas_obs(atlas_h5ad: Path) -> pd.DataFrame:
    import scanpy as sc

    atlas_b = sc.read_h5ad(Path(atlas_h5ad), backed="r")
    try:
        return atlas_b.obs.copy()
    finally:
        atlas_b.file.close()


def read_atlas_obsm(atlas_h5ad: Path, key: str) -> np.ndarray | None:
    import h5py

    with h5py.File(Path(atlas_h5ad), "r") as fh:
        if "obsm" not in fh or key not in fh["obsm"]:
            return None
        return np.asarray(fh["obsm"][key][:], dtype=np.float32)


def resolve_atlas_latent_key(atlas_h5ad: Path, preferred: str = "X_scpoli") -> str | None:
    for key in (preferred, "X_scpoli", "X_emb"):
        arr = read_atlas_obsm(atlas_h5ad, key)
        if arr is not None:
            return key
    return None


def load_atlas_subsample_h5py(
    h5ad_path: Path,
    row_idx: np.ndarray,
    obs_sub: pd.DataFrame,
    barcodes: np.ndarray,
    ref_genes: list[str],
):
    """Load atlas CSR counts for row indices via h5py (backed subset is unreliable)."""
    import anndata as ad
    import h5py

    row_idx = np.asarray(row_idx, dtype=np.int64)
    with h5py.File(h5ad_path, "r") as fh:
        if "counts" in fh.get("layers", {}):
            g = fh["layers"]["counts"]
        else:
            g = fh["X"]
        shape = tuple(int(x) for x in g.attrs["shape"])
        indptr = np.asarray(g["indptr"][:], dtype=np.int64)
        starts = indptr[row_idx]
        ends = indptr[row_idx + 1]
        nnz = int((ends - starts).sum())
        data = np.empty(nnz, dtype=np.float32)
        indices = np.empty(nnz, dtype=np.int32)
        new_indptr = np.zeros(len(row_idx) + 1, dtype=np.int64)
        pos = 0
        for i, (s, e) in enumerate(zip(starts, ends, strict=True)):
            n = int(e - s)
            if n:
                data[pos : pos + n] = g["data"][s:e]
                indices[pos : pos + n] = g["indices"][s:e]
            pos += n
            new_indptr[i + 1] = pos
        x_sub = sp.csr_matrix((data, indices, new_indptr), shape=(len(row_idx), shape[1]))
        var_key = "features" if "features" in fh["var"] else "_index"
        genes = [x.decode() if isinstance(x, bytes) else str(x) for x in fh["var"][var_key][:]]

    atlas_sub = ad.AnnData(X=x_sub, obs=obs_sub.copy(), var=pd.DataFrame(index=pd.Index(genes)))
    atlas_sub.obs_names = pd.Index(barcodes.astype(str))
    return align_to_reference_genes(atlas_sub, ref_genes)
