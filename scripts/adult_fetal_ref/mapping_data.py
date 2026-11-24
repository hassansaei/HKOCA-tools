"""Load adult/fetal reference + organoid latents for wkNN mapping / presence."""
from __future__ import annotations

import h5py
import numpy as np

from lib_config import PATHS, WORKSPACE, read_h5_obs_col

HCKA_FULL = (
    WORKSPACE
    / "human_kidney_reference_projection"
    / "Single_cell_human_kidney_cell_atlas.h5ad"
)
ORG_ON_HCKA = (
    WORKSPACE
    / "human_kidney_reference_projection"
    / "04_projection"
    / "organoid_projected_HCKA_healthy_adjacent.h5ad"
)


def load_adult() -> tuple:
    """Return (ref_lat, ref_umap, ref_labels, org_lat, org_l2, label_key)."""
    with h5py.File(HCKA_FULL, "r") as f:
        barcodes = read_h5_obs_col(f, "_index")
        disease = read_h5_obs_col(f, "disease")
        site = np.array(
            [str(x).strip().lower() for x in read_h5_obs_col(f, "sampled_site_condition")],
            dtype=object,
        )
        broad = read_h5_obs_col(f, "broad_annotation")
        umap_full = np.asarray(f["obsm"]["X_scANVI_umap"][:], dtype=np.float32)
    keep = (disease == "normal") & np.isin(site, ["healthy", "adjacent"])
    ref_umap = umap_full[keep]
    ref_broad = broad[keep]
    ref_bc = barcodes[keep]

    with h5py.File(PATHS["adult_meta_broad"], "r") as f:
        lat_bc = read_h5_obs_col(f, "_index")
        ref_lat = np.asarray(f["obsm"]["X_scANVI"][:], dtype=np.float32)
    pos = {b: i for i, b in enumerate(ref_bc)}
    order = np.array([pos[b] for b in lat_bc], dtype=np.int64)
    ref_umap = ref_umap[order]
    ref_broad = ref_broad[order]

    with h5py.File(ORG_ON_HCKA, "r") as fo:
        org_lat = np.asarray(fo["obsm"]["X_scANVI"][:], dtype=np.float32)
        org_l2 = read_h5_obs_col(fo, "Level_2")
    return ref_lat, ref_umap, ref_broad, org_lat, org_l2, "broad_annotation"


def load_fetal() -> tuple:
    """Return (ref_lat, ref_umap, ref_labels, org_lat, org_l2, label_key)."""
    with h5py.File(PATHS["fetal_harmony"], "r") as fh:
        harm_bc = read_h5_obs_col(fh, "_index")
        harm_umap = np.asarray(fh["obsm"]["X_umap"][:], dtype=np.float32)

    with h5py.File(PATHS["fetal_umap"], "r") as fs:
        scanvi_lat = np.asarray(fs["obsm"]["X_scANVI"][:], dtype=np.float32)
        scanvi_ibc = read_h5_obs_col(fs, "integrated_barcode")
        scanvi_broad = read_h5_obs_col(fs, "harmonized_celltype")

    harm_pos = {b: i for i, b in enumerate(harm_bc)}
    scanvi_to_harm = np.array([harm_pos[b] for b in scanvi_ibc], dtype=np.int64)
    ref_umap = harm_umap[scanvi_to_harm]

    with h5py.File(PATHS["org_on_fetal"], "r") as fo:
        org_lat = np.asarray(fo["obsm"]["X_scANVI"][:], dtype=np.float32)
        org_l2 = read_h5_obs_col(fo, "Level_2")
    return scanvi_lat, ref_umap, scanvi_broad, org_lat, org_l2, "harmonized_celltype"
