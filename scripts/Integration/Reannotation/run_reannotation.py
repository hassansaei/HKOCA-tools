#!/usr/bin/env python3
"""Apply post-integration cluster reannotation to the scPoli atlas.

Reads the integrated object and the audited before/after label table, applies
the Level_1/2/3 updates, and writes a new reannotated object. The source
integrated file is never overwritten.
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path("/data-master/pure-workspace/labss/hmami/atlas_shortcut_new_start")
AUDIT = ROOT / "level1_lineage_audit"
INV = AUDIT / "investigations"
REAN = ROOT / "reannotation"

SRC_INT = ROOT / "scpoli_integrated_leiden_res2.0_clean.h5ad"
OUT_INT = ROOT / "scpoli_integrated_leiden_res2.0_clean_reannotated.h5ad"
AFTER_LAB = INV / "00_umap_overlap_retune" / "tables" / "after_labels_barcode.tsv.gz"
C2_SUB = INV / "cache" / "nephron_c2_subclusters.csv"

LINEAGES = [
    "Nephron",
    "Stromal_cells",
    "Off_target_cells",
    "Endothelial_cells",
    "Ureteric_Epithelium",
    "undifferentiated_pluripotent_cells",
]

# Cluster / subcluster remaps (ordered by cell count).
CLUSTER_REANNOTATIONS: list[dict] = [
    {
        "lineage_before": "Stromal_cells",
        "cluster": "13",
        "n_cells": 14413,
        "to_Level_1": "Nephron",
        "to_Level_2": "Nephron_progenitor_cells",
        "to_Level_3": "Nephron_progenitor_cells",
    },
    {
        "lineage_before": "Off_target_cells",
        "cluster": "10",
        "n_cells": 8252,
        "to_Level_1": "Nephron",
        "to_Level_2": "Proliferative_nephron",
        "to_Level_3": "Proliferative_nephron",
    },
    {
        "lineage_before": "Nephron",
        "cluster": "2_sub0",
        "n_cells": 7439,
        "to_Level_1": "Stromal_cells",
        "to_Level_2": "S2",
        "to_Level_3": "S2",
    },
    {
        "lineage_before": "Nephron",
        "cluster": "14",
        "n_cells": 6143,
        "to_Level_1": "Nephron",
        "to_Level_2": "Proliferative_nephron",
        "to_Level_3": "Proliferative_nephron",
    },
    {
        "lineage_before": "Nephron",
        "cluster": "2_sub1",
        "n_cells": 5581,
        "to_Level_1": "Nephron",
        "to_Level_2": "Proliferative_nephron",
        "to_Level_3": "Proliferative_nephron",
    },
    {
        "lineage_before": "Nephron",
        "cluster": "2_sub2",
        "n_cells": 5237,
        "to_Level_1": "Nephron",
        "to_Level_2": "Nephron_progenitor_cells",
        "to_Level_3": "Nephron_progenitor_cells",
    },
    {
        "lineage_before": "Stromal_cells",
        "cluster": "20",
        "n_cells": 4305,
        "to_Level_1": "Nephron",
        "to_Level_2": "Nephron_progenitor_cells",
        "to_Level_3": "Nephron_progenitor_cells",
    },
    {
        "lineage_before": "Nephron",
        "cluster": "19",
        "n_cells": 3969,
        "to_Level_1": "Nephron",
        "to_Level_2": "Proliferative_nephron",
        "to_Level_3": "Proliferative_nephron",
    },
    {
        "lineage_before": "Nephron",
        "cluster": "25",
        "n_cells": 2006,
        "to_Level_1": "Off_target_cells",
        "to_Level_2": "Off_target_cells",
        "to_Level_3": "Off_target_cells",
    },
    {
        "lineage_before": "Off_target_cells",
        "cluster": "29",
        "n_cells": 1448,
        "to_Level_1": "Stromal_cells",
        "to_Level_2": "S2",
        "to_Level_3": "S2",
    },
    {
        "lineage_before": "Stromal_cells",
        "cluster": "26",
        "n_cells": 110,
        "to_Level_1": "Ureteric_Epithelium",
        "to_Level_2": "Collecting_duct",
        "to_Level_3": "Collecting_duct",
    },
    {
        "lineage_before": "undifferentiated_pluripotent_cells",
        "cluster": "32",
        "n_cells": 93,
        "to_Level_1": "Endothelial_cells",
        "to_Level_2": "Endothelial_cells",
        "to_Level_3": "Endothelial_cells",
    },
    {
        "lineage_before": "Off_target_cells",
        "cluster": "34",
        "n_cells": 37,
        "to_Level_1": "Nephron",
        "to_Level_2": "Proximal_tubule",
        "to_Level_3": "Dev_PT",
    },
    {
        "lineage_before": "Off_target_cells",
        "cluster": "35",
        "n_cells": 27,
        "to_Level_1": "Endothelial_cells",
        "to_Level_2": "Endothelial_cells",
        "to_Level_3": "Endothelial_cells",
    },
]


def load_lineage_leiden_and_c2(barcodes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    parts = []
    for lin in LINEAGES:
        obj = next((AUDIT / lin / "objects").glob("*_leiden1.5.h5ad"))
        a = ad.read_h5ad(obj, backed="r")
        parts.append(
            pd.DataFrame(
                {
                    "barcode": a.obs_names.astype(str),
                    "lineage_before": lin,
                    "leiden_res1.5": a.obs["leiden_res1.5"].astype(str).to_numpy(),
                }
            )
        )
        a.file.close()
    meta = pd.concat(parts, ignore_index=True).set_index("barcode").reindex(barcodes)

    c2 = pd.read_csv(C2_SUB)
    c2["barcode"] = c2["barcode"].astype(str)
    sub_col = "c2_subcluster" if "c2_subcluster" in c2.columns else "c2_sub"
    c2_map = (
        c2.set_index("barcode")[sub_col]
        .astype(str)
        .str.replace("c2_sub", "", regex=False)
        .to_dict()
    )
    c2_sub = np.array([c2_map.get(b, "NA") for b in barcodes], dtype=object)
    return (
        meta["lineage_before"].astype(str).to_numpy(),
        meta["leiden_res1.5"].astype(str).to_numpy(),
        c2_sub,
    )


def build_rule_per_barcode(lin_b, cl_b, c2_sub, b1, a1, b2, a2, b3, a3) -> np.ndarray:
    n = len(b1)
    rule = np.full(n, "unchanged", dtype=object)
    changed = (b1 != a1) | (b2 != a2) | (b3 != a3)

    for sub in ("0", "1", "2"):
        m = (lin_b == "Nephron") & (cl_b == "2") & (c2_sub == sub) & changed
        if m.any():
            rule[m] = f"Nephron_c2_sub{sub}->{a1[m][0]}/{a2[m][0]}"

    for entry in CLUSTER_REANNOTATIONS:
        cl = entry["cluster"]
        if "_" in cl:
            continue
        lin = entry["lineage_before"]
        m = (lin_b == lin) & (cl_b == cl) & changed
        if m.any():
            rule[m] = f"{lin}_c{cl}->{a1[m][0]}/{a2[m][0]}"

    rem = changed & (rule == "unchanged")
    if rem.any():
        for i in np.flatnonzero(rem):
            rule[i] = f"{lin_b[i]}_c{cl_b[i]}->{a1[i]}/{a2[i]}"
    return rule


def main() -> None:
    (REAN / "tables").mkdir(parents=True, exist_ok=True)

    lab = pd.read_csv(AFTER_LAB, sep="\t")
    lab["barcode"] = lab["barcode"].astype(str)
    lab_map = lab.set_index("barcode")

    print(f"Loading integrated h5ad: {SRC_INT}")
    adata = ad.read_h5ad(SRC_INT)
    barcodes = adata.obs_names.astype(str).to_numpy()
    aligned = lab_map.reindex(barcodes)
    if aligned["Level_1_after"].isna().any():
        raise ValueError("Missing barcodes in AFTER_LAB")

    b1 = aligned["Level_1_before"].astype(str).to_numpy()
    b2 = aligned["Level_2_before"].astype(str).to_numpy()
    b3 = adata.obs["Level_3"].astype(str).to_numpy()
    a1 = aligned["Level_1_after"].astype(str).to_numpy()
    a2 = aligned["Level_2_after"].astype(str).to_numpy()
    a3_from_tab = aligned["Level_3_after"].astype(str).to_numpy()

    lin_b, cl_b, c2_sub = load_lineage_leiden_and_c2(barcodes)
    changed_l12 = (b1 != a1) | (b2 != a2)
    a3 = b3.copy()
    a3[changed_l12] = a3_from_tab[changed_l12]
    rule = build_rule_per_barcode(lin_b, cl_b, c2_sub, b1, a1, b2, a2, b3, a3)
    ch1 = b1 != a1
    ch2 = b2 != a2
    ch3 = b3 != a3
    ch_any = ch1 | ch2 | ch3

    adata.obs["Level_1_before_reannotation"] = adata.obs["Level_1"].astype(str)
    adata.obs["Level_2_before_reannotation"] = adata.obs["Level_2"].astype(str)
    adata.obs["Level_3_before_reannotation"] = b3
    if "Level_3_latest" in adata.obs:
        adata.obs["Level_3_latest_before_reannotation"] = adata.obs[
            "Level_3_latest"
        ].astype(str)

    adata.obs["Level_1"] = pd.Categorical(a1)
    adata.obs["Level_2"] = pd.Categorical(a2)
    adata.obs["Level_3"] = pd.Categorical(a3)
    adata.obs["reannotation_changed"] = ch_any
    adata.obs["reannotation_rule"] = pd.Categorical(rule)
    adata.obs["nephron_c2_subcluster"] = pd.Categorical(c2_sub)

    adata.uns["reannotation_metadata"] = {
        "promoted_at": datetime.now().isoformat(timespec="seconds"),
        "source_integrated_h5ad": str(SRC_INT),
        "n_cells_total": int(adata.n_obs),
        "n_cells_changed_L1": int(ch1.sum()),
        "n_cells_changed_L2": int(ch2.sum()),
        "n_cells_changed_L3": int(ch3.sum()),
        "n_cells_changed_any": int(ch_any.sum()),
    }

    print(f"Writing reannotated h5ad → {OUT_INT}")
    adata.write_h5ad(OUT_INT, compression="gzip")

    per_bc = pd.DataFrame(
        {
            "barcode": barcodes,
            "lineage_before": lin_b,
            "leiden_res1.5_in_lineage": cl_b,
            "nephron_c2_subcluster": c2_sub,
            "Level_1_before": b1,
            "Level_2_before": b2,
            "Level_3_before": b3,
            "Level_1_after": a1,
            "Level_2_after": a2,
            "Level_3_after": a3,
            "changed_L1": ch1,
            "changed_L2": ch2,
            "changed_L3": ch3,
            "changed_any": ch_any,
            "reannotation_rule": rule,
        }
    )
    per_bc.to_csv(REAN / "tables" / "reannotation_per_barcode.tsv.gz", sep="\t", index=False)

    pd.DataFrame(CLUSTER_REANNOTATIONS).to_csv(
        REAN / "tables" / "reannotation_per_cluster_summary.csv", index=False
    )
    print(
        f"DONE: L1={ch1.sum():,} L2={ch2.sum():,} L3={ch3.sum():,} "
        f"any={ch_any.sum():,} → {OUT_INT}"
    )


if __name__ == "__main__":
    main()
