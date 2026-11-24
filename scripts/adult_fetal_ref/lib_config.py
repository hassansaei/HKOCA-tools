"""Shared configuration and utilities for organoid_vs_primary_atlas."""
from __future__ import annotations

import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# Workspace & output roots
WORKSPACE = Path("/data-master/pure-workspace/labss/hmami/atlas_shortcut_new_start")
ROOT = WORKSPACE / "organoid_vs_primary_atlas"

# Input objects (FULL references WITH immune — no "no_immune" objects)
PATHS = {
    # 600k adult healthy+adjacent reference (598,928 cells)
    "adult_counts": WORKSPACE
    / "organoid_on_primary_projection/01_reference/adult_counts_shared.h5ad",
    "adult_umap": WORKSPACE
    / "organoid_on_primary_projection/01_reference/adult_scANVI_latent_umap.h5ad",
    "adult_meta_broad": WORKSPACE
    / "human_kidney_reference_projection/01_reference/HCKA_healthy_adjacent_reference_scANVI_latent.h5ad",
    # Fetal weeks 9-17 integrated reference (47,999 cells)
    "fetal_counts": WORKSPACE
    / "organoid_on_primary_projection/01_reference/fetal_counts_shared.h5ad",
    "fetal_umap": WORKSPACE
    / "organoid_on_primary_projection/01_reference/fetal_scANVI_latent_umap.h5ad",
    "fetal_harmony": WORKSPACE
    / "fetal_new_3_datasets/integrated_weeks_9_17/fetal_weeks_9_17_harmony_integrated.h5ad",
    # Organoid atlas & projections onto adult and fetal references (847,399 cells)
    "organoid_counts": WORKSPACE
    / "organoid_on_primary_projection/02_query/organoid_counts_shared.h5ad",
    "org_on_adult": WORKSPACE
    / "organoid_on_primary_projection/04_projection/organoid_on_adult_umap.h5ad",
    "org_on_fetal": WORKSPACE
    / "organoid_on_primary_projection/04_projection/organoid_on_fetal_umap.h5ad",
    # External gene resources
    "tf_csv": ROOT
    / "scripts/03_similarity_score/TF_names_Lambert_v1.01.csv",
    "gmt_dir": Path("/home/hmami/.cache/gseapy"),
}

# Clean analysis subfolders
DIRS = {
    "00_fig": ROOT / "00_reference_umaps/figures",
    "00_tab": ROOT / "00_reference_umaps/tables",
    "01_fig": ROOT / "01_presence_coverage/figures",
    "01_tab": ROOT / "01_presence_coverage/tables",
    "02_fig": ROOT / "02_differential_expression/figures",
    "02_tab": ROOT / "02_differential_expression/tables",
    "02_path_fig": ROOT / "02_differential_expression/pathways/figures",
    "02_path_tab": ROOT / "02_differential_expression/pathways/tables",
    "02_gsea_fig": ROOT / "02_differential_expression/gsea/figures",
    "02_gsea_tab": ROOT / "02_differential_expression/gsea/tables",
    "03_fig": ROOT / "03_similarity_score/figures",
    "03_tab": ROOT / "03_similarity_score/tables",
}

RSCRIPT = Path("/data-master/pure-workspace/hmami/conda_envs/atlas/bin/Rscript")
EDGER_R = ROOT / "scripts/r/edger_qlf_org_vs_primary.R"

# Protocols merged into the 3 major protocols for DEG / protocol comparisons
PROTOCOLS_3 = ["Morizane", "Takasato", "Freedman"]

PROTOCOL_COLORS = {
    "Morizane": "#0072B2",
    "Takasato": "#D55E00",
    "Freedman": "#009E73",
    "All organoids": "#333333",
    "Other": "#999999",
}

REF_COLORS = {
    "fetal": "#0072B2",
    "adult": "#D55E00",
}

# Maximally distinct organoid Level_2 colors (no near-shades)
# Maximally distinct Level_2 colors (no near-orange / near-red pairs)
ORG_L2_COLORS = {
    "Nephron_progenitor_cells": "#1A237E",
    "Podocytes": "#E41A1C",
    "Proximal_tubule": "#FF7F00",
    "Distal_tubule": "#377EB8",
    "LOH": "#4DAF4A",
    "Collecting_duct": "#984EA3",
    "PEC-like": "#00BCD4",
    "Proliferative_nephron": "#F781BF",
    "Proliferative_stroma": "#A6D854",
    "S1": "#8B4513",
    "S2": "#6A3D9A",
    "S3": "#FFD92F",
    "Endothelial_cells": "#212121",
    "Off_target_cells": "#9E9E9E",
    "undifferentiated_pluripotent_cells": "#80DEEA",
}

# Adult broad types collapsed into one "immune cells" row for presence plots
ADULT_IMMUNE_TYPES = {
    "CD8 T cells",
    "CD4 T cells",
    "T cells",
    "B cells",
    "Plasma cells",
    "Macrophages",
    "Monocytes",
    "NK cells",
    "cDCs",
    "pDCs",
    "Mast cells",
    "Neutrophils",
    "RBC",
    "Platelets",
}

# Fetal types collapsed into "immune cells"
FETAL_IMMUNE_TYPES = {
    "Immune",
    "Megakaryocyte",
    "Erythroid",
}

# Exact 9 organoid cell types requested for DE & similarity
TARGET_CELLTYPES_9 = [
    "Distal_tubule",
    "Nephron_progenitor_cells",
    "PEC-like",
    "Podocytes",
    "Proliferative_nephron",
    "Proliferative_stroma",
    "Proximal_tubule",
    "S1",
    "S2",
]

# Organoid Level_2 -> (fetal harmonized_celltype list, adult fine_annotation list)
PAIRS_L2 = {
    "Distal_tubule": (
        ["Distal_tubule_LOH", "Connecting_tubule"],
        ["DCT1", "DCT2"],
    ),
    "Nephron_progenitor_cells": (
        ["NPC"],
        [],
    ),
    "PEC-like": (
        ["PEC"],
        ["PEC"],
    ),
    "Podocytes": (
        ["Podocytes", "Early_podocyte"],
        ["POD"],
    ),
    "Proliferative_nephron": (
        ["NPC", "PTA", "Renal_vesicle_CSB", "SSB"],
        [],
    ),
    "Proliferative_stroma": (
        ["Fibroblast", "Stroma_progenitor", "Myofibroblast"],
        [
            "Fibroblasts_NEGR1+",
            "Fibroblasts_INHBA+",
            "Fibroblasts_adventitial",
            "Fibroblasts_cortex",
            "Medullary_fibroblasts",
            "Pericytes",
            "MC",
            "Myofibroblasts",
        ],
    ),
    "Proximal_tubule": (
        ["Proximal_tubule"],
        ["PTS1", "PTS2", "PTS3", "PTS1/S2", "PTS2/S3"],
    ),
    "S1": (
        ["Fibroblast", "Stroma_progenitor", "Myofibroblast"],
        [
            "Fibroblasts_NEGR1+",
            "Fibroblasts_INHBA+",
            "Fibroblasts_adventitial",
            "Fibroblasts_cortex",
            "Medullary_fibroblasts",
            "Pericytes",
            "MC",
            "Myofibroblasts",
        ],
    ),
    "S2": (
        ["Fibroblast", "Stroma_progenitor", "Myofibroblast"],
        [
            "Fibroblasts_NEGR1+",
            "Fibroblasts_INHBA+",
            "Fibroblasts_adventitial",
            "Fibroblasts_cortex",
            "Medullary_fibroblasts",
            "Pericytes",
            "MC",
            "Myofibroblasts",
        ],
    ),
}

# Thresholds
CONF_THRESHOLD = 0.50
PRESENCE_MIN_CELLS = 100
PRESENCE_MIN_FRAC = 0.001  # 0.1%
REF_MIN_CELLS = 50

MIN_CELLS_PB = 10
MIN_COUNTS_PB = 1000
MIN_REPLICATES = 3
DEG_FDR = 0.05
DEG_LFC = 0.5

# Publication figure style (no bold, large readable fonts, dpi=300)
FIG_RC = {
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "legend.title_fontsize": 13,
    "figure.titlesize": 18,
    "font.weight": "normal",
    "axes.titleweight": "normal",
    "axes.labelweight": "normal",
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def ensure_dirs() -> None:
    for p in DIRS.values():
        p.mkdir(parents=True, exist_ok=True)


def apply_fig_style() -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update(FIG_RC)


def merge_protocol(proto: str) -> str:
    p = str(proto).strip()
    if p in ("Morizane", "Morizane_modified"):
        return "Morizane"
    if p in ("Takasato", "Takasato_modified"):
        return "Takasato"
    if p == "Freedman":
        return "Freedman"
    return "Other"


def _dec(x) -> str:
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode()
    return str(x)


def read_h5_obs_col(f: h5py.File, key: str) -> np.ndarray:
    g = f["obs"][key]
    if isinstance(g, h5py.Group) and "categories" in g:
        cats = np.array([_dec(x) for x in g["categories"][:]], dtype=object)
        codes = np.asarray(g["codes"][:])
        return np.array([cats[c] if int(c) >= 0 else "NA" for c in codes], dtype=object)
    arr = g[:]
    if arr.dtype.kind in "SUO":
        return np.array([_dec(x) for x in arr], dtype=object)
    return np.asarray(arr)


def parse_age_numeric(val) -> float:
    """Parse organoid age ('d26' -> 26.0) or adult age_range ('50-59' -> 54.5, '72' -> 72.0)."""
    if val is None or (isinstance(val, float) and not np.isfinite(val)):
        return np.nan
    s = str(val).strip().lower()
    if s in {"", "nan", "none", "na"}:
        return np.nan
    m_range = re.match(r"^(\d+(?:\.\d+)?)\s*[-–]\s*(\d+(?:\.\d+)?)$", s)
    if m_range:
        return 0.5 * (float(m_range.group(1)) + float(m_range.group(2)))
    m_num = re.search(r"(\d+(?:\.\d+)?)", s)
    return float(m_num.group(1)) if m_num else np.nan


def load_human_tfs() -> list[str]:
    df = pd.read_csv(PATHS["tf_csv"])
    col = "gene" if "gene" in df.columns else df.columns[0]
    return sorted(df[col].dropna().astype(str).str.strip().unique().tolist())


def build_adult_fine_to_broad() -> dict[str, str]:
    """Map each adult fine_annotation to its dominant non-NA broad_annotation."""
    with h5py.File(PATHS["adult_umap"], "r") as f1, h5py.File(
        PATHS["adult_meta_broad"], "r"
    ) as f2:
        fine = read_h5_obs_col(f1, "fine_annotation")
        broad = read_h5_obs_col(f2, "broad_annotation")
    df = pd.DataFrame({"fine": fine, "broad": broad})
    non_na = df[df["broad"] != "NA"]
    mapping = (
        non_na.groupby("fine")["broad"]
        .agg(lambda s: s.value_counts().index[0])
        .to_dict()
    )
    # Explicit fallbacks for the 4 fine labels whose cells had broad_annotation == 'NA'
    fallbacks = {
        "aCD-PC": "PC",
        "Mono-Mφ": "Monocytes",
        "dEC-cap": "EC",
        "PAC": "Papillary epithelium",
    }
    for k, v in fallbacks.items():
        mapping.setdefault(k, v)
    return mapping


def build_fetal_author_to_broad() -> dict[str, str]:
    """1-to-1 map from fetal author_annotation to harmonized_celltype."""
    with h5py.File(PATHS["fetal_umap"], "r") as f:
        auth = read_h5_obs_col(f, "author_annotation")
        harm = read_h5_obs_col(f, "harmonized_celltype")
    df = pd.DataFrame({"auth": auth, "harm": harm})
    return (
        df.groupby("auth")["harm"]
        .agg(lambda s: s.value_counts().index[0])
        .to_dict()
    )
