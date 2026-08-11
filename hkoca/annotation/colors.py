"""Fixed cell-type color palettes for Snapseed Level_1 / Level_2 / Level_3 UMAPs."""

from __future__ import annotations

FALLBACK_COLOR = "#BBBBBB"

LEVEL1_COLOR_MAP: dict[str, str] = {
    "Nephron": "#0072B2",
    "Ureteric_Epithelium": "#6A51A3",
    "Stromal_cells": "#D55E00",
    "Endothelial_cells": "#009E73",
    "Off_target_cells": "#737373",
    "undifferentiated_pluripotent_cells": "#F0E442",
}

LEVEL2_COLOR_MAP: dict[str, str] = {
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
    "Off_target_cells": "#737373",
    "undifferentiated_pluripotent_cells": "#F0E442",
}

LEVEL3_COLOR_MAP: dict[str, str] = {
    "Precursor_podocytes": "#1F4E99",
    "Mature_podocytes": "#438CE3",
    "PT (PEC-like)": "#146A8F",
    "Dev_PT": "#1D91C0",
    "PT": "#29AADB",
    "Dev_DT": "#0086B3",
    "DT": "#33CCFF",
    "Endodermal_cells": "#B49C1B",
    "Muscle_progenitor_cells": "#A59430",
    "Neural_progenitor_cells": "#F05270",
    "Glial_cells": "#CFED24",
    "Collecting_duct": "#544082",
    "Nephron_progenitor_cells": "#081D58",
    "LOH": "#008080",
    "Proliferative_tubule": "#41B6C4",
    "Ureteric_Bud": "#3F007D",
    "Stromal_progenitors": "#7A0177",
    "Proliferative_stroma": "#C51B8A",
    "Mesangial_cells": "#F768A1",
    "Pericytes": "#E31A1C",
    "Interstitial_Fibroblasts": "#FF7F00",
    "Activated Fibroblasts": "#FDBF6F",
    "Early_endothelium": "#00441B",
    "Late_endothelium": "#238B45",
    "Off_target_cells": "#252525",
    "undifferentiated_pluripotent_cells": "#F0E442",
}

# Level_latest can hold labels from any hierarchy depth.
_COMBINED_COLOR_MAP: dict[str, str] = {
    **LEVEL1_COLOR_MAP,
    **LEVEL2_COLOR_MAP,
    **LEVEL3_COLOR_MAP,
}


def _normalize_label(label: str) -> str:
    return str(label).strip()


def lookup_color(label: str, palette: dict[str, str]) -> str:
    """Resolve a label to a hex color; try space/underscore variants."""
    key = _normalize_label(label)
    if key in palette:
        return palette[key]
    underscored = key.replace(" ", "_")
    if underscored in palette:
        return palette[underscored]
    spaced = key.replace("_", " ")
    if spaced in palette:
        return palette[spaced]
    return FALLBACK_COLOR


def palette_for_obs_key(color_key: str) -> dict[str, str] | None:
    """Return the fixed palette for a Snapseed obs column, or None for Leiden/etc."""
    key = str(color_key)
    if key.startswith("Level_1_") or key == "Level_1":
        return LEVEL1_COLOR_MAP
    if key.startswith("Level_2_") or key == "Level_2":
        return LEVEL2_COLOR_MAP
    if key.startswith("Level_3_") or key == "Level_3":
        return LEVEL3_COLOR_MAP
    if key.startswith("Level_latest_") or key in {"Level_latest", "Level_3_latest"}:
        return _COMBINED_COLOR_MAP
    return None
