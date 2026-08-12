"""Fixed cell-type color palettes for Snapseed Level_1 / Level_2 / Level_3 plots."""

from __future__ import annotations

from functools import lru_cache

import yaml

from hkoca.config import celltype_colors_path

FALLBACK_COLOR = "#BBBBBB"


@lru_cache(maxsize=1)
def _load_palettes_yaml() -> dict:
    with open(celltype_colors_path(), encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}
    return raw


def _palette_from_yaml(key: str) -> dict[str, str]:
    raw = _load_palettes_yaml()
    block = raw.get(key) or {}
    return {str(k): str(v) for k, v in block.items()}


@lru_cache(maxsize=1)
def level1_color_map() -> dict[str, str]:
    return _palette_from_yaml("level1")


@lru_cache(maxsize=1)
def level2_color_map() -> dict[str, str]:
    return _palette_from_yaml("level2")


@lru_cache(maxsize=1)
def level3_color_map() -> dict[str, str]:
    return _palette_from_yaml("level3")


@lru_cache(maxsize=1)
def combined_color_map() -> dict[str, str]:
    merged: dict[str, str] = {}
    merged.update(level1_color_map())
    merged.update(level2_color_map())
    merged.update(level3_color_map())
    return merged


def fallback_color() -> str:
    raw = _load_palettes_yaml()
    return str(raw.get("fallback") or FALLBACK_COLOR)


# Backward-compatible module-level names (lazy via functions at import time).
LEVEL1_COLOR_MAP = level1_color_map()
LEVEL2_COLOR_MAP = level2_color_map()
LEVEL3_COLOR_MAP = level3_color_map()
_COMBINED_COLOR_MAP = combined_color_map()


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
    return fallback_color()


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
