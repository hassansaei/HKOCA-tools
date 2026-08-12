"""Packaged HKOCA config files (markers, colors, module defaults)."""

from __future__ import annotations

from importlib import resources
from pathlib import Path


def _config_path(name: str) -> Path:
    return Path(resources.files("hkoca.config").joinpath(name)).resolve()


def celltype_colors_path() -> Path:
    return _config_path("celltype_colors.yaml")


def annotation_config_path() -> Path:
    return _config_path("annotation.config.yaml")


def snapseed_markers_path() -> Path:
    return _config_path("snapseed_markers_v4.yaml")
