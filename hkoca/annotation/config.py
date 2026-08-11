"""Load and resolve Snapseed annotation config / marker YAML paths."""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Any

import yaml

DEFAULT_RESOLUTIONS: tuple[float, ...] = (0.4, 0.6, 1.0)

DEFAULT_PARAMS: dict[str, Any] = {
    "seed": 42,
    "scanpy_verbosity": 1,
    "resolutions": list(DEFAULT_RESOLUTIONS),
    "hvg_n_top_genes": 3000,
    "pca_n_comps": 50,
    "neighbors_n_neighbors": 20,
    "neighbors_n_pcs": 30,
    "scale_max_value": 10,
    "normalize_target_sum": 10000,
    "save_plots": False,
    "dpi": 150,
    "skip_existing": True,
}


def packaged_markers_path() -> Path:
    return Path(
        resources.files("hkoca.annotation.data").joinpath("snapseed_markers_v4.yaml")
    ).resolve()


def packaged_config_path() -> Path:
    return Path(
        resources.files("hkoca.annotation.data").joinpath("annotation.config.yaml")
    ).resolve()


def load_yaml(path: Path | str) -> Any:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def load_marker_hierarchy(path: Path | str | None = None) -> dict:
    markers_path = Path(path) if path else packaged_markers_path()
    data = load_yaml(markers_path)
    if not isinstance(data, dict) or not data:
        raise ValueError(f"Marker YAML is empty or invalid: {markers_path}")
    # Drop comment-only / non-dict leaves used as notes.
    cleaned = {k: v for k, v in data.items() if isinstance(v, dict)}
    if not cleaned:
        raise ValueError(f"No hierarchical marker entries found in: {markers_path}")
    return cleaned


def export_markers(dest: Path | str, *, source: Path | str | None = None) -> Path:
    """Copy marker YAML to an editable location (creates parent dirs)."""
    src = Path(source) if source else packaged_markers_path()
    dest_path = Path(dest).expanduser().resolve()
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    dest_path.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
    return dest_path


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _as_float_list(value: Any) -> list[float]:
    if value is None:
        return list(DEFAULT_RESOLUTIONS)
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, str):
        parts = [p.strip() for p in value.replace(";", ",").split(",") if p.strip()]
        return [float(p) for p in parts]
    return [float(v) for v in value]


def resolve_path(value: Any, base: Path) -> Path | None:
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    p = Path(str(value)).expanduser()
    if not p.is_absolute():
        p = (base / p).resolve()
    return p


def load_annotation_config(
    config_path: Path | str | None = None,
    *,
    working_dir: Path | str | None = None,
) -> dict[str, Any]:
    """Return merged config dict with resolved paths and parameters."""
    cfg_path = Path(config_path).expanduser().resolve() if config_path else None
    raw: dict[str, Any] = {}
    if cfg_path and cfg_path.is_file():
        loaded = load_yaml(cfg_path) or {}
        if not isinstance(loaded, dict):
            raise ValueError(f"Config must be a mapping: {cfg_path}")
        raw = loaded
        base = cfg_path.parent
    else:
        base = Path(working_dir).expanduser().resolve() if working_dir else Path.cwd()

    if working_dir:
        base = Path(working_dir).expanduser().resolve()

    paths_raw = dict(raw.get("paths") or {})
    params = {**DEFAULT_PARAMS, **dict(raw.get("parameters") or {})}
    params["resolutions"] = _as_float_list(params.get("resolutions"))
    params["save_plots"] = _as_bool(params.get("save_plots"), False)
    params["skip_existing"] = _as_bool(params.get("skip_existing"), True)

    markers = resolve_path(paths_raw.get("markers"), base) or packaged_markers_path()
    output_dir = resolve_path(paths_raw.get("output_dir"), base) or (base / "annotation_results")
    annotated_subdir = str(paths_raw.get("annotated_subdir") or "annotated")
    clustered_subdir = str(paths_raw.get("clustered_subdir") or "clustered")
    figures_subdir = str(paths_raw.get("figures_subdir") or "figures")

    return {
        "config_path": cfg_path,
        "base_dir": base,
        "input": resolve_path(paths_raw.get("input"), base),
        "markers": markers,
        "output_dir": output_dir,
        "annotated_dir": output_dir / annotated_subdir,
        "clustered_dir": output_dir / clustered_subdir,
        "figures_dir": output_dir / figures_subdir,
        "parameters": params,
    }
