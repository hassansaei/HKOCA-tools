"""Load projection module config."""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Any

import yaml


def packaged_config_path() -> Path:
    return Path(resources.files("hkoca.projection").joinpath("projection.config.yaml")).resolve()


def packaged_prototype_map_path() -> Path:
    return Path(
        resources.files("hkoca.projection.data").joinpath(
            "atlas_prototype_to_Level3Integrated_map.tsv"
        )
    ).resolve()


def load_yaml(path: Path | str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config must be a mapping: {path}")
    return data


def resolve_path(value: Any, base: Path) -> Path | None:
    if value is None or (isinstance(value, str) and not str(value).strip()):
        return None
    p = Path(str(value)).expanduser()
    if not p.is_absolute():
        p = (base / p).resolve()
    return p


def load_projection_config(
    config_path: Path | str | None = None,
    *,
    working_dir: Path | str | None = None,
) -> dict[str, Any]:
    base = Path(working_dir).expanduser().resolve() if working_dir else Path.cwd()
    raw: dict[str, Any] = {}
    if config_path:
        cfg_file = Path(config_path).expanduser().resolve()
        if cfg_file.is_file():
            raw = load_yaml(cfg_file)
            base = cfg_file.parent
        elif not Path(config_path).is_absolute():
            packaged = packaged_config_path()
            if packaged.is_file():
                raw = load_yaml(packaged)
    else:
        packaged = packaged_config_path()
        if packaged.is_file():
            raw = load_yaml(packaged)

    paths = raw.get("paths") or {}
    reference = raw.get("reference") or {}
    query = raw.get("query") or {}
    scpoli = raw.get("scpoli") or {}
    outputs = raw.get("outputs") or {}
    plotting = raw.get("plotting") or {}
    prototype_map = resolve_path(reference.get("prototype_map"), base) or packaged_prototype_map_path()

    return {
        "config_path": str(Path(config_path).resolve()) if config_path else str(packaged_config_path()),
        "atlas_h5ad": resolve_path(paths.get("atlas_h5ad"), base),
        "query_path": resolve_path(paths.get("query") or paths.get("query_h5ad"), base),
        "model_dir": resolve_path(paths.get("model_dir"), base),
        "output_dir": resolve_path(paths.get("output_dir"), base) or (base / "projection_results").resolve(),
        "query_batch_key": query.get("batch_key") or "sample_id",
        "query_label_column": query.get("label_column"),
        "counts_layer": query.get("counts_layer") or "counts",
        "celltype_key": reference.get("celltype_key") or "Level_3_Integrated",
        "condition_key": reference.get("condition_key") or "sample_id",
        "unknown_label": reference.get("unknown_label") or "Unknown",
        "atlas_umap_key": reference.get("atlas_umap_key") or "X_umap_scpoli",
        "atlas_latent_key": reference.get("atlas_latent_key") or "X_scpoli",
        "prototype_map": prototype_map,
        "n_epochs": int(scpoli.get("n_epochs") or 50),
        "pretrain_epochs": int(scpoli.get("pretrain_epochs") or 40),
        "eta": float(scpoli.get("eta") or 10),
        "seed": int(scpoli.get("seed") or 42),
        "projected_subdir": outputs.get("projected_subdir") or "projected_obj",
        "figures_subdir": outputs.get("figures_subdir") or "figures",
        "tables_subdir": outputs.get("tables_subdir") or "tables",
        "logs_subdir": outputs.get("logs_subdir") or "logs",
        "models_subdir": outputs.get("models_subdir") or "models",
        "save_plots": bool(plotting.get("save_plots", True)),
        "dpi": int(plotting.get("dpi") or 150),
        "skip_existing": bool(plotting.get("skip_existing", True)),
        "joint_umap": bool(plotting.get("joint_umap", False)),
        "atlas_umap_subsample": int(plotting.get("atlas_umap_subsample") or 80000),
        "knn_neighbors": int(plotting.get("knn_neighbors") or 30),
        "atlas_bg_max": int(plotting.get("atlas_bg_max") or 120000),
    }
