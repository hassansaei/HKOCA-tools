"""Load projection module config."""

from __future__ import annotations

from importlib import resources
from pathlib import Path
from typing import Any

import yaml

DEFAULT_LABEL_PRIORITY: tuple[str, ...] = (
    "celltype_final",
    "Level_3",
    "Level_3_latest",
    "Level_2",
    "Level_1",
)


def packaged_config_path() -> Path:
    return Path(resources.files("hkoca.projection").joinpath("projection.config.yaml")).resolve()


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
    projection = raw.get("projection") or {}
    outputs = raw.get("outputs") or {}
    plotting = raw.get("plotting") or {}

    label_columns = reference.get("label_columns")
    if isinstance(label_columns, str):
        label_columns = [c.strip() for c in label_columns.split(",") if c.strip()]
    elif label_columns is not None:
        label_columns = list(label_columns)

    return {
        "config_path": str(Path(config_path).resolve()) if config_path else str(packaged_config_path()),
        "atlas_h5ad": resolve_path(paths.get("atlas_h5ad"), base),
        "query_h5ad": resolve_path(paths.get("query_h5ad"), base),
        "output_dir": resolve_path(paths.get("output_dir"), base) or (base / "projection_results").resolve(),
        "reference_label_columns": label_columns,
        "reference_batch_key": reference.get("batch_key"),
        "query_batch_key": query.get("batch_key") or "sample_id",
        "query_label_column": query.get("label_column"),
        "normalize_target_sum": float(projection.get("normalize_target_sum") or 10000),
        "hvg_n_top_genes": int(projection.get("hvg_n_top_genes") or 3000),
        "pca_n_comps": int(projection.get("pca_n_comps") or 50),
        "neighbors_n_neighbors": int(projection.get("neighbors_n_neighbors") or 20),
        "neighbors_n_pcs": int(projection.get("neighbors_n_pcs") or 30),
        "seed": int(projection.get("seed") or 42),
        "projected_subdir": outputs.get("projected_subdir") or "projected_obj",
        "figures_subdir": outputs.get("figures_subdir") or "figures",
        "tables_subdir": outputs.get("tables_subdir") or "tables",
        "logs_subdir": outputs.get("logs_subdir") or "logs",
        "save_plots": bool(plotting.get("save_plots", True)),
        "dpi": int(plotting.get("dpi") or 150),
        "skip_existing": bool(plotting.get("skip_existing", True)),
    }
