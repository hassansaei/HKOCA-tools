"""Pipeline configuration and sample_info CSV helpers."""

from __future__ import annotations

import configparser
import os
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any, Literal

CONFIG_FILENAME = "pipeline.config"

PipelineStage = Literal["cellbender", "qc_filter", "annotation", "integration"]
CellBenderMode = Literal["h5", "mtx"]

PIPELINE_ONLY_COLS = {"sample_dir", "run_cellbender", "cellbender_mode"}


@dataclass(frozen=True)
class PipelineConfig:
    working_dir: str
    sample_info_csv: str
    gtf_file: str
    output_root: str
    run_cellbender: bool
    cellbender_mode: CellBenderMode
    from_stage: PipelineStage
    to_stage: PipelineStage
    cellbender_config: str
    harmonize_config: str
    transgenes: tuple[str, ...]
    qc_config: str
    qc_output: str
    qc_output_subdir: str
    qc_run: str
    annotation_output_subdir: str
    annotation_config: str
    annotation_markers: str
    integration_output_subdir: str
    integration_config: str
    integration_methods: str
    config_path: str


def default_config_path(cli_config: str | None = None) -> str:
    if cli_config and cli_config.strip():
        return cli_config.strip()
    cwd = os.path.join(os.getcwd(), CONFIG_FILENAME)
    if os.path.isfile(cwd):
        return cwd
    pkg_cfg = os.path.join(os.path.dirname(os.path.abspath(__file__)), CONFIG_FILENAME)
    if os.path.isfile(pkg_cfg):
        return pkg_cfg
    return cwd


def template_csv_path() -> Path:
    return Path(resources.files("hkoca.pipeline").joinpath("sample_info_template.csv")).resolve()


def _as_bool(value: str | bool | None, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip()
    if not text:
        return default
    return text.lower() in {"1", "true", "t", "yes", "y"}


def _parse_list(raw: str | None) -> tuple[str, ...]:
    if not raw or not str(raw).strip():
        return ()
    return tuple(s.strip() for s in str(raw).split(",") if s.strip())


def load_config_file(config_path: str) -> configparser.ConfigParser:
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    if os.path.isfile(config_path):
        cfg.read(config_path, encoding="utf-8")
    return cfg


def _get(cfg: configparser.ConfigParser, section: str, key: str, default: str = "") -> str:
    if cfg.has_option(section, key):
        return cfg.get(section, key).strip()
    return default


def resolve_config(
    config_path: str,
    *,
    working_dir: str | None = None,
    sample_info_csv: str | None = None,
    gtf_file: str | None = None,
    output_root: str | None = None,
    run_cellbender: bool | None = None,
    cellbender_mode: str | None = None,
    from_stage: str | None = None,
    to_stage: str | None = None,
    cellbender_config: str | None = None,
    harmonize_config: str | None = None,
    transgenes: str | None = None,
    qc_config: str | None = None,
    qc_output: str | None = None,
    qc_output_subdir: str | None = None,
    qc_run: str | None = None,
) -> PipelineConfig:
    cfg = load_config_file(config_path)

    def path_val(cli_val: str | None, section: str, key: str) -> str:
        if cli_val is not None and str(cli_val).strip():
            return str(cli_val).strip()
        return _get(cfg, section, key)

    wd = path_val(working_dir, "paths", "working_dir") or os.getcwd()
    if not os.path.isabs(wd):
        wd = os.path.abspath(wd)

    def abs_path(value: str) -> str:
        if not value:
            return ""
        if os.path.isabs(value):
            return os.path.normpath(value)
        return os.path.normpath(os.path.join(wd, value))

    cb_mode = (cellbender_mode or _get(cfg, "stages", "cellbender_mode", "h5") or "h5").lower()
    if cb_mode not in ("h5", "mtx"):
        raise ValueError(f"Invalid cellbender_mode: {cb_mode} (choose h5 or mtx)")

    from_st = (from_stage or _get(cfg, "stages", "from_stage", "cellbender") or "cellbender").lower()
    to_st = (to_stage or _get(cfg, "stages", "to_stage", "integration") or "integration").lower()

    run_cb = _as_bool(_get(cfg, "stages", "run_cellbender", "true"), True)
    if run_cellbender is not None:
        run_cb = run_cellbender

    if not run_cb and from_st == "cellbender":
        from_st = "qc_filter"

    return PipelineConfig(
        working_dir=wd,
        sample_info_csv=abs_path(path_val(sample_info_csv, "paths", "sample_info_csv")),
        gtf_file=abs_path(path_val(gtf_file, "paths", "gtf_file")),
        output_root=abs_path(path_val(output_root, "paths", "output_root")),
        run_cellbender=run_cb,
        cellbender_mode=cb_mode,  # type: ignore[arg-type]
        from_stage=from_st,  # type: ignore[arg-type]
        to_stage=to_st,  # type: ignore[arg-type]
        cellbender_config=path_val(cellbender_config, "cellbender", "config"),
        harmonize_config=path_val(harmonize_config, "harmonize", "config"),
        transgenes=_parse_list(transgenes or _get(cfg, "harmonize", "transgenes")),
        qc_config=path_val(qc_config, "qc", "config"),
        qc_output=abs_path(path_val(qc_output, "qc", "output_dir")),
        qc_output_subdir=path_val(qc_output_subdir, "qc", "output_subdir") or "qc_filter",
        qc_run=path_val(qc_run, "qc", "run") or "all",
        annotation_output_subdir=_get(cfg, "annotation", "output_subdir", "annotation"),
        annotation_config=path_val(None, "annotation", "config"),
        annotation_markers=path_val(None, "annotation", "markers"),
        integration_output_subdir=_get(cfg, "integration", "output_subdir", "integration"),
        integration_config=path_val(None, "integration", "config"),
        integration_methods=_get(cfg, "integration", "methods", "harmony,rpca,cca") or "harmony,rpca,cca",
        config_path=config_path,
    )


def validate_config(cfg: PipelineConfig) -> None:
    missing: list[str] = []
    if not cfg.sample_info_csv:
        missing.append("sample_info_csv (--csv)")
    if not cfg.gtf_file:
        missing.append("gtf_file (--gtf)")
    if not cfg.output_root:
        missing.append("output_root (--output)")
    if missing:
        raise ValueError("Missing required pipeline inputs: " + ", ".join(missing))
    if not os.path.isfile(cfg.sample_info_csv):
        raise FileNotFoundError(f"Sample info CSV not found: {cfg.sample_info_csv}")
    if not os.path.isfile(cfg.gtf_file):
        raise FileNotFoundError(f"GTF file not found: {cfg.gtf_file}")


def load_sample_info(csv_path: str) -> Any:
    from hkoca.qc_filter.harmonize.pipeline import load_metadata_csv

    return load_metadata_csv(csv_path)


def resolve_sample_dir(row: Any, working_dir: str) -> str:
    explicit = str(row.get("sample_dir", "") or "").strip()
    if explicit:
        path = explicit
    else:
        data_path = str(row["data_path"]).strip()
        if os.path.isfile(data_path) or data_path.lower().endswith((".h5", ".h5ad", ".hdf5")):
            path = os.path.dirname(data_path)
        else:
            path = data_path
    if working_dir and not os.path.isabs(path):
        path = os.path.normpath(os.path.join(working_dir, path))
    return path


def sample_runs_cellbender(row: Any, *, run_cellbender_stage: bool) -> bool:
    if not run_cellbender_stage:
        return False
    return _as_bool(row.get("run_cellbender", "true"), True)


def harmonize_data_path(
    row: Any,
    *,
    working_dir: str,
    output_suffix: str,
    used_cellbender: bool,
) -> str:
    if used_cellbender:
        sample_dir = resolve_sample_dir(row, working_dir)
        sample_id = str(row["sample_id"]).strip()
        return os.path.join(sample_dir, f"{sample_id}{output_suffix}")
    path = str(row["data_path"]).strip()
    if working_dir and not os.path.isabs(path):
        path = os.path.normpath(os.path.join(working_dir, path))
    return path


def write_harmonize_csv(
    df: Any,
    dest_path: str,
    *,
    working_dir: str,
    output_suffix: str,
    cellbender_sample_ids: set[str],
) -> str:
    import pandas as pd

    from hkoca.qc_filter.harmonize.pipeline import MANDATORY_META_COLS

    drop_cols = PIPELINE_ONLY_COLS | {"skip", "output_dir"}
    rows: list[dict] = []
    for _, row in df.iterrows():
        record = {col: row[col] for col in df.columns if col not in drop_cols}
        sample_id = str(row["sample_id"]).strip()
        record["data_path"] = harmonize_data_path(
            row,
            working_dir=working_dir,
            output_suffix=output_suffix,
            used_cellbender=sample_id in cellbender_sample_ids,
        )
        if "file_prefix" in df.columns:
            record["file_prefix"] = row.get("file_prefix", "")
        rows.append(record)

    out_df = pd.DataFrame(rows)
    for col in MANDATORY_META_COLS:
        if col not in out_df.columns:
            out_df[col] = ""
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    out_df.to_csv(dest_path, index=False)
    return dest_path


def qc_output_dir(cfg: PipelineConfig) -> str:
    if cfg.qc_output:
        return cfg.qc_output
    return os.path.join(cfg.output_root, cfg.qc_output_subdir)


def stages_in_range(from_stage: PipelineStage, to_stage: PipelineStage) -> tuple[PipelineStage, ...]:
    order = ("cellbender", "qc_filter", "annotation", "integration")
    start = order.index(from_stage)
    end = order.index(to_stage)
    if start > end:
        raise ValueError(f"from_stage ({from_stage}) must come before to_stage ({to_stage})")
    return order[start : end + 1]
