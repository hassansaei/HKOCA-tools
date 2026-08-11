"""Detect completed pipeline outputs so finished steps can be skipped."""

from __future__ import annotations

import glob
import os
from typing import Any

from hkoca.pipeline.config import PipelineConfig, qc_output_dir, resolve_sample_dir, sample_runs_cellbender


def _nonempty_file(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _study_names(df) -> list[str]:
    return sorted({str(row["study"]).strip() for _, row in df.iterrows() if str(row.get("study", "")).strip()})


def harmonize_rds_path(output_root: str, study: str) -> str:
    return os.path.join(output_root, "rds", f"{study}_harmonized.rds")


def harmonize_complete(output_root: str, df, *, require_rds: bool = True) -> tuple[bool, list[str]]:
    """Return (all_done, missing_paths)."""
    missing: list[str] = []
    for study in _study_names(df):
        if require_rds:
            path = harmonize_rds_path(output_root, study)
        else:
            path = os.path.join(output_root, "harmonized", f"{study}_harmonized.h5ad")
        if not _nonempty_file(path):
            missing.append(path)
    return (len(missing) == 0, missing)


def cellbender_output_path(sample_dir: str, sample_id: str, output_suffix: str) -> str:
    return os.path.join(sample_dir, f"{sample_id}{output_suffix}")


def cellbender_complete(
    cfg: PipelineConfig,
    df,
    *,
    output_suffix: str,
) -> tuple[bool, list[str], set[str]]:
    """Return (all_done, missing_paths, sample_ids_with_outputs)."""
    missing: list[str] = []
    done_ids: set[str] = set()
    for _, row in df.iterrows():
        if not sample_runs_cellbender(row, run_cellbender_stage=cfg.run_cellbender):
            continue
        sample_id = str(row["sample_id"]).strip()
        sample_dir = resolve_sample_dir(row, cfg.working_dir)
        path = cellbender_output_path(sample_dir, sample_id, output_suffix)
        if _nonempty_file(path):
            done_ids.add(sample_id)
        else:
            missing.append(path)
    return (len(missing) == 0, missing, done_ids)


def _glob_nonempty(pattern: str) -> list[str]:
    return [p for p in glob.glob(pattern) if _nonempty_file(p)]


def qc_doublet_complete(qc_output: str, df) -> tuple[bool, list[str]]:
    """Doublet stage done if each study has a matching *_singlets.rds."""
    dbl_dir = os.path.join(qc_output, "doublet_filtered_rds")
    missing: list[str] = []
    for study in _study_names(df):
        matches = _glob_nonempty(os.path.join(dbl_dir, f"{study}*_singlets.rds"))
        if not matches:
            missing.append(os.path.join(dbl_dir, f"{study}*_singlets.rds"))
    return (len(missing) == 0, missing)


def qc_filter_complete(qc_output: str, df) -> tuple[bool, list[str]]:
    """QC filtering done if each study has a matching *_filtered.rds under qc_filtered_rds."""
    filt_dir = os.path.join(qc_output, "qc_filtered_rds")
    missing: list[str] = []
    for study in _study_names(df):
        matches = _glob_nonempty(os.path.join(filt_dir, f"{study}*_filtered.rds"))
        if not matches:
            missing.append(os.path.join(filt_dir, f"{study}*_filtered.rds"))
    return (len(missing) == 0, missing)


def qc_stage_complete(cfg: PipelineConfig, df) -> tuple[bool, str]:
    """Whether the configured QC run (all/qc/doublet) looks finished."""
    qc_out = qc_output_dir(cfg)
    run = (cfg.qc_run or "all").lower()
    reasons: list[str] = []

    if run in ("all", "doublet"):
        ok, missing = qc_doublet_complete(qc_out, df)
        if not ok:
            reasons.append(f"doublet missing: {missing[0]}")
    if run in ("all", "qc"):
        ok, missing = qc_filter_complete(qc_out, df)
        if not ok:
            reasons.append(f"qc-filter missing: {missing[0]}")

    if reasons:
        return False, "; ".join(reasons)
    return True, "QC outputs present"
