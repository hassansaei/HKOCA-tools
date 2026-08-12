"""Detect completed pipeline outputs so finished steps can be skipped."""

from __future__ import annotations

import glob
import os

from hkoca.pipeline.config import PipelineConfig, qc_output_dir, resolve_sample_dir, sample_runs_cellbender
from hkoca.pipeline.paths import (
    annotation_output_dir,
    discover_study_qc_filtered_h5ad,
    expected_annotated_h5ad,
    integration_method_rds,
    integration_prepared_rds,
)


def _nonempty_file(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _glob_nonempty(pattern: str) -> list[str]:
    return sorted(p for p in glob.glob(pattern) if _nonempty_file(p))


def _study_names(df) -> list[str]:
    return sorted({str(row["study"]).strip() for _, row in df.iterrows() if str(row.get("study", "")).strip()})


def _parse_methods(methods: str) -> tuple[str, ...]:
    return tuple(m.strip().lower() for m in methods.split(",") if m.strip())


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


def annotation_complete(cfg: PipelineConfig, df) -> tuple[bool, str]:
    """Annotation done when each study has a matching *_annotated.h5ad."""
    qc_out = qc_output_dir(cfg)
    ann_root = annotation_output_dir(cfg)
    missing: list[str] = []

    for study in _study_names(df):
        filtered_h5ad = discover_study_qc_filtered_h5ad(qc_out, study)
        if not filtered_h5ad:
            missing.append(f"{study}: QC filtered h5ad missing under h5ad_converted/")
            continue
        annotated = expected_annotated_h5ad(ann_root, filtered_h5ad)
        if not _nonempty_file(annotated):
            missing.append(annotated)

    if missing:
        return False, missing[0]
    return True, "annotation outputs present"


def integration_prep_complete(integration_dir: str) -> bool:
    return _nonempty_file(integration_prepared_rds(integration_dir))


def integration_run_complete(integration_dir: str, methods: str) -> tuple[bool, list[str]]:
    missing: list[str] = []
    for method in _parse_methods(methods):
        path = integration_method_rds(integration_dir, method)
        if not _nonempty_file(path):
            missing.append(path)
    return (len(missing) == 0, missing)


def integration_stage_complete(
    cfg: PipelineConfig,
    df,
    *,
    methods: str,
) -> tuple[bool, str]:
    from hkoca.pipeline.paths import integration_output_dir

    studies = _study_names(df)
    n_studies = len(studies)
    missing: list[str] = []

    for study in studies:
        int_dir = integration_output_dir(cfg, study, n_studies=n_studies)
        if not integration_prep_complete(int_dir):
            missing.append(integration_prepared_rds(int_dir))
            continue
        ok, method_missing = integration_run_complete(int_dir, methods)
        if not ok:
            missing.extend(method_missing)

    if missing:
        return False, missing[0]
    return True, "integration outputs present"
