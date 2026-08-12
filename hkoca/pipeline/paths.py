"""Resolve stage input/output paths for the HKOCA pipeline."""

from __future__ import annotations

import glob
import os
from pathlib import Path

from hkoca.pipeline.config import PipelineConfig, qc_output_dir


def _nonempty_file(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _glob_nonempty(pattern: str) -> list[str]:
    return sorted(p for p in glob.glob(pattern) if _nonempty_file(p))


def annotation_output_dir(cfg: PipelineConfig) -> str:
    return os.path.join(cfg.output_root, cfg.annotation_output_subdir)


def integration_base_dir(cfg: PipelineConfig) -> str:
    return os.path.join(cfg.output_root, cfg.integration_output_subdir)


def integration_output_dir(cfg: PipelineConfig, study: str, *, n_studies: int) -> str:
    base = integration_base_dir(cfg)
    if n_studies <= 1:
        return base
    return os.path.join(base, study)


def qc_h5ad_dir(qc_output: str) -> str:
    return os.path.join(qc_output, "h5ad_converted")


def qc_filtered_rds_dir(qc_output: str) -> str:
    return os.path.join(qc_output, "qc_filtered_rds")


def discover_study_qc_filtered_h5ad(qc_output: str, study: str) -> str | None:
    matches = _glob_nonempty(os.path.join(qc_h5ad_dir(qc_output), f"{study}*_filtered.h5ad"))
    return matches[0] if matches else None


def discover_study_qc_filtered_rds(qc_output: str, study: str) -> str | None:
    matches = _glob_nonempty(os.path.join(qc_filtered_rds_dir(qc_output), f"{study}*_filtered.rds"))
    return matches[0] if matches else None


def expected_annotated_h5ad(annotation_root: str, qc_h5ad_path: str) -> str:
    stem = Path(qc_h5ad_path).stem
    return os.path.join(annotation_root, "annotated_obj", f"{stem}_annotated.h5ad")


def integration_prepared_rds(integration_dir: str) -> str:
    return os.path.join(integration_dir, "prep", "sct_prepared.rds")


def integration_method_rds(integration_dir: str, method: str) -> str:
    return os.path.join(integration_dir, "objects", f"integrated_{method}.rds")


def collect_study_artifacts(cfg: PipelineConfig, df) -> list[dict[str, str]]:
    from hkoca.pipeline.checkpoints import _study_names

    qc_out = qc_output_dir(cfg)
    ann_root = annotation_output_dir(cfg)
    studies = _study_names(df)
    n_studies = len(studies)
    rows: list[dict[str, str]] = []

    for study in studies:
        filtered_h5ad = discover_study_qc_filtered_h5ad(qc_out, study)
        filtered_rds = discover_study_qc_filtered_rds(qc_out, study)
        if not filtered_rds:
            raise FileNotFoundError(
                f"Missing QC filtered RDS for study '{study}' under {qc_filtered_rds_dir(qc_out)}"
            )
        annotated_h5ad = ""
        if filtered_h5ad:
            ann_path = expected_annotated_h5ad(ann_root, filtered_h5ad)
            if _nonempty_file(ann_path):
                annotated_h5ad = ann_path
        rows.append(
            {
                "study": study,
                "filtered_h5ad": filtered_h5ad or "",
                "filtered_rds": filtered_rds,
                "annotated_h5ad": annotated_h5ad,
                "annotation_root": ann_root,
                "integration_dir": integration_output_dir(cfg, study, n_studies=n_studies),
            }
        )
    return rows
