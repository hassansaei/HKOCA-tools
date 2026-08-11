"""Run HKOCA pipeline stages in order."""

from __future__ import annotations

import logging
import os

from hkoca.pipeline.config import (
    PipelineConfig,
    load_sample_info,
    qc_output_dir,
    resolve_sample_dir,
    sample_runs_cellbender,
    stages_in_range,
    write_harmonize_csv,
)

logger = logging.getLogger("hkoca.pipeline")


def run_cellbender_stage(
    cfg: PipelineConfig,
    df,
    *,
    dry_run: bool = False,
    skip_existing: bool = False,
) -> tuple[int, set[str]]:
    from hkoca.cellbender.config import default_config_path, resolve_config as resolve_cellbender_config
    from hkoca.cellbender.runner import run_jobs

    if not cfg.run_cellbender:
        logger.info("CellBender stage skipped (--skip-cellbender or run_cellbender=false).")
        return 0, set()

    samples: list[str] = []
    for _, row in df.iterrows():
        if not sample_runs_cellbender(row, run_cellbender_stage=True):
            continue
        sample_dir = resolve_sample_dir(row, cfg.working_dir)
        samples.append(sample_dir)

    if not samples:
        logger.info("No samples marked for CellBender; skipping stage.")
        return 0, set()

    cb_cfg = resolve_cellbender_config(
        default_config_path(cfg.cellbender_config or None),
        working_dir=cfg.working_dir,
        samples_dir=None,
        samples=tuple(dict.fromkeys(samples)),
    )
    rc = run_jobs(
        cb_cfg,
        cfg.cellbender_mode,
        dry_run=dry_run,
        skip_existing=skip_existing,
    )
    ran_ids = {str(row["sample_id"]).strip() for _, row in df.iterrows() if sample_runs_cellbender(row, run_cellbender_stage=True)}
    return rc, ran_ids if rc == 0 else set()


def run_qc_filter_stage(
    cfg: PipelineConfig,
    harmonize_csv: str,
    *,
    dry_run: bool = False,
    verbose: bool = False,
) -> int:
    from hkoca.qc_filter import main as qc_filter_main

    qc_output = qc_output_dir(cfg)
    argv = [
        "--qc-output",
        qc_output,
        "--run",
        cfg.qc_run,
    ]
    if cfg.qc_config:
        argv.extend(["--qc-config", cfg.qc_config])
    if dry_run:
        argv.append("--dry-run")
    if verbose:
        argv.append("-v")

    argv.extend(
        [
            "--csv",
            harmonize_csv,
            "--gtf",
            cfg.gtf_file,
            "--output",
            cfg.output_root,
            "-w",
            cfg.working_dir,
        ]
    )
    if cfg.harmonize_config:
        argv.extend(["--config", cfg.harmonize_config])
    if cfg.transgenes:
        argv.extend(["--transgenes", ",".join(cfg.transgenes)])

    logger.info("QC-filter argv: %s", " ".join(argv))
    return int(qc_filter_main(argv) or 0)


def run_annotation_stage(cfg: PipelineConfig, *, dry_run: bool = False) -> int:
    from hkoca.annotation import status_message

    out_dir = os.path.join(cfg.output_root, cfg.annotation_output_subdir)
    logger.info("STEP: annotation -> %s", out_dir)
    if dry_run:
        logger.info("[dry-run] would run annotation (not implemented yet)")
        return 0
    logger.warning("Annotation stage is not implemented yet; skipping.")
    for line in status_message().splitlines():
        logger.info("  %s", line)
    return 0


def run_integration_stage(cfg: PipelineConfig, *, dry_run: bool = False) -> int:
    from hkoca.integration import status_message

    out_dir = os.path.join(cfg.output_root, cfg.integration_output_subdir)
    logger.info("STEP: integration -> %s", out_dir)
    if dry_run:
        logger.info("[dry-run] would run integration (not implemented yet)")
        return 0
    logger.warning("Integration stage is not implemented yet; skipping.")
    for line in status_message().splitlines():
        logger.info("  %s", line)
    return 0


def run_pipeline(
    cfg: PipelineConfig,
    *,
    dry_run: bool = False,
    skip_existing: bool = False,
    verbose: bool = False,
) -> int:
    stages = stages_in_range(cfg.from_stage, cfg.to_stage)
    logger.info("Pipeline stages: %s", " -> ".join(stages))
    logger.info("Sample info CSV: %s", cfg.sample_info_csv)
    logger.info("Output root: %s", cfg.output_root)

    df = load_sample_info(cfg.sample_info_csv)
    work_dir = os.path.join(cfg.output_root, ".hkoca")
    harmonize_csv = os.path.join(work_dir, "harmonize_metadata.csv")
    cellbender_sample_ids: set[str] = set()

    total = len(stages)
    for idx, stage in enumerate(stages, start=1):
        logger.info("=" * 64)
        logger.info("STAGE %d/%d: %s", idx, total, stage)

        if stage == "cellbender":
            rc, cellbender_sample_ids = run_cellbender_stage(
                cfg,
                df,
                dry_run=dry_run,
                skip_existing=skip_existing,
            )
            if rc != 0:
                logger.error("CellBender stage failed; stopping pipeline.")
                return rc
            continue

        if stage == "qc_filter":
            from hkoca.cellbender.config import default_config_path, resolve_config as resolve_cellbender_config

            cb_cfg = resolve_cellbender_config(default_config_path(cfg.cellbender_config or None))
            write_harmonize_csv(
                df,
                harmonize_csv,
                working_dir=cfg.working_dir,
                output_suffix=cb_cfg.output_suffix,
                cellbender_sample_ids=cellbender_sample_ids,
            )
            logger.info("Harmonize metadata written: %s", harmonize_csv)
            rc = run_qc_filter_stage(cfg, harmonize_csv, dry_run=dry_run, verbose=verbose)
            if rc != 0:
                logger.error("QC-filter stage failed; stopping pipeline.")
                return rc
            continue

        if stage == "annotation":
            rc = run_annotation_stage(cfg, dry_run=dry_run)
            if rc != 0:
                return rc
            continue

        if stage == "integration":
            rc = run_integration_stage(cfg, dry_run=dry_run)
            if rc != 0:
                return rc
            continue

    logger.info("=" * 64)
    logger.info("Pipeline complete. Outputs under: %s", cfg.output_root)
    return 0
