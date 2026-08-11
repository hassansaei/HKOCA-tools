"""Run HKOCA pipeline stages in order with resume support."""

from __future__ import annotations

import logging
import os

from hkoca.pipeline.checkpoints import (
    cellbender_complete,
    harmonize_complete,
    qc_stage_complete,
)
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


def _cellbender_suffix(cfg: PipelineConfig) -> str:
    from hkoca.cellbender.config import default_config_path, resolve_config as resolve_cellbender_config

    return resolve_cellbender_config(
        default_config_path(cfg.cellbender_config or None)
    ).output_suffix


def run_cellbender_stage(
    cfg: PipelineConfig,
    df,
    *,
    dry_run: bool = False,
    resume: bool = True,
    force: bool = False,
) -> tuple[int, set[str]]:
    from hkoca.cellbender.config import default_config_path, resolve_config as resolve_cellbender_config
    from hkoca.cellbender.runner import run_jobs

    if not cfg.run_cellbender:
        logger.info("CellBender stage skipped (--skip-cellbender or run_cellbender=false).")
        return 0, set()

    suffix = _cellbender_suffix(cfg)
    done, missing, done_ids = cellbender_complete(cfg, df, output_suffix=suffix)

    if resume and not force and done:
        logger.info(
            "CellBender stage already complete (%d sample output(s)); skipping.",
            len(done_ids),
        )
        return 0, done_ids

    samples: list[str] = []
    for _, row in df.iterrows():
        if not sample_runs_cellbender(row, run_cellbender_stage=True):
            continue
        sample_dir = resolve_sample_dir(row, cfg.working_dir)
        samples.append(sample_dir)

    if not samples:
        logger.info("No samples marked for CellBender; skipping stage.")
        return 0, set()

    if missing and resume and not force:
        logger.info(
            "CellBender resume: %d sample(s) still missing; will skip existing outputs.",
            len(missing),
        )

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
        skip_existing=resume and not force,
    )
    if rc != 0:
        return rc, set()

    _, _, done_ids = cellbender_complete(cfg, df, output_suffix=suffix)
    return 0, done_ids


def run_qc_filter_stage(
    cfg: PipelineConfig,
    df,
    harmonize_csv: str,
    *,
    dry_run: bool = False,
    verbose: bool = False,
    resume: bool = True,
    force: bool = False,
) -> int:
    from hkoca.qc_filter import main as qc_filter_main

    if resume and not force:
        ok, reason = qc_stage_complete(cfg, df)
        if ok:
            logger.info("QC-filter stage already complete (%s); skipping.", reason)
            return 0
        logger.info("QC-filter resume: %s", reason)

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
    if force:
        argv.append("--force-overwrite")
    # qc-filter skips existing by default; do not pass a redundant flag

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
        logger.info("[dry-run] would run Snapseed annotation via hkoca annotation run")
        return 0
    logger.warning(
        "Pipeline annotation stage not auto-wired yet; run "
        "`hkoca annotation run --input <h5ad> --output-dir %s` "
        "(default resolutions 0.4,0.6,1.0).",
        out_dir,
    )
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
    resume: bool = True,
    force: bool = False,
    verbose: bool = False,
) -> int:
    stages = stages_in_range(cfg.from_stage, cfg.to_stage)
    logger.info("Pipeline stages: %s", " -> ".join(stages))
    logger.info("Sample info CSV: %s", cfg.sample_info_csv)
    logger.info("Output root: %s", cfg.output_root)
    logger.info("Resume: %s | force: %s", resume and not force, force)

    df = load_sample_info(cfg.sample_info_csv)
    work_dir = os.path.join(cfg.output_root, ".hkoca")
    os.makedirs(work_dir, exist_ok=True)
    harmonize_csv = os.path.join(work_dir, "harmonize_metadata.csv")
    cellbender_sample_ids: set[str] = set()
    suffix = _cellbender_suffix(cfg)

    total = len(stages)
    for idx, stage in enumerate(stages, start=1):
        logger.info("=" * 64)
        logger.info("STAGE %d/%d: %s", idx, total, stage)

        if stage == "cellbender":
            rc, cellbender_sample_ids = run_cellbender_stage(
                cfg,
                df,
                dry_run=dry_run,
                resume=resume,
                force=force,
            )
            if rc != 0:
                logger.error("CellBender stage failed; stopping pipeline.")
                return rc
            continue

        if stage == "qc_filter":
            # If CellBender was not in this run range, still prefer filtered H5
            # paths when those outputs already exist.
            if not cellbender_sample_ids and cfg.run_cellbender:
                _, _, cellbender_sample_ids = cellbender_complete(
                    cfg, df, output_suffix=suffix
                )

            if resume and not force:
                harm_ok, harm_missing = harmonize_complete(cfg.output_root, df)
                if harm_ok:
                    logger.info("Harmonize outputs already present under %s/rds/", cfg.output_root)
                else:
                    logger.info(
                        "Harmonize resume: %d study output(s) missing (e.g. %s)",
                        len(harm_missing),
                        harm_missing[0],
                    )

            write_harmonize_csv(
                df,
                harmonize_csv,
                working_dir=cfg.working_dir,
                output_suffix=suffix,
                cellbender_sample_ids=cellbender_sample_ids,
            )
            logger.info("Harmonize metadata written: %s", harmonize_csv)
            rc = run_qc_filter_stage(
                cfg,
                df,
                harmonize_csv,
                dry_run=dry_run,
                verbose=verbose,
                resume=resume,
                force=force,
            )
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
