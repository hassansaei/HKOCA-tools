"""Run HKOCA pipeline stages in order with resume support."""

from __future__ import annotations

import logging
import os

from hkoca.pipeline.checkpoints import (
    annotation_complete,
    cellbender_complete,
    harmonize_complete,
    integration_prep_complete,
    integration_run_complete,
    integration_stage_complete,
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
from hkoca.pipeline.paths import (
    annotation_output_dir,
    collect_study_artifacts,
    discover_study_qc_filtered_h5ad,
    expected_annotated_h5ad,
    integration_prepared_rds,
    qc_h5ad_dir,
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


def run_annotation_stage(
    cfg: PipelineConfig,
    df,
    *,
    dry_run: bool = False,
    resume: bool = True,
    force: bool = False,
) -> int:
    from hkoca.annotation.config import DEFAULT_RESOLUTIONS, load_annotation_config
    from hkoca.annotation.runner import run_annotation_batch
    from hkoca.config import snapseed_markers_path

    ann_root = annotation_output_dir(cfg)
    qc_out = qc_output_dir(cfg)

    if resume and not force:
        ok, reason = annotation_complete(cfg, df)
        if ok:
            logger.info("Annotation stage already complete (%s); skipping.", reason)
            return 0
        logger.info("Annotation resume: %s", reason)

    inputs: list[str] = []
    for study in sorted({str(row["study"]).strip() for _, row in df.iterrows()}):
        filtered_h5ad = discover_study_qc_filtered_h5ad(qc_out, study)
        if not filtered_h5ad:
            logger.error(
                "Missing QC filtered h5ad for study '%s' under %s. "
                "Ensure qc-filter ran with h5ad conversion enabled.",
                study,
                qc_h5ad_dir(qc_out),
            )
            return 1
        inputs.append(filtered_h5ad)

    logger.info("Annotation inputs (%d):", len(inputs))
    for path in inputs:
        logger.info("  %s -> %s", path, expected_annotated_h5ad(ann_root, path))

    if dry_run:
        logger.info("[dry-run] would annotate %d h5ad file(s) under %s", len(inputs), ann_root)
        return 0

    ann_cfg = load_annotation_config(cfg.annotation_config or None, working_dir=cfg.working_dir)
    params = dict(ann_cfg["parameters"])
    if not params.get("resolutions"):
        params["resolutions"] = list(DEFAULT_RESOLUTIONS)
    if force:
        params["skip_existing"] = False

    markers = cfg.annotation_markers or str(snapseed_markers_path())
    annotated_dir = os.path.join(ann_root, "annotated_obj")
    clustered_dir = os.path.join(ann_root, "clustered")
    figures_dir = os.path.join(ann_root, "figures")
    os.makedirs(annotated_dir, exist_ok=True)

    run_annotation_batch(
        inputs,
        markers_path=markers,
        resolutions=params["resolutions"],
        annotated_dir=annotated_dir,
        clustered_dir=clustered_dir,
        figures_dir=figures_dir,
        parameters=params,
        force=force,
    )
    ok, reason = annotation_complete(cfg, df)
    if not ok:
        logger.error("Annotation stage finished but outputs are incomplete: %s", reason)
        return 1
    logger.info("Annotation stage completed successfully under %s", ann_root)
    return 0


def run_integration_stage(
    cfg: PipelineConfig,
    df,
    *,
    dry_run: bool = False,
    resume: bool = True,
    force: bool = False,
) -> int:
    from hkoca.integration.integration_runner import run_methods, run_prep

    if resume and not force:
        ok, reason = integration_stage_complete(cfg, df, methods=cfg.integration_methods)
        if ok:
            logger.info("Integration stage already complete (%s); skipping.", reason)
            return 0
        logger.info("Integration resume: %s", reason)

    try:
        artifacts = collect_study_artifacts(cfg, df)
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1

    for row in artifacts:
        study = row["study"]
        int_dir = row["integration_dir"]
        filtered_rds = row["filtered_rds"]
        annotated_h5ad = row["annotated_h5ad"] or None
        prepared_rds = integration_prepared_rds(int_dir)

        logger.info("Integration study: %s", study)
        logger.info("  filtered RDS : %s", filtered_rds)
        logger.info("  annotated    : %s", annotated_h5ad or "(none)")
        logger.info("  output dir   : %s", int_dir)

        if dry_run:
            logger.info("[dry-run] would run integration prep + run for study %s", study)
            continue

        if resume and not force and integration_prep_complete(int_dir):
            logger.info("Integration prep already complete for %s; skipping prep.", study)
        else:
            extra = (
                ["--transgenes", ",".join(cfg.transgenes)] if cfg.transgenes else None
            )
            rc = run_prep(
                input_rds=filtered_rds,
                output_dir=int_dir,
                annotated_h5ad=annotated_h5ad,
                config=cfg.integration_config or None,
                force_overwrite=force,
                dry_run=False,
                extra_args=extra,
            )
            if rc != 0:
                logger.error("Integration prep failed for study '%s'.", study)
                return rc

        if resume and not force:
            ok, missing = integration_run_complete(int_dir, cfg.integration_methods)
            if ok:
                logger.info("Integration run already complete for %s; skipping run.", study)
                continue
            logger.info(
                "Integration run resume for %s: missing %d method output(s).",
                study,
                len(missing),
            )

        extra = (
            ["--transgenes", ",".join(cfg.transgenes)] if cfg.transgenes else None
        )
        rc = run_methods(
            prepared_rds=prepared_rds,
            output_dir=int_dir,
            methods=cfg.integration_methods,
            config=cfg.integration_config or None,
            force_overwrite=force,
            dry_run=False,
            extra_args=extra,
        )
        if rc != 0:
            logger.error("Integration run failed for study '%s'.", study)
            return rc

    if dry_run:
        return 0

    ok, reason = integration_stage_complete(cfg, df, methods=cfg.integration_methods)
    if not ok:
        logger.error("Integration stage finished but outputs are incomplete: %s", reason)
        return 1
    logger.info("Integration stage completed successfully.")
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
    logger.info(
        "Projection is a separate module; run `hkoca projection map` after integration."
    )

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
            rc = run_annotation_stage(
                cfg,
                df,
                dry_run=dry_run,
                resume=resume,
                force=force,
            )
            if rc != 0:
                logger.error("Annotation stage failed; stopping pipeline.")
                return rc
            continue

        if stage == "integration":
            rc = run_integration_stage(
                cfg,
                df,
                dry_run=dry_run,
                resume=resume,
                force=force,
            )
            if rc != 0:
                logger.error("Integration stage failed; stopping pipeline.")
                return rc
            continue

    logger.info("=" * 64)
    logger.info("Pipeline complete. Outputs under: %s", cfg.output_root)
    return 0
