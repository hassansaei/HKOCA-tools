"""Build and execute CellBender remove-background commands."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import Literal

from hkoca.cellbender.config import (
    CellBenderConfig,
    resolve_sample_dir,
    sample_id_from_path,
)

logger = logging.getLogger("hkoca.cellbender")

InputMode = Literal["h5", "mtx"]


@dataclass(frozen=True)
class SampleJob:
    sample_id: str
    sample_dir: str
    input_path: str
    output_path: str


def find_cellbender() -> str:
    exe = shutil.which("cellbender")
    if exe is None:
        raise FileNotFoundError(
            "cellbender executable not found on PATH. "
            "Activate the hkoca_cellbender conda env "
            "(conda/environment_cellbender.yaml)."
        )
    return exe


def build_jobs(cfg: CellBenderConfig, mode: InputMode) -> list[SampleJob]:
    if not cfg.samples:
        raise ValueError(
            "No samples specified. Set [paths] samples in the config "
            "or pass --samples / --samples-file."
        )

    jobs: list[SampleJob] = []
    for sample in cfg.samples:
        sample_dir = resolve_sample_dir(cfg, sample)
        sid = sample_id_from_path(sample_dir)
        if mode == "h5":
            input_path = os.path.join(sample_dir, cfg.h5_filename)
        else:
            input_path = os.path.join(sample_dir, cfg.mtx_dirname)
        output_path = os.path.join(sample_dir, f"{sid}{cfg.output_suffix}")
        jobs.append(
            SampleJob(
                sample_id=sid,
                sample_dir=sample_dir,
                input_path=input_path,
                output_path=output_path,
            )
        )
    return jobs


def build_command(exe: str, job: SampleJob, cfg: CellBenderConfig) -> list[str]:
    p = cfg.params
    cmd = [
        exe,
        "remove-background",
        "--input",
        job.input_path,
        "--output",
        job.output_path,
        "--mode",
        p.mode,
        "--epochs",
        str(p.epochs),
        "--posterior-batch-size",
        str(p.posterior_batch_size),
        "--total-droplets-included",
        str(p.total_droplets_included),
        "--learning-rate",
        str(p.learning_rate),
        "--cpu-threads",
        str(p.cpu_threads),
    ]
    if p.cuda:
        cmd.append("--cuda")
    return cmd


def validate_job(job: SampleJob, mode: InputMode) -> None:
    if not os.path.isdir(job.sample_dir):
        raise FileNotFoundError(f"Sample directory not found: {job.sample_dir}")
    if mode == "h5":
        if not os.path.isfile(job.input_path):
            raise FileNotFoundError(f"H5 input not found: {job.input_path}")
        return

    if not os.path.isdir(job.input_path):
        raise FileNotFoundError(f"MTX directory not found: {job.input_path}")

    names = {f.lower() for f in os.listdir(job.input_path)}
    has_matrix = any(n.startswith("matrix.mtx") for n in names)
    has_barcodes = any(n.startswith("barcodes.tsv") for n in names)
    has_features = any(
        n.startswith("features.tsv") or n.startswith("genes.tsv") for n in names
    )
    if not (has_matrix and has_barcodes and has_features):
        raise FileNotFoundError(
            f"MTX triplet incomplete in {job.input_path}. "
            "Expected matrix.mtx[.gz], barcodes.tsv[.gz], and "
            "features.tsv[.gz] or genes.tsv[.gz]."
        )


def run_jobs(
    cfg: CellBenderConfig,
    mode: InputMode,
    *,
    dry_run: bool = False,
    skip_existing: bool = False,
) -> int:
    exe = "cellbender" if dry_run else find_cellbender()
    jobs = build_jobs(cfg, mode)
    failed = 0

    logger.info("CellBender executable: %s", exe)
    logger.info("Config: %s", cfg.config_path)
    logger.info("Mode: %s | samples: %d", mode, len(jobs))
    logger.info(
        "Params: epochs=%s droplets=%s lr=%s threads=%s cuda=%s",
        cfg.params.epochs,
        cfg.params.total_droplets_included,
        cfg.params.learning_rate,
        cfg.params.cpu_threads,
        cfg.params.cuda,
    )

    for job in jobs:
        try:
            validate_job(job, mode)
        except FileNotFoundError as exc:
            logger.error("[%s] %s", job.sample_id, exc)
            failed += 1
            continue

        if skip_existing and os.path.isfile(job.output_path) and os.path.getsize(job.output_path) > 0:
            logger.info("[%s] skip existing output: %s", job.sample_id, job.output_path)
            continue

        cmd = build_command(exe, job, cfg)
        logger.info("[%s] %s", job.sample_id, " ".join(cmd))
        if dry_run:
            continue

        os.makedirs(job.sample_dir, exist_ok=True)
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            logger.error("[%s] cellbender failed (exit %s)", job.sample_id, result.returncode)
            failed += 1
        else:
            logger.info("[%s] wrote %s", job.sample_id, job.output_path)

    if failed:
        logger.error("Finished with %d failed sample(s).", failed)
        return 1
    logger.info("All samples completed successfully.")
    return 0
