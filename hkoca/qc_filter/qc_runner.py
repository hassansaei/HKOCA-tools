"""Invoke the bundled QC_scdbl_Combined.R engine."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from importlib import resources
from pathlib import Path

logger = logging.getLogger("hkoca.qc_filter")


def r_script_path() -> Path:
    return Path(resources.files("hkoca.qc_filter.r").joinpath("QC_scdbl_Combined.R")).resolve()


def packaged_qc_config() -> Path:
    return Path(resources.files("hkoca.qc_filter.r").joinpath("qc_config.dcf")).resolve()


def find_rscript() -> str:
    exe = shutil.which("Rscript")
    if exe is None:
        raise FileNotFoundError(
            "Rscript not found on PATH. Activate the QC conda env "
            "(conda/environment_qc.yaml) or an env that provides R + Seurat."
        )
    return exe


def default_qc_config() -> Path:
    cwd = Path.cwd() / "qc_config.dcf"
    if cwd.is_file():
        return cwd.resolve()
    return packaged_qc_config()


def build_qc_command(
    *,
    rds_dir: str,
    output_dir: str,
    config: str | None = None,
    stage: str = "all",
    recursive: bool = True,
    rds_pattern: str = r"_harmonized\.rds$",
    force_overwrite: bool = False,
    extra_args: list[str] | None = None,
) -> list[str]:
    script = str(r_script_path())
    cfg = config or str(default_qc_config())
    cmd = [
        find_rscript(),
        script,
        "--config",
        cfg,
        "--stage",
        stage,
        "--rds_dir",
        os.path.abspath(rds_dir),
        "--output_dir",
        os.path.abspath(output_dir),
        "--rds_pattern",
        rds_pattern,
    ]
    if recursive:
        cmd.extend(["--recursive_discovery", "TRUE"])
    if force_overwrite:
        cmd.append("--force_overwrite")
    if extra_args:
        cmd.extend(extra_args)
    return cmd


def run_qc(
    *,
    rds_dir: str,
    output_dir: str,
    config: str | None = None,
    stage: str = "all",
    recursive: bool = True,
    rds_pattern: str = r"_harmonized\.rds$",
    force_overwrite: bool = False,
    dry_run: bool = False,
    extra_args: list[str] | None = None,
) -> int:
    if not dry_run and not os.path.isdir(rds_dir):
        logger.error("RDS directory does not exist: %s", rds_dir)
        return 1

    try:
        cmd = build_qc_command(
            rds_dir=rds_dir,
            output_dir=output_dir,
            config=config,
            stage=stage,
            recursive=recursive,
            rds_pattern=rds_pattern,
            force_overwrite=force_overwrite,
            extra_args=extra_args,
        )
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1

    logger.info("QC command: %s", " ".join(cmd))
    if dry_run:
        return 0

    os.makedirs(output_dir, exist_ok=True)
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        logger.error("QC R script failed (exit %s)", result.returncode)
        return result.returncode
    logger.info("QC stage completed. Outputs under: %s", output_dir)
    return 0
