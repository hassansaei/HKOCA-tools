"""Invoke the bundled QC_scdbl_Combined.R engine.

QC R packages (scDblFinder, r-anndata, ...) live in the ``sc_qc_pipeline``
conda env (``conda/environment_qc.yaml``). Harmonization usually runs in
``hkoca_harmonize``, so this module resolves the QC env's ``Rscript`` and
runs the subprocess with that env on ``PATH`` / ``CONDA_PREFIX`` — no manual
``conda activate`` mid-pipeline.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from importlib import resources
from pathlib import Path

from hkoca.conda_env import resolve_env_prefix, subprocess_env_for_prefix

logger = logging.getLogger("hkoca.qc_filter")

# Matches ``name:`` in conda/environment_qc.yaml and docker/Dockerfile.
DEFAULT_QC_ENV = "sc_qc_pipeline"


def r_script_path() -> Path:
    return Path(resources.files("hkoca.qc_filter.r").joinpath("QC_scdbl_Combined.R")).resolve()


def packaged_qc_config() -> Path:
    return Path(resources.files("hkoca.qc_filter.r").joinpath("qc_config.dcf")).resolve()


def resolve_qc_env_prefix() -> Path | None:
    """Return the QC conda env prefix if it exists on this machine."""
    env_name = os.environ.get("HKOCA_QC_ENV", DEFAULT_QC_ENV).strip() or DEFAULT_QC_ENV
    return resolve_env_prefix(env_name, "Rscript")


def find_rscript() -> str:
    """Locate Rscript for QC, preferring ``sc_qc_pipeline`` over PATH."""
    override = os.environ.get("HKOCA_RSCRIPT", "").strip()
    if override:
        if not Path(override).is_file():
            raise FileNotFoundError(f"HKOCA_RSCRIPT is set but not a file: {override}")
        return str(Path(override).resolve())

    qc_prefix = resolve_qc_env_prefix()
    if qc_prefix is not None:
        rscript = qc_prefix / "bin" / "Rscript"
        logger.info("Using QC conda env Rscript: %s", rscript)
        return str(rscript)

    exe = shutil.which("Rscript")
    if exe is None:
        env_name = os.environ.get("HKOCA_QC_ENV", DEFAULT_QC_ENV)
        raise FileNotFoundError(
            "Rscript not found. Create the QC conda env and ensure it is "
            f"discoverable (expected env name '{env_name}' from "
            "conda/environment_qc.yaml), or set HKOCA_RSCRIPT to an Rscript "
            "that provides scDblFinder + anndata."
        )

    logger.warning(
        "QC env '%s' not found; using Rscript on PATH (%s). "
        "If QC fails with missing packages, create conda/environment_qc.yaml "
        "or set HKOCA_QC_ENV / HKOCA_RSCRIPT.",
        os.environ.get("HKOCA_QC_ENV", DEFAULT_QC_ENV),
        exe,
    )
    return exe


def _subprocess_env_for_rscript(rscript: str) -> dict[str, str]:
    """Run R under the same conda prefix as ``rscript`` (reticulate + libs)."""
    prefix = Path(rscript).resolve().parent.parent
    return subprocess_env_for_prefix(prefix)


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
    result = subprocess.run(
        cmd,
        check=False,
        env=_subprocess_env_for_rscript(cmd[0]),
    )
    if result.returncode != 0:
        logger.error("QC R script failed (exit %s)", result.returncode)
        return result.returncode
    logger.info("QC stage completed. Outputs under: %s", output_dir)
    return 0
