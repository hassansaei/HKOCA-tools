"""Invoke the bundled integration R scripts (Seurat / SCT prep)."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from importlib import resources
from pathlib import Path

from hkoca.conda_env import resolve_env_prefix, subprocess_env_for_prefix

from hkoca.config import celltype_colors_path

logger = logging.getLogger("hkoca.integration")

DEFAULT_INTEGRATION_ENV = "hkoca_integration"


def r_script_path(stage: str = "prep") -> Path:
    scripts = {
        "prep": "integration_prep.R",
        "run": "integration_methods.R",
    }
    name = scripts.get(stage)
    if name is None:
        raise ValueError(f"Unknown integration stage: {stage}")
    return Path(resources.files("hkoca.integration.r").joinpath(name)).resolve()


def packaged_config() -> Path:
    return Path(resources.files("hkoca.integration.r").joinpath("integration.config.dcf")).resolve()


def resolve_integration_env_prefix() -> Path | None:
    env_name = os.environ.get("HKOCA_INTEGRATION_ENV", DEFAULT_INTEGRATION_ENV).strip()
    if not env_name:
        env_name = DEFAULT_INTEGRATION_ENV
    return resolve_env_prefix(env_name, "Rscript")


def find_rscript() -> str:
    override = os.environ.get("HKOCA_RSCRIPT", "").strip()
    if override:
        path = Path(override)
        if not path.is_file():
            raise FileNotFoundError(f"HKOCA_RSCRIPT is set but not a file: {override}")
        return str(path.resolve())

    prefix = resolve_integration_env_prefix()
    if prefix is not None:
        rscript = prefix / "bin" / "Rscript"
        logger.info("Using integration conda env Rscript: %s", rscript)
        return str(rscript)

    exe = shutil.which("Rscript")
    if exe is None:
        env_name = os.environ.get("HKOCA_INTEGRATION_ENV", DEFAULT_INTEGRATION_ENV)
        raise FileNotFoundError(
            f"Rscript not found. Create conda env '{env_name}' from "
            "conda/environment_integration.yaml, or set HKOCA_INTEGRATION_ENV / HKOCA_RSCRIPT."
        )

    logger.warning(
        "Integration env '%s' not found; using Rscript on PATH (%s).",
        os.environ.get("HKOCA_INTEGRATION_ENV", DEFAULT_INTEGRATION_ENV),
        exe,
    )
    return exe


def _subprocess_env_for_rscript(rscript: str) -> dict[str, str]:
    prefix = Path(rscript).resolve().parent.parent
    env = subprocess_env_for_prefix(prefix)
    env["HKOCA_CELLTYPE_COLORS"] = str(celltype_colors_path())
    return env


def default_config_path() -> Path:
    cwd = Path.cwd() / "integration.config.dcf"
    if cwd.is_file():
        return cwd.resolve()
    return packaged_config()


def build_prep_command(
    *,
    input_rds: str,
    output_dir: str,
    annotated_h5ad: str | None = None,
    config: str | None = None,
    force_overwrite: bool = False,
    extra_args: list[str] | None = None,
) -> list[str]:
    cmd = [
        find_rscript(),
        str(r_script_path("prep")),
        "--config",
        str(config or default_config_path()),
        "--celltype_colors_yaml",
        str(celltype_colors_path()),
        "--input_rds",
        os.path.abspath(input_rds),
        "--output_dir",
        os.path.abspath(output_dir),
    ]
    if annotated_h5ad:
        cmd.extend(["--annotated_h5ad", os.path.abspath(annotated_h5ad)])
    if force_overwrite:
        cmd.append("--force_overwrite")
    if extra_args:
        cmd.extend(extra_args)
    return cmd


def build_run_command(
    *,
    prepared_rds: str,
    output_dir: str,
    methods: str | None = None,
    config: str | None = None,
    force_overwrite: bool = False,
    extra_args: list[str] | None = None,
) -> list[str]:
    from hkoca.config import celltype_colors_path, snapseed_markers_path

    cmd = [
        find_rscript(),
        str(r_script_path("run")),
        "--config",
        str(config or default_config_path()),
        "--celltype_colors_yaml",
        str(celltype_colors_path()),
        "--markers_yaml",
        str(snapseed_markers_path()),
        "--prepared_rds",
        os.path.abspath(prepared_rds),
        "--output_dir",
        os.path.abspath(output_dir),
    ]
    if methods:
        cmd.extend(["--methods", methods])
    if force_overwrite:
        cmd.append("--force_overwrite")
    if extra_args:
        cmd.extend(extra_args)
    return cmd


def run_methods(
    *,
    prepared_rds: str,
    output_dir: str,
    methods: str | None = None,
    config: str | None = None,
    force_overwrite: bool = False,
    dry_run: bool = False,
    extra_args: list[str] | None = None,
) -> int:
    if not dry_run and not os.path.isfile(prepared_rds):
        logger.error("Prepared RDS does not exist: %s", prepared_rds)
        return 1

    try:
        cmd = build_run_command(
            prepared_rds=prepared_rds,
            output_dir=output_dir,
            methods=methods,
            config=config,
            force_overwrite=force_overwrite,
            extra_args=extra_args,
        )
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1

    logger.info("Integration methods command: %s", " ".join(cmd))
    if dry_run:
        return 0

    os.makedirs(output_dir, exist_ok=True)
    result = subprocess.run(cmd, check=False, env=_subprocess_env_for_rscript(cmd[0]))
    if result.returncode != 0:
        logger.error("Integration methods failed (exit %s).", result.returncode)
        return result.returncode
    logger.info("Integration methods completed successfully. Outputs under: %s", output_dir)
    return 0


def run_prep(
    *,
    input_rds: str,
    output_dir: str,
    annotated_h5ad: str | None = None,
    config: str | None = None,
    force_overwrite: bool = False,
    dry_run: bool = False,
    extra_args: list[str] | None = None,
) -> int:
    if not dry_run and not os.path.isfile(input_rds):
        logger.error("Input RDS does not exist: %s", input_rds)
        return 1
    if annotated_h5ad and not dry_run and not os.path.isfile(annotated_h5ad):
        logger.error("Annotated h5ad does not exist: %s", annotated_h5ad)
        return 1

    try:
        cmd = build_prep_command(
            input_rds=input_rds,
            output_dir=output_dir,
            annotated_h5ad=annotated_h5ad,
            config=config,
            force_overwrite=force_overwrite,
            extra_args=extra_args,
        )
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1

    logger.info("Integration prep command: %s", " ".join(cmd))
    if dry_run:
        return 0

    os.makedirs(output_dir, exist_ok=True)
    result = subprocess.run(cmd, check=False, env=_subprocess_env_for_rscript(cmd[0]))
    if result.returncode != 0:
        logger.error("Integration prep failed (exit %s).", result.returncode)
        return result.returncode
    logger.info("Integration prep completed successfully. Outputs under: %s", output_dir)
    return 0
