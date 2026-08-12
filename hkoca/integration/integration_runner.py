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
from hkoca.env_check import log_env_problems, verify_stage_env

logger = logging.getLogger("hkoca.integration")

DEFAULT_INTEGRATION_ENV = "hkoca_integration"
DEFAULT_ANNOTATION_ENV = "hkoca_harmonize"


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


def resolve_annotation_python() -> str | None:
    env_name = os.environ.get("HKOCA_ANNOTATION_ENV", DEFAULT_ANNOTATION_ENV).strip()
    if not env_name:
        env_name = DEFAULT_ANNOTATION_ENV
    prefix = resolve_env_prefix(env_name, "python")
    if prefix is None:
        return None
    py = prefix / "bin" / "python"
    return str(py) if py.is_file() else None


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

    ann_py = resolve_annotation_python()
    if ann_py:
        env["HKOCA_ANNOTATION_PYTHON"] = ann_py
        env["RETICULATE_PYTHON"] = ann_py
        logger.info("Snapseed re-annotation Python: %s", ann_py)
    else:
        logger.warning(
            "Annotation env not found (%s). Snapseed re-annotation will be skipped in integration run.",
            os.environ.get("HKOCA_ANNOTATION_ENV", DEFAULT_ANNOTATION_ENV),
        )

    return env


def _verify_integration_env(stage: str = "integration", *, require_annotation: bool = False) -> int:
    problems = verify_stage_env(stage)
    if stage == "run":
        ann_env = os.environ.get("HKOCA_ANNOTATION_ENV", DEFAULT_ANNOTATION_ENV).strip() or DEFAULT_ANNOTATION_ENV
        ann_prefix = resolve_env_prefix(ann_env, "python")
        if ann_prefix is None:
            msg = (
                f"annotation: conda env '{ann_env}' not found. "
                "Snapseed re-annotation requires hkoca_harmonize (or set HKOCA_ANNOTATION_ENV)."
            )
            if require_annotation:
                problems.append(msg)
            else:
                logger.warning(msg)
        else:
            ann_problems = verify_stage_env("harmonize")
            for msg in ann_problems:
                if require_annotation:
                    problems.append(f"annotation: {msg}")
                else:
                    logger.warning("annotation: %s", msg)
    if problems:
        log_env_problems(stage, problems)
        return 1
    return 0


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
    remove_intermediate: bool = False,
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
    if remove_intermediate:
        cmd.append("--remove_intermediate")
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
    remove_intermediate: bool = False,
    extra_args: list[str] | None = None,
) -> list[str]:
    from hkoca.config import snapseed_markers_path

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
    if remove_intermediate:
        cmd.append("--remove_intermediate")
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
    remove_intermediate: bool = False,
    dry_run: bool = False,
    extra_args: list[str] | None = None,
) -> int:
    if not dry_run and not os.path.isfile(prepared_rds):
        logger.error("Prepared RDS does not exist: %s", prepared_rds)
        return 1

    if not dry_run and _verify_integration_env("run") != 0:
        return 1

    try:
        cmd = build_run_command(
            prepared_rds=prepared_rds,
            output_dir=output_dir,
            methods=methods,
            config=config,
            force_overwrite=force_overwrite,
            remove_intermediate=remove_intermediate,
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
    remove_intermediate: bool = False,
    dry_run: bool = False,
    extra_args: list[str] | None = None,
) -> int:
    if not dry_run and not os.path.isfile(input_rds):
        logger.error("Input RDS does not exist: %s", input_rds)
        return 1
    if annotated_h5ad and not dry_run and not os.path.isfile(annotated_h5ad):
        logger.error("Annotated h5ad does not exist: %s", annotated_h5ad)
        return 1

    if not dry_run and _verify_integration_env("integration") != 0:
        return 1

    try:
        cmd = build_prep_command(
            input_rds=input_rds,
            output_dir=output_dir,
            annotated_h5ad=annotated_h5ad,
            config=config,
            force_overwrite=force_overwrite,
            remove_intermediate=remove_intermediate,
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
