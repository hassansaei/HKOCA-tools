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
DEFAULT_ANNOTATION_ENV = "hkoca_harmonize"
DEFAULT_METHODS = ("harmony", "rpca", "cca")
DEFAULT_TRANSGENES = ("AAV", "EGFP", "mCherry", "GFP", "eGFP", "LK03_eGFP")


def r_script_path(stage: str = "prep") -> Path:
    scripts = {
        "prep": "integration_prep.R",
        "run": "integration_methods.R",
        "export_benchmark": "export_benchmark.R",
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
    env["HKOCA_TRANSGENES"] = os.environ.get("HKOCA_TRANSGENES", "").strip() or ",".join(
        DEFAULT_TRANSGENES
    )
    return env


def parse_methods(methods: str | None) -> list[str]:
    raw = methods or ",".join(DEFAULT_METHODS)
    out = [m.strip().lower() for m in raw.split(",") if m.strip()]
    return [m for m in out if m in DEFAULT_METHODS] or list(DEFAULT_METHODS)


def method_rds_path(output_dir: str, method: str) -> Path:
    return Path(output_dir) / "objects" / f"integrated_{method}.rds"


def all_method_rds_exist(output_dir: str, methods: list[str]) -> bool:
    for method in methods:
        path = method_rds_path(output_dir, method)
        if not path.is_file() or path.stat().st_size <= 0:
            return False
    return True


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
    from hkoca.integration.benchmark import benchmark_complete

    methods_list = parse_methods(methods)
    if not dry_run and not os.path.isfile(prepared_rds):
        logger.error("Prepared RDS does not exist: %s", prepared_rds)
        return 1

    try:
        cmd = build_run_command(
            prepared_rds=prepared_rds,
            output_dir=output_dir,
            methods=",".join(methods_list),
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
    rc = result.returncode
    if rc != 0:
        logger.error("Integration methods failed (exit %s).", rc)
        return rc
    logger.info("Integration methods completed successfully. Outputs under: %s", output_dir)

    if dry_run:
        return 0

    if not all_method_rds_exist(output_dir, methods_list):
        logger.error("Not all requested method RDS exist; skipping scIB benchmark.")
        return 1

    if not force_overwrite and benchmark_complete(output_dir):
        logger.info("Skipping scIB benchmark; results already exist under %s/benchmark.", output_dir)
        return 0

    return run_scib_stage(
        prepared_rds=prepared_rds,
        output_dir=output_dir,
        methods=methods_list,
        force_overwrite=force_overwrite,
    )


def run_scib_stage(
    *,
    prepared_rds: str,
    output_dir: str,
    methods: list[str],
    force_overwrite: bool = False,
) -> int:
    from hkoca.integration.benchmark import scib_install_hint

    try:
        rscript = find_rscript()
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1

    py = resolve_annotation_python() or shutil.which("python")
    if not py:
        logger.warning(
            "No Python interpreter found for scIB benchmark; skipping. "
            "Set HKOCA_ANNOTATION_ENV to enable benchmark metrics."
        )
        return 0

    env = os.environ.copy()
    prefix = Path(py).resolve().parent.parent
    if (prefix / "bin" / "python").is_file():
        env = subprocess_env_for_prefix(prefix)

    probe = subprocess.run(
        [py, "-c", "import scanpy, scib"],
        check=False,
        env=env,
        capture_output=True,
        text=True,
    )
    if probe.returncode != 0:
        env_name = os.environ.get("HKOCA_ANNOTATION_ENV", DEFAULT_ANNOTATION_ENV).strip() or DEFAULT_ANNOTATION_ENV
        logger.warning(
            "Skipping optional scIB benchmark (%s). Integration RDS outputs are unchanged.",
            scib_install_hint(env_name),
        )
        return 0

    export_cmd = [
        rscript,
        str(r_script_path("export_benchmark")),
        "--output_dir",
        os.path.abspath(output_dir),
        "--prepared_rds",
        os.path.abspath(prepared_rds),
        "--methods",
        ",".join(methods),
    ]
    if force_overwrite:
        export_cmd.append("--force_overwrite")

    logger.info("Exporting embeddings for scIB: %s", " ".join(export_cmd))
    result = subprocess.run(export_cmd, check=False, env=_subprocess_env_for_rscript(rscript))
    if result.returncode != 0:
        logger.error("Benchmark embedding export failed (exit %s).", result.returncode)
        return result.returncode

    bench_cmd = [
        py,
        "-m",
        "hkoca.integration.benchmark",
        "--output-dir",
        os.path.abspath(output_dir),
        "--methods",
        ",".join(methods),
    ]
    logger.info("Running scIB benchmark: %s", " ".join(bench_cmd))
    result = subprocess.run(bench_cmd, check=False, env=env)
    if result.returncode != 0:
        logger.warning(
            "Optional scIB benchmark failed (exit %s). Integration RDS outputs are unchanged. %s",
            result.returncode,
            scib_install_hint(
                os.environ.get("HKOCA_ANNOTATION_ENV", DEFAULT_ANNOTATION_ENV).strip()
                or DEFAULT_ANNOTATION_ENV
            ),
        )
        return 0
    logger.info("scIB benchmark complete: %s/benchmark", output_dir)
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
