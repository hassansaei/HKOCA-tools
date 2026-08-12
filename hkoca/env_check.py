"""Verify conda environments expose required executables and imports."""

from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass

from hkoca.conda_env import resolve_env_prefix, subprocess_env_for_prefix

logger = logging.getLogger("hkoca.env_check")


@dataclass(frozen=True)
class EnvSpec:
    env_var: str
    default_name: str
    executable: str
    python_modules: tuple[str, ...] = ()
    r_packages: tuple[str, ...] = ()


ENV_SPECS: dict[str, EnvSpec] = {
    "harmonize": EnvSpec(
        env_var="HKOCA_HARMONIZE_ENV",
        default_name="hkoca_harmonize",
        executable="python",
        python_modules=("scanpy", "anndata", "snapseed", "igraph", "leidenalg"),
    ),
    "qc": EnvSpec(
        env_var="HKOCA_QC_ENV",
        default_name="sc_qc_pipeline",
        executable="Rscript",
        r_packages=(
            "Seurat", "MASS", "ggplot2", "patchwork", "cowplot", "scales",
            "reshape2", "ggridges", "ggrepel", "SingleCellExperiment",
            "scDblFinder", "yaml", "dplyr", "jsonlite", "digest", "anndata",
            "reticulate",
        ),
    ),
    "integration": EnvSpec(
        env_var="HKOCA_INTEGRATION_ENV",
        default_name="hkoca_integration",
        executable="Rscript",
        r_packages=(
            "Seurat", "ggplot2", "patchwork", "clustree", "glmGamPoi",
            "cluster", "digest", "jsonlite", "yaml", "dplyr", "tidyr",
            "harmony", "reticulate",
        ),
    ),
    "cellbender": EnvSpec(
        env_var="HKOCA_CELLBENDER_ENV",
        default_name="hkoca_cellbender",
        executable="cellbender",
        python_modules=("cellbender", "torch"),
    ),
}


def resolve_env_name(stage: str) -> str:
    spec = ENV_SPECS[stage]
    import os

    name = os.environ.get(spec.env_var, spec.default_name).strip()
    return name or spec.default_name


def verify_stage_env(stage: str, *, strict: bool = True) -> list[str]:
    """Return list of problems; empty list means OK."""
    if stage not in ENV_SPECS:
        return [f"Unknown stage: {stage}"]

    spec = ENV_SPECS[stage]
    env_name = resolve_env_name(stage)
    prefix = resolve_env_prefix(env_name, spec.executable)
    problems: list[str] = []

    if prefix is None:
        msg = (
            f"{stage}: conda env '{env_name}' not found or missing "
            f"bin/{spec.executable}. Create from conda/environment_*.yaml."
        )
        if strict:
            problems.append(msg)
        else:
            logger.warning(msg)
        return problems

    exe = prefix / "bin" / spec.executable
    env = subprocess_env_for_prefix(prefix)

    if spec.python_modules:
        mod_expr = ", ".join(repr(m) for m in spec.python_modules)
        code = (
            "import importlib, sys; missing=[]\n"
            f"for m in [{mod_expr}]:\n"
            "    try:\n"
            "        importlib.import_module(m)\n"
            "    except Exception:\n"
            "        missing.append(m)\n"
            "if missing:\n"
            "    sys.exit('missing python modules: ' + ', '.join(missing))\n"
        )
        result = subprocess.run(
            [str(exe), "-c", code],
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            detail = (result.stderr or result.stdout or "").strip()
            problems.append(f"{stage}: {detail or 'python import check failed'}")

    if spec.r_packages:
        pkg_expr = ", ".join(repr(p) for p in spec.r_packages)
        code = (
            "missing <- c()\n"
            f"for (p in c({pkg_expr})) {{\n"
            "  if (!requireNamespace(p, quietly=TRUE)) missing <- c(missing, p)\n"
            "}\n"
            "if (length(missing)) stop('missing R packages: ', paste(missing, collapse=', '))\n"
        )
        result = subprocess.run(
            [str(exe), "-e", code],
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            detail = (result.stderr or result.stdout or "").strip()
            problems.append(f"{stage}: {detail or 'R package check failed'}")

    return problems


def log_env_problems(stage: str, problems: list[str]) -> None:
    for msg in problems:
        logger.error(msg)
    if problems:
        logger.error(
            "Fix the '%s' environment (see conda/*.yaml) or set %s.",
            stage,
            ENV_SPECS[stage].env_var,
        )
