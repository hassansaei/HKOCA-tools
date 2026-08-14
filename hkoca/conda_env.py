"""Resolve sibling conda env prefixes for stage subprocesses."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

# Names must match ``name:`` in conda/environment_*.yaml
ENV_CELLBENDER = "hkoca_cellbender"
ENV_HARMONIZE = "hkoca_harmonize"
ENV_QC = "sc_qc_pipeline"
ENV_INTEGRATION = "hkoca_integration"
ENV_PROJECTION = "hkoca_projection"


def conda_roots() -> list[Path]:
    """Candidate conda install roots (base prefixes that contain ``envs/``)."""
    roots: list[Path] = []
    seen: set[str] = set()

    def add(path: Path | None) -> None:
        if path is None:
            return
        try:
            resolved = path.resolve()
        except OSError:
            return
        key = str(resolved)
        if key in seen:
            return
        seen.add(key)
        roots.append(resolved)

    explicit = os.environ.get("CONDA_ROOT") or os.environ.get("MAMBA_ROOT_PREFIX")
    if explicit:
        add(Path(explicit))

    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        p = Path(prefix)
        if p.parent.name == "envs":
            add(p.parent.parent)
        else:
            add(p)

    conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if conda_exe:
        add(Path(conda_exe).resolve().parent.parent)

    add(Path("/opt/conda"))
    return roots


def resolve_env_prefix(env_name: str, executable: str) -> Path | None:
    """Return env prefix if ``envs/<env_name>/bin/<executable>`` exists."""
    current = os.environ.get("CONDA_PREFIX")
    if current:
        cur = Path(current)
        if cur.name == env_name and (cur / "bin" / executable).is_file():
            return cur.resolve()

    for root in conda_roots():
        candidate = root / "envs" / env_name
        if (candidate / "bin" / executable).is_file():
            return candidate.resolve()
    return None


def subprocess_env_for_prefix(prefix: Path) -> dict[str, str]:
    """Run a subprocess with ``prefix/bin`` on PATH and CONDA_PREFIX set."""
    env = os.environ.copy()
    bin_dir = prefix / "bin"
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["CONDA_PREFIX"] = str(prefix)
    env.pop("PYTHONHOME", None)

    py = bin_dir / "python"
    if py.is_file():
        env["RETICULATE_PYTHON"] = str(py)
    return env


def resolve_python(env_name: str) -> str | None:
    prefix = resolve_env_prefix(env_name, "python")
    if prefix is None:
        return None
    py = prefix / "bin" / "python"
    return str(py) if py.is_file() else None


def wire_harmonize_python(env: dict[str, str]) -> dict[str, str]:
    """Point reticulate / integration prep at the harmonize env Python (scanpy/anndata)."""
    ann_env = os.environ.get("HKOCA_ANNOTATION_ENV", ENV_HARMONIZE).strip() or ENV_HARMONIZE
    ann_py = resolve_python(ann_env)
    if ann_py:
        env["HKOCA_ANNOTATION_PYTHON"] = ann_py
        env["RETICULATE_PYTHON"] = ann_py
    return env


def r_env_for_rscript(rscript: str, *, harmonize_python: bool = False) -> dict[str, str]:
    """Build subprocess env for an Rscript under its conda prefix."""
    prefix = Path(rscript).resolve().parent.parent
    env = subprocess_env_for_prefix(prefix)
    if harmonize_python:
        env = wire_harmonize_python(env)
    return env


def python_env(env_name: str, *, need_hkoca: bool = False) -> tuple[str, dict[str, str]]:
    """Return (python, env) for a named conda env, optionally wiring hkoca on PYTHONPATH."""
    prefix = resolve_env_prefix(env_name, "python")
    if prefix is None:
        raise FileNotFoundError(
            f"Conda env '{env_name}' not found. Create it from conda/environment_*.yaml."
        )
    python = str(prefix / "bin" / "python")
    env = subprocess_env_for_prefix(prefix)
    if need_hkoca:
        env = ensure_hkoca_on_pythonpath(env, python)
    return python, env


def hkoca_source_root() -> Path:
    """Directory that contains the ``hkoca`` package (repo root for editable installs)."""
    import hkoca

    return Path(hkoca.__file__).resolve().parent.parent


def ensure_hkoca_on_pythonpath(env: dict[str, str], python: str) -> dict[str, str]:
    """Ensure ``hkoca`` is importable in a stage-specific conda env subprocess."""
    probe = subprocess.run(
        [python, "-c", "import hkoca"],
        env=env,
        capture_output=True,
        text=True,
    )
    if probe.returncode == 0:
        return env
    root = str(hkoca_source_root())
    existing = env.get("PYTHONPATH", "").strip()
    env["PYTHONPATH"] = root if not existing else f"{root}{os.pathsep}{existing}"
    return env
