"""Resolve sibling conda env prefixes for stage subprocesses."""

from __future__ import annotations

import os
import shutil
from pathlib import Path


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
