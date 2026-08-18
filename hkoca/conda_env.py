"""Resolve sibling conda env prefixes for stage subprocesses."""

from __future__ import annotations

import os
import shutil
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


def _is_site_packages(path: Path) -> bool:
    return path.name in {"site-packages", "dist-packages"}


def _looks_like_hkoca_root(path: Path) -> bool:
    return (path / "hkoca" / "__init__.py").is_file() and not _is_site_packages(path)


def _sanitize_ld_library_path(existing: str, prefix: Path) -> str:
    """Prefer ``prefix/lib`` (conda libstdc++) and drop other env lib dirs."""
    prefix_s = str(prefix.resolve())
    parts = [str(prefix / "lib")]
    for item in existing.split(os.pathsep):
        if not item:
            continue
        if item in parts:
            continue
        if "/envs/" in item.replace("\\", "/") and not item.startswith(prefix_s):
            continue
        parts.append(item)
    return os.pathsep.join(parts)


def subprocess_env_for_prefix(prefix: Path) -> dict[str, str]:
    """Run a subprocess with ``prefix/bin`` on PATH and CONDA_PREFIX set."""
    env = os.environ.copy()
    bin_dir = prefix / "bin"
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env["CONDA_PREFIX"] = str(prefix)
    env["CONDA_DEFAULT_ENV"] = prefix.name
    env.pop("PYTHONHOME", None)
    env["PYTHONNOUSERSITE"] = "1"
    env["LD_LIBRARY_PATH"] = _sanitize_ld_library_path(
        env.get("LD_LIBRARY_PATH", ""), prefix
    )
    # Do not inherit another env's site-packages via PYTHONPATH.
    env.pop("PYTHONPATH", None)

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
            f"Conda env '{env_name}' not found. It must already exist "
            f"(conda/environment_*.yaml); the pipeline will not create or install it."
        )
    python = str(prefix / "bin" / "python")
    env = subprocess_env_for_prefix(prefix)
    if need_hkoca:
        env = ensure_hkoca_on_pythonpath(env, python)
    return python, env


def _conda_prefix_from_path(path: Path) -> Path | None:
    for parent in path.resolve().parents:
        if (parent / "conda-meta").is_dir():
            return parent
    return None


def _isolated_hkoca_share(pkg_dir: Path) -> Path:
    """PYTHONPATH dir that contains only ``hkoca``, never a full site-packages tree."""
    explicit = os.environ.get("HKOCA_SHARE", "").strip()
    if explicit:
        share = Path(explicit).expanduser().resolve()
    else:
        prefix = _conda_prefix_from_path(pkg_dir)
        if prefix is not None:
            share = prefix / "share" / "hkoca-pythonpath"
        else:
            share = Path.home() / ".cache" / "hkoca" / "pythonpath"
    share.mkdir(parents=True, exist_ok=True)
    link = share / "hkoca"
    target = pkg_dir.resolve()
    if link.exists() or link.is_symlink():
        try:
            if link.resolve() != target:
                link.unlink()
                link.symlink_to(target)
        except OSError:
            if not link.is_dir():
                raise
    else:
        link.symlink_to(target)
    return share


def hkoca_source_root() -> Path:
    """Directory to put on PYTHONPATH so other envs can ``import hkoca``.

    Must be a repo-style root (contains ``hkoca/``) or an isolated share with a
    ``hkoca`` symlink. Never a ``site-packages`` directory: that would leak
    numpy/torch from hkoca_harmonize into hkoca_projection.
    """
    override = os.environ.get("HKOCA_ROOT", "").strip()
    if override:
        root = Path(override).expanduser().resolve()
        if _looks_like_hkoca_root(root):
            return root

    import hkoca

    pkg = Path(hkoca.__file__).resolve().parent
    parent = pkg.parent
    if _looks_like_hkoca_root(parent):
        return parent

    for candidate in (Path("/opt/HKOCA-tools"), Path("/opt/hkoca")):
        if _looks_like_hkoca_root(candidate):
            return candidate.resolve()

    return _isolated_hkoca_share(pkg)


def wire_hkoca_pythonpath(env: dict[str, str]) -> dict[str, str]:
    """Set PYTHONPATH to the shared hkoca root only (no other env site-packages)."""
    root = str(hkoca_source_root())
    env["PYTHONPATH"] = root
    env["HKOCA_ROOT"] = root
    return env


def ensure_hkoca_on_pythonpath(env: dict[str, str], python: str) -> dict[str, str]:
    """Ensure ``hkoca`` is importable in a stage-specific conda env subprocess."""
    return wire_hkoca_pythonpath(env)


def projection_subprocess_env() -> tuple[str, dict[str, str]]:
    """Python + env for the hkoca_projection conda env."""
    env_name = os.environ.get("HKOCA_PROJECTION_ENV", ENV_PROJECTION).strip() or ENV_PROJECTION
    python, env = python_env(env_name, need_hkoca=True)
    return python, env
