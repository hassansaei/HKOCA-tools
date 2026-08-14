"""Verify HKOCA stage conda envs and cross-env tool access."""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from dataclasses import dataclass

from hkoca.conda_env import (
    ENV_CELLBENDER,
    ENV_HARMONIZE,
    ENV_INTEGRATION,
    ENV_PROJECTION,
    ENV_QC,
    ensure_hkoca_on_pythonpath,
    probe_projection_subprocess,
    resolve_env_prefix,
    resolve_python,
    subprocess_env_for_prefix,
    wire_harmonize_python,
)

logger = logging.getLogger("hkoca.env_verify")


@dataclass(frozen=True)
class EnvCheckResult:
    stage: str
    env: str
    ok: bool
    detail: str


def _probe_import(python: str, env: dict[str, str], statement: str) -> tuple[bool, str]:
    proc = subprocess.run(
        [python, "-c", statement],
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode == 0:
        return True, "ok"
    err = (proc.stderr or proc.stdout or "").strip()
    return False, err.splitlines()[-1] if err else "import failed"


def _check_env_exists(stage: str, env_name: str, executable: str) -> EnvCheckResult:
    prefix = resolve_env_prefix(env_name, executable)
    if prefix is None:
        return EnvCheckResult(
            stage,
            env_name,
            False,
            f"missing env (no {executable} in {env_name})",
        )
    return EnvCheckResult(stage, env_name, True, str(prefix))


def _check_python_modules(
    stage: str,
    env_name: str,
    modules: list[str],
    *,
    need_hkoca: bool = False,
    python_path: str | None = None,
    env: dict[str, str] | None = None,
) -> list[EnvCheckResult]:
    if python_path is None:
        python = resolve_python(env_name)
        if python is None:
            return [EnvCheckResult(stage, env_name, False, "python not found")]
        prefix = resolve_env_prefix(env_name, "python")
        assert prefix is not None
        run_env = subprocess_env_for_prefix(prefix)
    else:
        python = python_path
        run_env = dict(env or os.environ)
    if need_hkoca:
        run_env = ensure_hkoca_on_pythonpath(run_env, python)
    results: list[EnvCheckResult] = []
    for mod in modules:
        ok, detail = _probe_import(python, run_env, f"import {mod}")
        label = f"{env_name} ({mod})"
        results.append(EnvCheckResult(stage, label, ok, detail))
    return results


def verify_environments(*, include_projection: bool = True) -> list[EnvCheckResult]:
    """Run checks for each pipeline stage and cross-env wiring."""
    results: list[EnvCheckResult] = []

    results.append(_check_env_exists("harmonize / qc-filter / annotation", ENV_HARMONIZE, "python"))
    results.extend(
        _check_python_modules(
            "harmonize / annotation",
            ENV_HARMONIZE,
            ["scanpy", "anndata"],
        )
    )
    results.extend(
        _check_python_modules(
            "harmonize (hkoca package)",
            ENV_HARMONIZE,
            ["hkoca"],
            need_hkoca=True,
        )
    )

    results.append(_check_env_exists("qc-filter (R)", ENV_QC, "Rscript"))

    results.append(_check_env_exists("cellbender", ENV_CELLBENDER, "cellbender"))

    results.append(_check_env_exists("integration (R)", ENV_INTEGRATION, "Rscript"))

    ann_py = resolve_python(os.environ.get("HKOCA_ANNOTATION_ENV", ENV_HARMONIZE).strip() or ENV_HARMONIZE)
    harm_prefix = resolve_env_prefix(ENV_HARMONIZE, "python")
    if ann_py is None or harm_prefix is None:
        results.append(
            EnvCheckResult(
                "integration (reticulate -> harmonize)",
                ENV_HARMONIZE,
                False,
                "harmonize python not found for annotated h5ad transfer",
            )
        )
    else:
        env = wire_harmonize_python(subprocess_env_for_prefix(harm_prefix))
        ok, detail = _probe_import(ann_py, env, "import scanpy")
        results.append(
            EnvCheckResult(
                "integration (reticulate -> harmonize)",
                f"{ENV_HARMONIZE} scanpy",
                ok,
                detail,
            )
        )

    results.extend(
        _check_python_modules(
            "integration (scIB benchmark)",
            ENV_HARMONIZE,
            ["scib"],
        )
    )
    results.extend(
        _check_python_modules(
            "integration (scIB hkoca module)",
            ENV_HARMONIZE,
            ["hkoca.integration.benchmark"],
            need_hkoca=True,
        )
    )

    if include_projection:
        ok, detail = probe_projection_subprocess()
        results.append(
            EnvCheckResult(
                "projection (hkoca_projection + shared hkoca)",
                ENV_PROJECTION,
                ok,
                detail,
            )
        )
        int_r = resolve_env_prefix(ENV_INTEGRATION, "Rscript")
        if int_r is None:
            results.append(
                EnvCheckResult(
                    "projection (query RDS convert)",
                    ENV_INTEGRATION,
                    False,
                    "Rscript not found for Seurat RDS -> h5ad",
                )
            )
        else:
            results.append(
                EnvCheckResult(
                    "projection (query RDS convert)",
                    ENV_INTEGRATION,
                    True,
                    str(int_r),
                )
            )

    return results


def format_report(results: list[EnvCheckResult]) -> str:
    width = max(len(r.stage) for r in results) if results else 20
    lines = [f"{'STAGE':<{width}}  STATUS  DETAIL", "-" * (width + 40)]
    for row in results:
        status = "OK" if row.ok else "FAIL"
        lines.append(f"{row.stage:<{width}}  {status:<5}  {row.detail}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    include_projection = "--no-projection" not in argv
    verbose = "-v" in argv or "--verbose" in argv

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(message)s",
    )

    results = verify_environments(include_projection=include_projection)
    report = format_report(results)
    print(report)

    failed = [r for r in results if not r.ok]
    if failed:
        print(
            f"\n{len(failed)} check(s) failed. Create missing envs from conda/environment_*.yaml. "
            "Install hkoca once in hkoca_harmonize (pip install -e .) and set HKOCA_ROOT "
            "to the repo root if stage envs cannot import hkoca.",
            file=sys.stderr,
        )
        return 1
    print("\nAll environment checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
