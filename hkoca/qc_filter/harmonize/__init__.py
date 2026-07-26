"""Gene-space harmonization (internal step of qc_filter).

Not a public module. Invoked automatically by ``hkoca qc-filter`` (step 1).
"""

from __future__ import annotations

from typing import Any

__all__ = ["main", "parse_args", "run_pipeline", "run_summary"]


def __getattr__(name: str) -> Any:
    if name in ("main", "parse_args"):
        from hkoca.qc_filter.harmonize import cli

        return getattr(cli, name)
    if name in ("run_pipeline", "run_summary"):
        from hkoca.qc_filter.harmonize import pipeline

        return getattr(pipeline, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
