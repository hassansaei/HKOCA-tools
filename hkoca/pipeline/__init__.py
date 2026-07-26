"""End-to-end A–Z pipeline orchestrator.

Chains public modules as they become available:

    cellbender → qc_filter (harmonize → doublets/QC) → annotation → integration
"""

from __future__ import annotations

STAGES = (
    "cellbender",
    "qc_filter",
    "annotation",
    "integration",
)


def status_message() -> str:
    lines = [
        "Full A–Z pipeline orchestrator (in progress).",
        "Planned stage order:",
        "  1. cellbender     ambient RNA removal (optional per sample)",
        "  2. qc_filter      harmonization → doublet detection → QC filtering",
        "  3. annotation     cell-type annotation",
        "  4. integration    batch integration / atlas projection",
        "",
        "Run a single module today with:",
        "  hkoca cellbender",
        "  hkoca qc-filter",
        "  hkoca annotation",
        "  hkoca integration",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """CLI entry for ``hkoca pipeline``."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="hkoca pipeline",
        description="Run the full NephAtlas / HKOCA analysis pipeline (A–Z)",
    )
    parser.add_argument(
        "--list-stages",
        action="store_true",
        help="Print planned pipeline stages and exit",
    )
    args = parser.parse_args(argv)

    if args.list_stages:
        for name in STAGES:
            print(name)
        return 0

    print(status_message())
    return 1


__all__ = ["STAGES", "main", "status_message"]
