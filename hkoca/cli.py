"""Top-level CLI for the HKOCA modular toolkit."""

from __future__ import annotations

import argparse
import importlib
import sys

from hkoca._version import __version__

_COMMANDS = {
    "pipeline": ("hkoca.pipeline", "main"),
    "cellbender": ("hkoca.cellbender", "main"),
    "qc-filter": ("hkoca.qc_filter", "main"),
    "annotation": ("hkoca.annotation", "main"),
    "integration": ("hkoca.integration", "main"),
}


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = argparse.ArgumentParser(
        prog="hkoca",
        description=(
            "NephAtlas / HKOCA toolkit - modular single-cell organoid analysis"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "modules:\n"
            "  pipeline      End-to-end A-Z analysis\n"
            "  cellbender   Ambient RNA removal (h5 | mtx)\n"
            "  qc-filter    Harmonize then doublets + QC (run|harmonize|qc)\n"
            "  annotation   Cell-type annotation\n"
            "  integration  Batch integration / atlas projection\n"
        ),
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"hkoca {__version__}",
    )
    sub = parser.add_subparsers(dest="command")
    sub.add_parser("pipeline", help="Run the full A-Z analysis pipeline")
    sub.add_parser(
        "cellbender",
        help="Run CellBender ambient RNA removal (h5 or mtx input)",
    )
    sub.add_parser(
        "qc-filter",
        help="Harmonize then doublet/QC (run | harmonize | qc)",
    )
    sub.add_parser("annotation", help="Cell-type annotation")
    sub.add_parser("integration", help="Batch integration / atlas projection")

    if not argv or argv[0] in ("-h", "--help", "-V", "--version"):
        parser.parse_args(argv)
        if not argv:
            parser.print_help()
        return 0

    command, rest = argv[0], argv[1:]
    if command not in _COMMANDS:
        parser.error(f"unknown command: {command}")

    module_name, attr = _COMMANDS[command]
    entry = getattr(importlib.import_module(module_name), attr)
    return int(entry(rest) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
