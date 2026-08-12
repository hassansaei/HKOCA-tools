"""Batch integration and atlas projection (Seurat / Harmony / ...).

Stage 1 (``prep``): load QC-filtered Seurat RDS, normalize RNA, SCTransform,
PCA elbow, clustree, and silhouette-based resolution selection.
"""

from __future__ import annotations

import argparse
import logging
import sys

from hkoca.integration.integration_runner import (
    default_config_path,
    packaged_config,
    r_script_path,
    run_prep,
)

logger = logging.getLogger("hkoca.integration")


def status_message() -> str:
    return (
        "Integration module (Seurat R env).\n"
        "  hkoca integration prep --input-rds filtered.rds --output-dir out/integration\n"
        "  hkoca integration prep --input-rds filtered.rds --annotated-h5ad annotated.h5ad --output-dir out/integration"
    )


def _setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    if logger.handlers:
        logger.handlers.clear()
    logger.setLevel(level)
    logger.propagate = False
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%H:%M:%S",
        )
    )
    logger.addHandler(handler)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hkoca integration",
        description="Batch integration (Seurat R) and atlas projection",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  hkoca integration prep \\\n"
            "      --input-rds qc_filter/qc_filtered_rds/sample_filtered.rds \\\n"
            "      --output-dir results/integration\n"
            "  hkoca integration prep \\\n"
            "      --input-rds qc_filter/qc_filtered_rds/sample_filtered.rds \\\n"
            "      --annotated-h5ad annotation/annotated_obj/sample_annotated.h5ad \\\n"
            "      --output-dir results/integration\n"
        ),
    )
    sub = parser.add_subparsers(dest="cmd")

    p_prep = sub.add_parser(
        "prep",
        help="SCT + PCA elbow + clustree + silhouette resolution selection",
    )
    p_prep.add_argument(
        "--input-rds",
        required=True,
        help="QC-filtered Seurat RDS (e.g. *_harmonized_singlets_filtered.rds)",
    )
    p_prep.add_argument(
        "--annotated-h5ad",
        default=None,
        help="Annotated h5ad from hkoca annotation (optional; used in later stages)",
    )
    p_prep.add_argument(
        "--output-dir",
        "-o",
        required=True,
        help="Integration output directory",
    )
    p_prep.add_argument(
        "--config",
        default=None,
        help="Path to integration.config.dcf",
    )
    p_prep.add_argument(
        "--force",
        "--force-overwrite",
        dest="force_overwrite",
        action="store_true",
        help="Re-run prep even when outputs exist",
    )
    p_prep.add_argument("--dry-run", action="store_true")
    p_prep.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument(
        "--print-script",
        action="store_true",
        help="Print path to integration_prep.R and exit",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print path to packaged integration.config.dcf and exit",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = _build_parser()

    if argv and argv[0] in ("-h", "--help"):
        parser.parse_args(argv)
        return 0

    if not argv:
        parser.print_help()
        print()
        print(status_message())
        return 0

    if argv[0] == "--print-script":
        print(r_script_path("prep"))
        return 0
    if argv[0] == "--print-config":
        print(packaged_config())
        return 0

    args = parser.parse_args(argv)
    _setup_logging(bool(getattr(args, "verbose", False)))

    if args.cmd == "prep":
        return run_prep(
            input_rds=args.input_rds,
            output_dir=args.output_dir,
            annotated_h5ad=args.annotated_h5ad,
            config=args.config,
            force_overwrite=bool(args.force_overwrite),
            dry_run=bool(args.dry_run),
        )

    parser.print_help()
    return 1


__all__ = ["main", "status_message"]
