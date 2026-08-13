"""Batch integration and atlas projection (Seurat / Harmony / ...).

Stage 1 (``prep``): SplitObject by sample, SCTransform each sample, merge,
PCA elbow, clustree, and silhouette-based resolution selection.
Stage 2 (``run``): IntegrateLayers on the prep object (no second SCTransform),
then scIB metrics to rank Harmony / RPCA / CCA. Existing method RDS and
benchmark tables are skipped unless ``--force-overwrite``.
"""

from __future__ import annotations

import argparse
import logging
import sys

from hkoca.integration.integration_runner import (
    default_config_path,
    packaged_config,
    r_script_path,
    run_methods,
    run_prep,
)

logger = logging.getLogger("hkoca.integration")


def status_message() -> str:
    return (
        "Integration module (Seurat R env).\n"
        "  Stage 1 (prep):  hkoca integration prep --input-rds filtered.rds --output-dir out/integration\n"
        "  Stage 2 (run):   hkoca integration run --prepared-rds out/integration/prep/sct_prepared.rds \\\n"
        "                       --output-dir out/integration [--harmony] [--rpca] [--cca]\n"
        "                       (omit flags to run all three methods; scIB ranks them at the end)"
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
        help="Per-sample SCT + merge + PCA elbow/clustree/silhouette",
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
    p_prep.add_argument(
        "--remove-intermediate",
        action="store_true",
        help="Delete heavy intermediate files after successful prep (input RDS only if under output-dir)",
    )
    p_prep.add_argument("--dry-run", action="store_true")
    p_prep.add_argument("-v", "--verbose", action="store_true")

    p_run = sub.add_parser(
        "run",
        help="Run Harmony / RPCA / CCA integration + re-annotation on the prepared SCT object",
    )
    p_run.add_argument(
        "--prepared-rds",
        required=True,
        help="SCT-prepared RDS from 'hkoca integration prep' (prep/sct_prepared.rds)",
    )
    p_run.add_argument(
        "--output-dir",
        "-o",
        required=True,
        help="Integration output directory (same root used for prep)",
    )
    p_run.add_argument(
        "--methods",
        default=None,
        help="Comma-separated list of methods to run (default: harmony,rpca,cca)",
    )
    p_run.add_argument("--harmony", action="store_true", help="Run Harmony integration")
    p_run.add_argument("--rpca",    action="store_true", help="Run RPCA integration")
    p_run.add_argument("--cca",     action="store_true", help="Run CCA integration")
    p_run.add_argument(
        "--config",
        default=None,
        help="Path to integration.config.dcf",
    )
    p_run.add_argument(
        "--force",
        "--force-overwrite",
        dest="force_overwrite",
        action="store_true",
        help="Re-run even when per-method RDS or scIB benchmark already exists",
    )
    p_run.add_argument(
        "--remove-intermediate",
        action="store_true",
        help="Delete prep/sct_prepared.rds after successful integration run",
    )
    p_run.add_argument("--dry-run", action="store_true")
    p_run.add_argument("-v", "--verbose", action="store_true")

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
            remove_intermediate=bool(getattr(args, "remove_intermediate", False)),
            dry_run=bool(args.dry_run),
        )

    if args.cmd == "run":
        # Explicit flags take priority; fall back to --methods; default is all three.
        selected = [m for m, flag in [("harmony", args.harmony), ("rpca", args.rpca), ("cca", args.cca)] if flag]
        methods_str = ",".join(selected) if selected else getattr(args, "methods", None)
        return run_methods(
            prepared_rds=args.prepared_rds,
            output_dir=args.output_dir,
            methods=methods_str,
            config=args.config,
            force_overwrite=bool(args.force_overwrite),
            remove_intermediate=bool(getattr(args, "remove_intermediate", False)),
            dry_run=bool(args.dry_run),
        )

    parser.print_help()
    return 1


__all__ = ["main", "status_message", "run_prep", "run_methods"]
