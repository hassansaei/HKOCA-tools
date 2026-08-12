"""Project query datasets onto a reference atlas (scanpy ingest)."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from hkoca.env_check import log_env_problems, verify_stage_env
from hkoca.projection.config import load_projection_config, packaged_config_path
from hkoca.projection.runner import project_query

logger = logging.getLogger("hkoca.projection")


def status_message() -> str:
    return (
        "Atlas projection (scanpy ingest on hkoca_harmonize env).\n"
        "  hkoca projection map --atlas atlas.h5ad --query query.h5ad --output-dir out/projection\n"
        "  hkoca projection map --config projection.config.yaml"
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


def _parse_label_columns(raw: str | None) -> list[str] | None:
    if raw is None or not str(raw).strip():
        return None
    return [c.strip() for c in str(raw).replace(";", ",").split(",") if c.strip()]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hkoca projection",
        description="Project query h5ad onto a reference atlas and generate comparison plots",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  hkoca projection map \\\n"
            "      --atlas reference/nephatlas.h5ad \\\n"
            "      --query results/annotation/annotated_obj/sample_annotated.h5ad \\\n"
            "      --output-dir results/projection\n"
            "  hkoca projection map --config my_projection.yaml \\\n"
            "      --query external_query.h5ad --output-dir results/projection \\\n"
            "      --query-label-column celltype --force\n"
        ),
    )
    sub = parser.add_subparsers(dest="cmd")

    p_map = sub.add_parser(
        "map",
        help="Map query cells onto atlas labels and UMAP (scanpy ingest)",
    )
    p_map.add_argument(
        "--atlas",
        "--atlas-h5ad",
        dest="atlas_h5ad",
        default=None,
        help="Reference atlas .h5ad (cell-type labels in obs)",
    )
    p_map.add_argument(
        "--query",
        "--query-h5ad",
        dest="query_h5ad",
        default=None,
        help="Query .h5ad (HKOCA pipeline or external)",
    )
    p_map.add_argument(
        "--output-dir",
        "-o",
        default=None,
        help="Output directory (default: projection_results or config paths.output_dir)",
    )
    p_map.add_argument(
        "--config",
        default=None,
        help="Path to projection.config.yaml",
    )
    p_map.add_argument(
        "--label-columns",
        default=None,
        help="Comma-separated atlas obs columns to transfer (default: auto-detect)",
    )
    p_map.add_argument(
        "--query-label-column",
        default=None,
        help="Query obs column for ground-truth comparison (confusion matrix)",
    )
    p_map.add_argument(
        "--query-batch-key",
        default=None,
        help="Query obs column for split UMAPs (default: sample_id)",
    )
    p_map.add_argument(
        "--force",
        "--force-overwrite",
        dest="force_overwrite",
        action="store_true",
        help="Re-run even when projected h5ad already exists",
    )
    p_map.add_argument("--dry-run", action="store_true")
    p_map.add_argument("-v", "--verbose", action="store_true")

    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print path to packaged projection.config.yaml and exit",
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

    if argv[0] == "--print-config":
        print(packaged_config_path())
        return 0

    args = parser.parse_args(argv)
    _setup_logging(bool(getattr(args, "verbose", False)))

    if args.cmd != "map":
        parser.print_help()
        return 1

    cfg = load_projection_config(args.config)
    atlas = args.atlas_h5ad or cfg.get("atlas_h5ad")
    query = args.query_h5ad or cfg.get("query_h5ad")
    output_dir = args.output_dir or cfg.get("output_dir")

    if not atlas:
        logger.error("Atlas h5ad required (--atlas or paths.atlas_h5ad in config).")
        return 1
    if not query:
        logger.error("Query h5ad required (--query or paths.query_h5ad in config).")
        return 1
    if not output_dir:
        logger.error("Output directory required (-o/--output-dir or paths.output_dir in config).")
        return 1

    label_columns = _parse_label_columns(args.label_columns) or cfg.get("reference_label_columns")
    query_label_column = args.query_label_column or cfg.get("query_label_column")
    query_batch_key = args.query_batch_key or cfg.get("query_batch_key") or "sample_id"

    logger.info("Atlas: %s", atlas)
    logger.info("Query: %s", query)
    logger.info("Output: %s", output_dir)

    if args.dry_run:
        logger.info("Dry run: projection map would run with hkoca_harmonize env.")
        return 0

    problems = verify_stage_env("harmonize")
    if problems:
        log_env_problems("harmonize", problems)
        return 1

    try:
        out_path = project_query(
            atlas_h5ad=atlas,
            query_h5ad=query,
            output_dir=output_dir,
            label_columns=label_columns,
            query_label_column=query_label_column,
            query_batch_key=query_batch_key,
            cfg=cfg,
            force=bool(args.force_overwrite),
        )
    except Exception as exc:
        logger.error("Projection failed: %s", exc)
        if args.verbose:
            logger.exception("Projection error")
        return 1

    logger.info("Projected query saved: %s", out_path)
    return 0


__all__ = ["main", "status_message", "project_query"]
