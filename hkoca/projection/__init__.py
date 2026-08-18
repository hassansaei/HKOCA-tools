"""Project query datasets onto the HKOCA atlas (scPoli surgery)."""

from __future__ import annotations

import argparse
import logging
import sys

from hkoca.projection.config import load_projection_config, packaged_config_path
from hkoca.projection.runner import project_query
from hkoca.projection.stack import projection_stack_available

logger = logging.getLogger("hkoca.projection")


def status_message() -> str:
    return (
        "Atlas projection (scPoli surgery in hkoca_projection; hkoca shared via PYTHONPATH).\n"
        "  hkoca projection map --query sct_prepared.rds --atlas atlas.h5ad "
        "--model-dir scPoli_Reference_Model --output-dir out/projection\n"
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


def _hop_to_projection_env(argv: list[str]) -> int | None:
    """Run map in hkoca_projection once. Return None if already in that Python."""
    import os
    import subprocess
    from pathlib import Path

    from hkoca.conda_env import projection_subprocess_env

    if os.environ.get("HKOCA_PROJECTION_INNER") == "1":
        return None
    try:
        python, env = projection_subprocess_env()
    except FileNotFoundError as exc:
        logger.error("%s", exc)
        return 1
    if Path(python).resolve() == Path(sys.executable).resolve():
        return None
    env["HKOCA_PROJECTION_INNER"] = "1"
    logger.info("Using conda env: %s", env.get("CONDA_PREFIX", ""))
    return subprocess.run([python, "-m", "hkoca.projection", *argv], env=env).returncode


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hkoca projection",
        description=(
            "Project a query dataset onto the HKOCA kidney organoid atlas "
            "with scPoli surgery and prototype label transfer"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  # Run from hkoca_harmonize; torch/scArches run in hkoca_projection automatically\n"
            "  hkoca projection map \\\n"
            "      --query results/integration/prep/sct_prepared.rds \\\n"
            "      --atlas reference/Master_Atlas_scPoli_Integrated_Reannotated_fullgenes.h5ad \\\n"
            "      --model-dir reference/scPoli_Reference_Model \\\n"
            "      --output-dir results/projection\n"
            "  hkoca projection map --config my_projection.yaml --force\n"
        ),
    )
    sub = parser.add_subparsers(dest="cmd")

    p_map = sub.add_parser(
        "map",
        help="Map query cells onto the atlas with scPoli surgery",
    )
    p_map.add_argument(
        "--query",
        "--query-h5ad",
        dest="query_path",
        default=None,
        help="Query .h5ad or Seurat .rds (integration prep sct_prepared.rds uses RNA counts)",
    )
    p_map.add_argument(
        "--atlas",
        "--atlas-h5ad",
        dest="atlas_h5ad",
        default=None,
        help="HKOCA atlas .h5ad (Level_*_Integrated labels, X_scpoli / X_umap_scpoli)",
    )
    p_map.add_argument(
        "--model-dir",
        dest="model_dir",
        default=None,
        help="scPoli reference model directory (model_params.pt, attr.pkl, var_names.csv)",
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
        "--query-label-column",
        default=None,
        help="Query obs column for optional confusion-matrix comparison",
    )
    p_map.add_argument(
        "--query-batch-key",
        default=None,
        help="Query obs column used as scPoli condition (default: sample_id)",
    )
    p_map.add_argument(
        "--n-epochs",
        type=int,
        default=None,
        help="Surgery epochs (default: 50)",
    )
    p_map.add_argument(
        "--pretrain-epochs",
        type=int,
        default=None,
        help="Surgery pretraining epochs (default: 40)",
    )
    p_map.add_argument(
        "--eta",
        type=float,
        default=None,
        help="scPoli eta (default: 10)",
    )
    p_map.add_argument(
        "--joint-umap",
        action="store_true",
        help="Also compute exploratory joint latent UMAP (atlas subsample + query)",
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
    query = args.query_path or cfg.get("query_path")
    model_dir = args.model_dir or cfg.get("model_dir")
    output_dir = args.output_dir or cfg.get("output_dir")

    if args.n_epochs is not None:
        cfg["n_epochs"] = int(args.n_epochs)
    if args.pretrain_epochs is not None:
        cfg["pretrain_epochs"] = int(args.pretrain_epochs)
    if args.eta is not None:
        cfg["eta"] = float(args.eta)
    if args.joint_umap:
        cfg["joint_umap"] = True

    if not atlas:
        logger.error("Atlas h5ad required (--atlas or paths.atlas_h5ad in config).")
        return 1
    if not query:
        logger.error("Query .h5ad or .rds required (--query or paths.query in config).")
        return 1
    if not model_dir:
        logger.error("scPoli model dir required (--model-dir or paths.model_dir in config).")
        return 1
    if not output_dir:
        logger.error("Output directory required (-o/--output-dir or paths.output_dir in config).")
        return 1

    query_label_column = args.query_label_column or cfg.get("query_label_column")
    query_batch_key = args.query_batch_key or cfg.get("query_batch_key") or "sample_id"

    logger.info("Query: %s", query)
    logger.info("Atlas: %s", atlas)
    logger.info("Model: %s", model_dir)
    logger.info("Output: %s", output_dir)

    if args.dry_run:
        logger.info("Dry run: projection map would run in hkoca_projection (scPoli surgery).")
        return 0

    if not projection_stack_available():
        hopped = _hop_to_projection_env(argv)
        if hopped is not None:
            return hopped

    try:
        out_path = project_query(
            atlas_h5ad=atlas,
            query_path=query,
            model_dir=model_dir,
            output_dir=output_dir,
            query_label_column=query_label_column,
            query_batch_key=query_batch_key,
            cfg=cfg,
            force=bool(args.force_overwrite),
            joint_umap=True if args.joint_umap else None,
        )
    except Exception as exc:
        logger.error("Projection failed: %s", exc)
        if args.verbose:
            logger.exception("Projection error")
        return 1

    logger.info("Projected query saved: %s", out_path)
    return 0


__all__ = ["main", "status_message", "project_query"]
