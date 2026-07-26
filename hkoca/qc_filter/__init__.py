"""QC-filter stage: gene-space harmonization then doublet detection + QC.

Order
-----
1. Harmonize gene space (CellRanger / H5AD / MTX -> Seurat RDS via --to-rds)
2. Doublet detection (scDblFinder) + QC threshold filtering (R)

Upstream CellBender is optional.

CLI
---
  hkoca qc-filter run ...          harmonize --to-rds, then QC R script
  hkoca qc-filter harmonize ...    harmonization only
  hkoca qc-filter qc ...           QC R script only (expects existing RDS)

Environments
------------
- Harmonization: ``conda/environment_harmonize.yaml`` (includes R for --to-rds)
- Full QC R stack: ``conda/environment_qc.yaml``
  For ``run``, use an env that provides both scanpy and the QC R packages,
  or run ``harmonize`` then ``qc`` in their respective envs.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from importlib import resources
from pathlib import Path

logger = logging.getLogger("hkoca.qc_filter")


def r_script_path() -> Path:
    """Return the filesystem path to ``QC_scdbl_Combined.R``."""
    return Path(resources.files("hkoca.qc_filter.r").joinpath("QC_scdbl_Combined.R")).resolve()


def config_template_path() -> Path:
    """Return the filesystem path to packaged ``qc_config.dcf``."""
    return Path(resources.files("hkoca.qc_filter.r").joinpath("qc_config.dcf")).resolve()


def _setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
        force=True,
    )


def status_message() -> str:
    script = r_script_path()
    cfg = config_template_path()
    return (
        "QC-filter module - always harmonize first, then doublet detection + QC.\n"
        "CellBender upstream is optional.\n"
        "\n"
        "Subcommands:\n"
        "  hkoca qc-filter run ...         harmonize (--to-rds) then QC R script\n"
        "  hkoca qc-filter harmonize ...   gene-space harmonization only\n"
        "  hkoca qc-filter qc ...          QC R script only (existing RDS)\n"
        "\n"
        f"R script : {script}\n"
        f"QC config: {cfg}\n"
        "\n"
        "Example (full stage):\n"
        "  hkoca qc-filter run \\\n"
        "    --csv meta.csv --gtf genes.gtf --output results \\\n"
        "    --qc-output results/qc_filter --stage all"
    )


def _split_run_argv(argv: list[str]) -> tuple[list[str], argparse.Namespace]:
    """Separate qc-filter run options from harmonize passthrough args."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--qc-config", default=None)
    parser.add_argument("--qc-output", default=None)
    parser.add_argument("--stage", default="all", choices=["all", "qc", "doublet"])
    parser.add_argument("--skip-harmonize", action="store_true")
    parser.add_argument("--skip-qc", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force-overwrite", action="store_true")
    parser.add_argument("--rds-pattern", default=r"_harmonized\.rds$")
    parser.add_argument("-v", "--verbose", action="store_true")

    qc_ns, harmonize_argv = parser.parse_known_args(argv)
    return harmonize_argv, qc_ns


def _resolve_harmonize_output(harmonize_argv: list[str]) -> str:
    from hkoca.qc_filter.harmonize.cli import parse_args
    from hkoca.qc_filter.harmonize.config import config_path, load_config, resolve_paths

    args = parse_args(harmonize_argv)
    cfg = load_config(config_path(args.config))
    paths = resolve_paths(args, cfg)
    output_root = paths["output_root"]
    if not output_root:
        raise ValueError(
            "Harmonize output_root is required. Pass --output or set "
            "[paths] output_root in harmonize.config."
        )
    working_dir = paths["working_dir"]
    if not os.path.isabs(output_root):
        output_root = os.path.normpath(os.path.join(working_dir, output_root))
    return output_root


def _ensure_to_rds(harmonize_argv: list[str]) -> list[str]:
    if "--to-rds" in harmonize_argv:
        return list(harmonize_argv)
    return [*harmonize_argv, "--to-rds"]


def _cmd_run(argv: list[str]) -> int:
    harmonize_argv, qc_opts = _split_run_argv(argv)
    _setup_logging(qc_opts.verbose)
    harmonize_argv = _ensure_to_rds(harmonize_argv)

    try:
        output_root = _resolve_harmonize_output(harmonize_argv)
    except SystemExit as exc:
        return int(exc.code or 1)
    except Exception as exc:
        logger.error("%s", exc)
        return 1

    qc_output = qc_opts.qc_output or os.path.join(output_root, "qc_filter")

    logger.info("QC-filter run order: 1) harmonize (--to-rds)  2) QC R script")
    logger.info("Harmonize output_root (RDS search root): %s", output_root)
    logger.info("QC output_dir: %s", qc_output)

    if not qc_opts.skip_harmonize:
        from hkoca.qc_filter.harmonize.cli import main as harmonize_main

        if qc_opts.dry_run:
            logger.info("[dry-run] would run: hkoca qc-filter harmonize %s", " ".join(harmonize_argv))
        else:
            logger.info("STEP 1/2: harmonization")
            rc = int(harmonize_main(harmonize_argv) or 0)
            if rc != 0:
                logger.error("Harmonization failed; skipping QC.")
                return rc
    else:
        logger.info("Skipping harmonization (--skip-harmonize)")

    if qc_opts.skip_qc:
        logger.info("Skipping QC (--skip-qc)")
        return 0

    from hkoca.qc_filter.qc_runner import run_qc

    logger.info("STEP 2/2: doublet detection + QC filtering")
    return run_qc(
        rds_dir=output_root,
        output_dir=qc_output,
        config=qc_opts.qc_config,
        stage=qc_opts.stage,
        recursive=True,
        rds_pattern=qc_opts.rds_pattern,
        force_overwrite=qc_opts.force_overwrite,
        dry_run=qc_opts.dry_run,
    )


def _cmd_qc(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="hkoca qc-filter qc",
        description="Run doublet detection + QC filtering on existing Seurat RDS",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to qc_config.dcf (default: CWD then package qc_config.dcf)",
    )
    parser.add_argument(
        "--rds-dir",
        default=None,
        help="Directory containing input .rds files (harmonize output_root)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="QC output base directory",
    )
    parser.add_argument(
        "--stage",
        default="all",
        choices=["all", "qc", "doublet"],
        help="QC stage to run (default: all)",
    )
    parser.add_argument(
        "--rds-pattern",
        default=r"_harmonized\.rds$",
        help="Regex for RDS basenames (default: _harmonized\\.rds$)",
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        help="Do not search rds-dir recursively",
    )
    parser.add_argument("--force-overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--print-script",
        action="store_true",
        help="Print path to QC_scdbl_Combined.R and exit",
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print path to packaged qc_config.dcf and exit",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    _setup_logging(args.verbose)

    if args.print_script:
        print(r_script_path())
        return 0
    if args.print_config:
        print(config_template_path())
        return 0

    if not args.rds_dir or not args.output_dir:
        parser.error("--rds-dir and --output-dir are required unless using --print-script/--print-config")

    from hkoca.qc_filter.qc_runner import run_qc

    return run_qc(
        rds_dir=args.rds_dir,
        output_dir=args.output_dir,
        config=args.config,
        stage=args.stage,
        recursive=not args.no_recursive,
        rds_pattern=args.rds_pattern,
        force_overwrite=args.force_overwrite,
        dry_run=args.dry_run,
    )


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = argparse.ArgumentParser(
        prog="hkoca qc-filter",
        description=(
            "Harmonization then doublet detection + QC filtering "
            "(optional CellBender upstream)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "order:\n"
            "  1. harmonize  (gene space + optional --to-rds)\n"
            "  2. qc         (scDblFinder + QC filters via R)\n"
        ),
    )
    sub = parser.add_subparsers(dest="subcommand")
    sub.add_parser(
        "run",
        help="Run harmonize (--to-rds) then QC R script",
    )
    sub.add_parser(
        "harmonize",
        help="Gene-space harmonization only (use --to-rds for Seurat RDS)",
    )
    sub.add_parser(
        "qc",
        help="Run QC R script only on existing RDS",
    )

    if not argv or argv[0] in ("-h", "--help"):
        if argv and argv[0] in ("-h", "--help"):
            parser.parse_args(argv)
            return 0
        print(status_message())
        return 0

    subcommand = argv[0]
    rest = argv[1:]

    if subcommand == "run":
        return _cmd_run(rest)

    if subcommand == "harmonize":
        from hkoca.qc_filter.harmonize.cli import main as harmonize_main

        return int(harmonize_main(rest) or 0)

    if subcommand == "qc":
        return _cmd_qc(rest)

    parser.error(f"unknown subcommand: {subcommand}")
    return 2


__all__ = [
    "config_template_path",
    "main",
    "r_script_path",
    "status_message",
]
