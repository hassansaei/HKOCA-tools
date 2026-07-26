"""QC-filter stage: gene-space harmonization then doublet detection + QC.

Always runs both steps in order:

1. Harmonize gene space (CellRanger / H5AD / MTX -> Seurat RDS via --to-rds)
2. Doublet detection (scDblFinder) + QC threshold filtering (R)

Upstream CellBender is optional.

CLI::

    hkoca qc-filter --csv meta.csv --gtf genes.gtf --output results \\
      --qc-output results/qc_filter --stage all

Environments
------------
Use an env that provides scanpy (harmonize) and the QC R stack, or install
``conda/environment_harmonize.yaml`` plus QC packages from
``conda/environment_qc.yaml``.
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


def _ensure_to_rds(harmonize_argv: list[str]) -> list[str]:
    if "--to-rds" in harmonize_argv:
        return list(harmonize_argv)
    return [*harmonize_argv, "--to-rds"]


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


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hkoca qc-filter",
        description=(
            "Harmonize gene space (--to-rds), then run doublet detection + QC. "
            "Both steps always run in that order."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "order (fixed):\n"
            "  1. harmonize  gene space + Seurat RDS (--to-rds forced)\n"
            "  2. qc         scDblFinder + QC filters (R)\n"
            "\n"
            "harmonize options (also via harmonize.config):\n"
            "  --config, -w/--working-dir, --csv, --gtf, --output,\n"
            "  --summary, --summary-only, --transgenes\n"
            "\n"
            "example:\n"
            "  hkoca qc-filter --csv meta.csv --gtf genes.gtf --output results \\\n"
            "    --qc-output results/qc_filter --stage all\n"
        ),
    )
    parser.add_argument(
        "--qc-config",
        default=None,
        help="Path to qc_config.dcf (default: CWD then package qc_config.dcf)",
    )
    parser.add_argument(
        "--qc-output",
        default=None,
        help="QC output directory (default: <harmonize-output>/qc_filter)",
    )
    parser.add_argument(
        "--stage",
        default="all",
        choices=["all", "qc", "doublet"],
        help="QC R stage to run after harmonization (default: all)",
    )
    parser.add_argument(
        "--rds-pattern",
        default=r"_harmonized\.rds$",
        help="RDS basename regex under harmonize output (default: _harmonized\\.rds$)",
    )
    parser.add_argument("--force-overwrite", action="store_true", help="Force QC reprocessing")
    parser.add_argument("--dry-run", action="store_true", help="Print planned steps without running")
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
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
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = _build_parser()

    if argv and argv[0] in ("-h", "--help"):
        parser.parse_args(argv)
        return 0

    qc_opts, harmonize_argv = parser.parse_known_args(argv)
    _setup_logging(qc_opts.verbose)

    if qc_opts.print_script:
        print(r_script_path())
        return 0
    if qc_opts.print_config:
        print(config_template_path())
        return 0

    # Reject legacy split subcommands so usage stays a single command.
    if harmonize_argv and harmonize_argv[0] in ("run", "harmonize", "qc"):
        logger.error(
            "Split subcommands were removed. Use: hkoca qc-filter [options]\n"
            "Harmonization always runs first, then the QC R script."
        )
        return 2

    harmonize_argv = _ensure_to_rds(harmonize_argv)

    try:
        output_root = _resolve_harmonize_output(harmonize_argv)
    except SystemExit as exc:
        return int(exc.code or 1)
    except Exception as exc:
        logger.error("%s", exc)
        return 1

    qc_output = qc_opts.qc_output or os.path.join(output_root, "qc_filter")

    logger.info("QC-filter: 1) harmonize (--to-rds)  2) QC R script")
    logger.info("Harmonize output_root (RDS search root): %s", output_root)
    logger.info("QC output_dir: %s", qc_output)

    from hkoca.qc_filter.harmonize.cli import main as harmonize_main
    from hkoca.qc_filter.qc_runner import run_qc

    if qc_opts.dry_run:
        logger.info("[dry-run] would run harmonize: %s", " ".join(harmonize_argv))
    else:
        logger.info("STEP 1/2: harmonization")
        rc = int(harmonize_main(harmonize_argv) or 0)
        if rc != 0:
            logger.error("Harmonization failed; QC not started.")
            return rc

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


__all__ = [
    "config_template_path",
    "main",
    "r_script_path",
]
