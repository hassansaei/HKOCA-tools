"""CLI for gene-space harmonization (no heavy imports until run)."""

from __future__ import annotations

import argparse
import logging
import os
import sys

from hkoca.qc_filter.harmonize.config import (
    config_path,
    load_config,
    resolve_paths,
    resolve_transgenes,
)

logger = logging.getLogger("harmonize_atlas")


def parse_args(argv: list[str] | None = None):
    p = argparse.ArgumentParser(
        prog="hkoca qc-filter",
        description="scRNA-seq atlas gene-space harmonization (qc-filter step 1)",
    )
    p.add_argument(
        "--config",
        default=None,
        help="Path to harmonize.config (default: CWD then package harmonize.config)",
    )
    p.add_argument(
        "--working-dir",
        "-w",
        dest="working_dir",
        default=None,
        help="Working directory for relative paths (overrides [paths] working_dir)",
    )
    p.add_argument(
        "--csv",
        default=None,
        help="Metadata CSV (overrides [paths] metadata_csv)",
    )
    p.add_argument(
        "--gtf",
        default=None,
        help="GRCh38 GTF path (overrides [paths] gtf_file)",
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output root directory (overrides [paths] output_root)",
    )
    p.add_argument("--summary", action="store_true", help="Generate summary plots after pipeline")
    p.add_argument("--summary-only", action="store_true", help="Skip pipeline, only generate plots")
    p.add_argument("--to-rds", action="store_true", help="Convert harmonized .h5ad to Seurat .rds")
    p.add_argument(
        "--skip-existing",
        action="store_true",
        default=True,
        help="Skip studies whose harmonized RDS (or h5ad) already exists (default: on)",
    )
    p.add_argument(
        "--force",
        "--force-overwrite",
        dest="force",
        action="store_true",
        help="Re-run even when harmonized outputs already exist",
    )
    p.add_argument(
        "--transgenes",
        default=None,
        help="Comma-separated transgenes (overrides [transgenes] names)",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    cfg_path = config_path(args.config)
    cfg = load_config(cfg_path)
    paths = resolve_paths(args, cfg)
    working_dir = paths["working_dir"]
    gtf_file = paths["gtf_file"]
    metadata_csv = paths["metadata_csv"]
    output_root = paths["output_root"]

    missing = []
    if not gtf_file:
        missing.append("GTF path (--gtf or [paths] gtf_file)")
    if not args.summary_only:
        if not metadata_csv:
            missing.append("Metadata CSV (--csv or [paths] metadata_csv)")
        if not output_root:
            missing.append("Output directory (--output or [paths] output_root)")
    elif not output_root:
        missing.append(
            "Output directory for --summary-only (--output or [paths] output_root)"
        )
    if missing:
        # Avoid importing scanpy just to report missing inputs / --help already exited
        print("ERROR: Missing required inputs: " + ", ".join(missing), file=sys.stderr)
        print(f"Config file: {cfg_path}", file=sys.stderr)
        return 1

    if metadata_csv and not os.path.isabs(metadata_csv):
        metadata_csv = os.path.normpath(os.path.join(working_dir, metadata_csv))
    if gtf_file and not os.path.isabs(gtf_file):
        gtf_file = os.path.normpath(os.path.join(working_dir, gtf_file))
    if output_root and not os.path.isabs(output_root):
        output_root = os.path.normpath(os.path.join(working_dir, output_root))

    transgene_names = resolve_transgenes(args, cfg)

    # Heavy deps loaded only when a real run is requested
    from hkoca.qc_filter.harmonize import pipeline as P

    P.setup_logging(working_dir)
    logger.info("Python         : %s", sys.version.split()[0])
    logger.info(
        "Config file    : %s%s",
        cfg_path,
        "" if os.path.isfile(cfg_path) else " (not found - using defaults)",
    )
    logger.info("=== Atlas Harmonization Initialized ===")
    logger.info("Working dir  : %s", working_dir)
    logger.info("GTF file     : %s", gtf_file)
    logger.info("Metadata CSV : %s", metadata_csv)
    logger.info("Output root  : %s", output_root)
    if transgene_names:
        logger.info("Transgenes added to reference: %s", sorted(transgene_names))

    if not args.summary_only:
        failed = P.run_pipeline(
            metadata_csv=metadata_csv,
            gtf_file=gtf_file,
            output_root=output_root,
            working_dir=working_dir,
            to_rds=args.to_rds,
            transgene_names=transgene_names,
            skip_existing=bool(args.skip_existing and not args.force),
        )
        if failed:
            return 1

    if args.summary or args.summary_only:
        logger.info("Generating summary plots...")
        P.run_summary(scan_path=output_root or ".", cfg=cfg)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
