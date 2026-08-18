"""End-to-end A-Z pipeline orchestrator.

Chains stage modules in order:

    cellbender (optional) -> qc_filter -> annotation -> integration
    -> projection (optional, GPU)

Projection is off by default. Pass ``--run-projection`` (or set
``run_projection = true`` in pipeline.config) on a GPU node. The stage
uses the ``hkoca_projection`` conda env. Atlas h5ad and scPoli
``--model-dir`` are required when projection is enabled.

The sample_info CSV drives CellBender sample selection and harmonization
metadata. Use ``--skip-cellbender`` or per-row ``run_cellbender=False`` to
bypass ambient RNA removal.

By default the pipeline resumes: finished stage outputs are detected under
``--output`` and those steps are skipped. Pass ``--force`` to re-run.

Stage I/O (per study):
  qc_filter   -> qc_filter/qc_filtered_rds/{study}*_filtered.rds
                 qc_filter/h5ad_converted/{study}*_filtered.h5ad
  annotation  -> annotation/annotated_obj/{stem}_annotated.h5ad
  integration -> integration/prep/sct_prepared.rds (single-study layout)
                 integration/nonintegrated/*.png
                 integration/{harmony,rpca,cca}/*.png
                 integration/objects/integrated_{method}.rds
                 (multi-study: integration/{study}/...)
  projection  -> projection/projected_obj/sct_prepared_projected.h5ad
                 (multi-study: projection/{study}/...)
"""

from __future__ import annotations

import argparse
import logging
import sys

from hkoca.pipeline.config import (
    default_config_path,
    resolve_config,
    template_csv_path,
    validate_config,
)
from hkoca.pipeline.runner import run_pipeline

STAGES = (
    "cellbender",
    "qc_filter",
    "annotation",
    "integration",
    "projection",
)

logger = logging.getLogger("hkoca.pipeline")


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
        prog="hkoca pipeline",
        description=(
            "Run the full NephAtlas / HKOCA analysis pipeline from a sample_info CSV."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "sample_info CSV must include harmonize columns:\n"
            "  data_path, sample_id, study, source, diff_protocol, sc_protocol,\n"
            "  sequencing, genome_build, Age, type\n"
            "\n"
            "Optional pipeline columns:\n"
            "  sample_dir       CellRanger sample directory (default: dirname of data_path)\n"
            "  run_cellbender   True/False per sample (default: True)\n"
            "  file_prefix      MTX prefix when data_path is a directory\n"
            "  skip             True to exclude a row\n"
            "\n"
            "resume (default):\n"
            "  Finished stage outputs under --output are skipped.\n"
            "  Use --force to re-run everything.\n"
            "\n"
            "stages (default stops at integration; projection is optional GPU):\n"
            "  cellbender -> qc_filter -> annotation -> integration\n"
            "  add projection with --run-projection --atlas ATLAS.h5ad --model-dir MODEL_DIR\n"
            "\n"
            "example:\n"
            "  hkoca pipeline --csv sample_info.csv --gtf genes.gtf --output /data/out\n"
            "  hkoca pipeline --csv sample_info.csv --gtf genes.gtf --output /data/out \\\n"
            "    --skip-cellbender --to-stage qc_filter\n"
            "  hkoca pipeline --csv sample_info.csv --gtf genes.gtf --output /data/out --force\n"
            "  hkoca pipeline --csv sample_info.csv --gtf genes.gtf --output /data/out \\\n"
            "    --run-projection --atlas atlas.h5ad --model-dir scPoli_Reference_Model\n"
            "  hkoca pipeline --csv sample_info.csv --gtf genes.gtf --output /data/out \\\n"
            "    --from-stage projection --atlas atlas.h5ad --model-dir scPoli_Reference_Model\n"
        ),
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to pipeline.config (default: CWD then package pipeline.config)",
    )
    parser.add_argument(
        "-w",
        "--working-dir",
        default=None,
        help="Base directory for relative paths in the sample_info CSV",
    )
    parser.add_argument(
        "--csv",
        default=None,
        help="Sample info CSV (required unless set in pipeline.config)",
    )
    parser.add_argument(
        "--gtf",
        default=None,
        help="GRCh38 GTF for gene-space harmonization (required unless in config)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Root output directory for all pipeline stages (required unless in config)",
    )
    parser.add_argument(
        "--skip-cellbender",
        action="store_true",
        help="Skip CellBender for all samples (use raw or pre-filtered inputs in data_path)",
    )
    parser.add_argument(
        "--cellbender-mode",
        choices=["h5", "mtx"],
        default=None,
        help="CellBender input type (default: h5)",
    )
    parser.add_argument(
        "--from-stage",
        choices=STAGES,
        default=None,
        help="First stage to run (default: cellbender, or qc_filter with --skip-cellbender)",
    )
    parser.add_argument(
        "--to-stage",
        choices=STAGES,
        default=None,
        help="Last stage to run (default: integration; projection is opt-in)",
    )
    parser.add_argument(
        "--run-projection",
        action="store_true",
        help="Include the optional GPU projection stage after integration",
    )
    parser.add_argument(
        "--atlas",
        "--atlas-h5ad",
        dest="atlas_h5ad",
        default=None,
        help="HKOCA atlas .h5ad (required with --run-projection / --to-stage projection)",
    )
    parser.add_argument(
        "--model-dir",
        dest="model_dir",
        default=None,
        help="scPoli reference model directory (required with projection)",
    )
    parser.add_argument(
        "--qc-output",
        default=None,
        help="QC output directory (default: <output>/qc_filter)",
    )
    parser.add_argument(
        "--transgenes",
        default=None,
        help=(
            "Comma-separated transgenes for harmonize/integration "
            "(overrides [harmonize] transgenes in pipeline.config)"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned steps without executing them",
    )
    parser.add_argument(
        "--force",
        "--force-overwrite",
        dest="force",
        action="store_true",
        help="Re-run stages even when outputs already exist",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Disable resume checks (same as --force for skip logic)",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help=argparse.SUPPRESS,  # legacy alias; resume is now the default
    )
    parser.add_argument(
        "--print-template",
        action="store_true",
        help="Print path to the packaged sample_info_template.csv and exit",
    )
    parser.add_argument(
        "--list-stages",
        action="store_true",
        help="Print pipeline stages and exit",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = _build_parser()
    args = parser.parse_args(argv)
    _setup_logging(args.verbose)

    if args.print_template:
        print(template_csv_path())
        return 0

    if args.list_stages:
        for name in STAGES:
            print(name)
        return 0

    from_stage = args.from_stage
    if from_stage is None and args.skip_cellbender:
        from_stage = "qc_filter"

    try:
        cfg = resolve_config(
            default_config_path(args.config),
            working_dir=args.working_dir,
            sample_info_csv=args.csv,
            gtf_file=args.gtf,
            output_root=args.output,
            run_cellbender=False if args.skip_cellbender else None,
            cellbender_mode=args.cellbender_mode,
            from_stage=from_stage,
            to_stage=args.to_stage,
            run_projection=True if args.run_projection else None,
            projection_atlas_h5ad=args.atlas_h5ad,
            projection_model_dir=args.model_dir,
            qc_output=args.qc_output,
            transgenes=args.transgenes,
        )
        validate_config(cfg)
    except (ValueError, FileNotFoundError) as exc:
        logger.error("%s", exc)
        return 1

    force = bool(args.force or args.no_resume)
    return run_pipeline(
        cfg,
        dry_run=args.dry_run,
        resume=not force,
        force=force,
        verbose=args.verbose,
    )


__all__ = ["STAGES", "main"]
