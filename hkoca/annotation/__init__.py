"""Cell-type annotation with Snapseed (hierarchical marker YAML).

Works on QC-filtered or post-integration ``.h5ad`` objects. By default runs
Leiden at resolutions 0.4, 0.6, and 1.0, then Snapseed on each.

Marker YAML is updatable::

    hkoca annotation markers export ./my_markers.yaml
    # edit my_markers.yaml
    hkoca annotation run --input cells.h5ad --output-dir out --markers ./my_markers.yaml
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from hkoca.annotation.config import (
    DEFAULT_RESOLUTIONS,
    export_markers,
    load_annotation_config,
    load_marker_hierarchy,
    packaged_config_path,
    packaged_markers_path,
)

logger = logging.getLogger("hkoca.annotation")


def status_message() -> str:
    return (
        "Snapseed annotation helper.\n"
        "  hkoca annotation markers export markers.yaml\n"
        "  hkoca annotation run --input qc_filtered.h5ad --output-dir annotation_results"
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


def _parse_resolutions(raw: str | None) -> list[float] | None:
    if raw is None or not str(raw).strip():
        return None
    parts = [p.strip() for p in str(raw).replace(";", ",").split(",") if p.strip()]
    return [float(p) for p in parts]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hkoca annotation",
        description="Snapseed hierarchical cell-type annotation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  hkoca annotation markers export ./snapseed_markers.yaml\n"
            "  hkoca annotation run --input sample_filtered.h5ad --output-dir out \\\n"
            "      --markers ./snapseed_markers.yaml --resolutions 0.4,0.6,1.0\n"
            "  hkoca annotation run --input-dir qc_filtered_h5ad/ --output-dir out\n"
        ),
    )
    sub = parser.add_subparsers(dest="cmd")

    p_run = sub.add_parser("run", help="Cluster + Snapseed-annotate h5ad input(s)")
    src = p_run.add_mutually_exclusive_group(required=False)
    src.add_argument("--input", "-i", help="Single .h5ad file or glob")
    src.add_argument("--input-dir", help="Directory of .h5ad files")
    p_run.add_argument(
        "--output-dir",
        "-o",
        default=None,
        help="Output root (annotated_obj/, clustered_obj/, plots/)",
    )
    p_run.add_argument(
        "--markers",
        "-m",
        default=None,
        help="Marker hierarchy YAML (default: packaged snapseed_markers_v4.yaml)",
    )
    p_run.add_argument(
        "--config",
        default=None,
        help="Optional annotation.config.yaml",
    )
    p_run.add_argument(
        "-w",
        "--working-dir",
        default=None,
        help="Base directory for relative paths in config",
    )
    p_run.add_argument(
        "--resolutions",
        default=None,
        help="Comma-separated Leiden resolutions (default: 0.4,0.6,1.0)",
    )
    p_run.add_argument(
        "--save-plots",
        dest="save_plots",
        action="store_true",
        default=None,
        help="Write UMAP PNGs under plots/ (default: on)",
    )
    p_run.add_argument(
        "--no-save-plots",
        dest="save_plots",
        action="store_false",
        help="Disable UMAP PNG writing",
    )
    p_run.add_argument(
        "--force",
        "--force-overwrite",
        dest="force",
        action="store_true",
        help="Re-run even if annotated outputs exist",
    )
    p_run.add_argument(
        "--manual-overrides",
        default=None,
        help="Optional JSON file: {cluster_id: {Level_1: label, ...}}",
    )
    p_run.add_argument("-v", "--verbose", action="store_true")

    p_markers = sub.add_parser("markers", help="Inspect / export / validate marker YAML")
    markers_sub = p_markers.add_subparsers(dest="markers_cmd")

    p_show = markers_sub.add_parser("show", help="Print top-level categories")
    p_show.add_argument("--markers", "-m", default=None, help="Marker YAML path")

    p_export = markers_sub.add_parser(
        "export",
        help="Copy packaged (or source) markers to an editable path",
    )
    p_export.add_argument(
        "dest",
        nargs="?",
        default="snapseed_markers.yaml",
        help="Destination YAML path (default: ./snapseed_markers.yaml)",
    )
    p_export.add_argument(
        "--from",
        dest="source",
        default=None,
        help="Source YAML (default: packaged snapseed_markers_v4.yaml)",
    )

    p_validate = markers_sub.add_parser("validate", help="Validate marker YAML loads")
    p_validate.add_argument("path", help="Marker YAML to validate")

    p_path = markers_sub.add_parser("path", help="Print packaged markers path")
    p_cfg = markers_sub.add_parser("config-path", help="Print packaged config path")

    # Silence unused for type checkers / argparse wiring.
    _ = (p_path, p_cfg)
    return parser


def _cmd_markers(args: argparse.Namespace) -> int:
    cmd = args.markers_cmd
    if cmd is None or cmd == "show":
        path = Path(args.markers) if getattr(args, "markers", None) else packaged_markers_path()
        hierarchy = load_marker_hierarchy(path)
        print(f"Markers: {path}")
        print(f"Top-level categories ({len(hierarchy)}):")
        for key in hierarchy:
            print(f"  - {key}")
        return 0
    if cmd == "export":
        dest = export_markers(args.dest, source=args.source)
        print(f"Wrote editable markers YAML: {dest}")
        print("Edit this file, then pass it with: hkoca annotation run --markers PATH")
        return 0
    if cmd == "validate":
        hierarchy = load_marker_hierarchy(args.path)
        print(f"OK: {args.path} ({len(hierarchy)} top-level categories)")
        return 0
    if cmd == "path":
        print(packaged_markers_path())
        return 0
    if cmd == "config-path":
        print(packaged_config_path())
        return 0
    print(status_message())
    return 1


def _cmd_run(args: argparse.Namespace) -> int:
    from hkoca.annotation.runner import annotate_dataset, discover_h5ad_inputs, run_annotation_batch

    cfg = load_annotation_config(args.config, working_dir=args.working_dir)
    params = dict(cfg["parameters"])

    markers = Path(args.markers).expanduser().resolve() if args.markers else cfg["markers"]
    output_root = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else cfg["output_dir"]
    )
    annotated_dir = output_root / Path(cfg["annotated_dir"]).name
    clustered_dir = output_root / Path(cfg["clustered_dir"]).name
    plots_dir = output_root / Path(cfg["plots_dir"]).name

    if args.input:
        inputs = discover_h5ad_inputs(args.input)
    elif args.input_dir:
        inputs = discover_h5ad_inputs(args.input_dir)
    elif cfg["input"] is not None:
        inputs = discover_h5ad_inputs(cfg["input"])
    else:
        logger.error("Provide --input / --input-dir or set paths.input in config.")
        return 2

    resolutions = _parse_resolutions(args.resolutions)
    if resolutions is not None:
        params["resolutions"] = resolutions
    elif not params.get("resolutions"):
        params["resolutions"] = list(DEFAULT_RESOLUTIONS)

    if args.save_plots is not None:
        params["save_plots"] = bool(args.save_plots)
    if args.force:
        params["skip_existing"] = False

    manual = None
    if args.manual_overrides:
        with open(args.manual_overrides, encoding="utf-8") as fh:
            manual = json.load(fh)

    logger.info("Marker file : %s", markers)
    logger.info("Input files : %d", len(inputs))
    logger.info("Output dir  : %s", output_root)
    logger.info("Resolutions : %s", params["resolutions"])

    if manual:
        outputs = []
        hierarchy = load_marker_hierarchy(markers)
        for inp in inputs:
            outputs.append(
                annotate_dataset(
                    inp,
                    marker_dict=hierarchy,
                    markers_path=markers,
                    resolutions=params["resolutions"],
                    output_dir=annotated_dir,
                    clustered_dir=clustered_dir,
                    figures_dir=plots_dir,
                    manual_annotations=manual,
                    seed=int(params["seed"]),
                    hvg_n_top_genes=int(params["hvg_n_top_genes"]),
                    pca_n_comps=int(params["pca_n_comps"]),
                    neighbors_n_neighbors=int(params["neighbors_n_neighbors"]),
                    neighbors_n_pcs=int(params["neighbors_n_pcs"]),
                    scale_max_value=float(params["scale_max_value"]),
                    normalize_target_sum=float(params["normalize_target_sum"]),
                    save_plots=bool(params["save_plots"]),
                    dpi=int(params["dpi"]),
                    skip_existing=bool(params["skip_existing"]),
                    force=bool(args.force),
                )
            )
    else:
        outputs = run_annotation_batch(
            inputs,
            markers_path=markers,
            resolutions=params["resolutions"],
            annotated_dir=annotated_dir,
            clustered_dir=clustered_dir,
            figures_dir=plots_dir,
            parameters=params,
            force=bool(args.force),
        )

    logger.info("Finished annotation for %d file(s).", len(outputs))
    return 0


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = _build_parser()
    if not argv:
        parser.print_help()
        print()
        print(status_message())
        return 0

    args = parser.parse_args(argv)
    verbose = bool(getattr(args, "verbose", False))
    _setup_logging(verbose)

    if args.cmd == "markers":
        return _cmd_markers(args)
    if args.cmd == "run":
        return _cmd_run(args)

    parser.print_help()
    return 1


__all__ = ["main", "status_message"]
