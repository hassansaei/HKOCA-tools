"""Ambient RNA removal stage (CellBender only).

Two input CLIs
--------------
``hkoca cellbender h5``   CellRanger raw ``.h5`` per sample
``hkoca cellbender mtx``  CellRanger MTX directory (matrix/features/barcodes)

Environment: ``conda/environment_cellbender.yaml``
Cluster image: ``conda/cellbender_v0.3.2_updated.def``
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Sequence

from hkoca.cellbender.config import default_config_path, resolve_config
from hkoca.cellbender.runner import InputMode, run_jobs

logger = logging.getLogger("hkoca.cellbender")


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stdout,
        force=True,
    )


def _parse_samples_arg(raw: str | None) -> tuple[str, ...] | None:
    if raw is None or not raw.strip():
        return None
    return tuple(s.strip() for s in raw.split(",") if s.strip())


def _read_samples_file(path: str) -> tuple[str, ...]:
    samples: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            samples.append(line)
    return tuple(samples)


def _add_shared_args(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--config",
        default=None,
        help="Path to cellbender.config (default: CWD then package cellbender.config)",
    )
    p.add_argument(
        "-w",
        "--working-dir",
        default=None,
        help="Base directory for relative paths",
    )
    p.add_argument(
        "--samples-dir",
        default=None,
        help="Parent directory containing one subdirectory per sample",
    )
    p.add_argument(
        "--samples",
        default=None,
        help="Comma-separated sample IDs or paths (overrides config)",
    )
    p.add_argument(
        "--samples-file",
        default=None,
        help="Text file with one sample ID/path per line",
    )
    p.add_argument(
        "--epochs",
        type=int,
        default=None,
        help="Training epochs (default from config: 100)",
    )
    p.add_argument(
        "--total-droplets-included",
        type=int,
        default=None,
        help="Total droplets included (default from config: 35000)",
    )
    p.add_argument(
        "--learning-rate",
        type=float,
        default=None,
        help="Learning rate (default from config: 5e-5)",
    )
    p.add_argument(
        "--cpu-threads",
        type=int,
        default=None,
        help="CPU threads for CellBender (default from config: 24)",
    )
    cuda = p.add_mutually_exclusive_group()
    cuda.add_argument(
        "--cuda",
        action="store_true",
        help="Force GPU mode",
    )
    cuda.add_argument(
        "--no-cuda",
        action="store_true",
        help="Disable GPU",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running CellBender",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip samples whose output H5 already exists",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Debug logging",
    )


def _samples_from_args(args: argparse.Namespace) -> tuple[str, ...] | None:
    file_samples = _read_samples_file(args.samples_file) if args.samples_file else ()
    cli_samples = _parse_samples_arg(args.samples) or ()
    merged = tuple(dict.fromkeys([*file_samples, *cli_samples]))
    return merged or None


def _run_mode(mode: InputMode, argv: Sequence[str] | None) -> int:
    prog = f"hkoca cellbender {mode}"
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Run CellBender remove-background on CellRanger "
            + ("raw H5 files" if mode == "h5" else "MTX directories")
        ),
    )
    _add_shared_args(parser)
    args = parser.parse_args(list(argv) if argv is not None else None)
    _setup_logging(args.verbose)

    cuda: bool | None
    if args.no_cuda:
        cuda = False
    elif args.cuda:
        cuda = True
    else:
        cuda = None

    cfg = resolve_config(
        default_config_path(args.config),
        working_dir=args.working_dir,
        samples_dir=args.samples_dir,
        samples=_samples_from_args(args),
        epochs=args.epochs,
        total_droplets_included=args.total_droplets_included,
        learning_rate=args.learning_rate,
        cpu_threads=args.cpu_threads,
        cuda=cuda,
    )
    return run_jobs(
        cfg,
        mode,
        dry_run=args.dry_run,
        skip_existing=args.skip_existing,
    )


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    parser = argparse.ArgumentParser(
        prog="hkoca cellbender",
        description="Ambient RNA removal with CellBender",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "input modes:\n"
            "  h5    CellRanger raw_feature_bc_matrix.h5 per sample\n"
            "  mtx   CellRanger MTX directory (matrix/features/barcodes)\n"
        ),
    )
    sub = parser.add_subparsers(dest="mode")
    sub.add_parser("h5", help="Run CellBender on CellRanger raw H5 inputs")
    sub.add_parser("mtx", help="Run CellBender on CellRanger MTX directories")

    if not argv or argv[0] in ("-h", "--help"):
        parser.parse_args(argv if argv else ["-h"])
        return 0

    mode = argv[0]
    rest = argv[1:]
    if mode not in ("h5", "mtx"):
        parser.error(f"unknown mode: {mode} (choose h5 or mtx)")
    return _run_mode(mode, rest)  # type: ignore[arg-type]


__all__ = ["main"]
