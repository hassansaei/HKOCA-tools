#!/usr/bin/env python3
"""Rewrite sample_info data_path after CellBender (same rules as hkoca pipeline)."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path


PIPELINE_ONLY_COLS = {"sample_dir", "run_cellbender", "cellbender_mode", "skip", "output_dir"}
MANDATORY_META_COLS = (
    "data_path",
    "sample_id",
    "study",
    "source",
    "diff_protocol",
    "sc_protocol",
    "sequencing",
    "genome_build",
    "Age",
    "type",
)


def _as_bool(value, default: bool = True) -> bool:
    if value is None or str(value).strip() == "":
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def cellbender_filtered_output_path(sample_dir: str, sample_id: str, output_suffix: str) -> str:
    base = os.path.join(sample_dir, f"{sample_id}{output_suffix}")
    stem, ext = os.path.splitext(base)
    return f"{stem}_filtered{ext}"


def resolve_sample_dir(row: dict, working_dir: str) -> str:
    explicit = str(row.get("sample_dir", "") or "").strip()
    if explicit:
        path = explicit
    else:
        data_path = str(row["data_path"]).strip()
        if os.path.isfile(data_path) or data_path.lower().endswith((".h5", ".h5ad", ".hdf5")):
            path = os.path.dirname(data_path)
        else:
            path = data_path
    if working_dir and not os.path.isabs(path):
        path = os.path.normpath(os.path.join(working_dir, path))
    return path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample-info", required=True, help="Original sample_info CSV")
    p.add_argument("--output", required=True, help="Harmonize metadata CSV to write")
    p.add_argument(
        "--cellbender-map",
        default=None,
        help="TSV: sample_id<TAB>absolute_path_to_filtered_h5 (from Nextflow CELLBENDER)",
    )
    p.add_argument("--working-dir", default=os.getcwd())
    p.add_argument("--output-suffix", default="_filtered.h5")
    p.add_argument(
        "--skip-cellbender",
        action="store_true",
        help="Do not rewrite paths from CellBender; keep CSV data_path",
    )
    args = p.parse_args()

    cb_map: dict[str, str] = {}
    if args.cellbender_map and os.path.isfile(args.cellbender_map):
        with open(args.cellbender_map, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                cb_map[parts[0].strip()] = os.path.abspath(parts[1].strip())

    with open(args.sample_info, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            print("ERROR: empty CSV", file=sys.stderr)
            return 1
        fieldnames = list(reader.fieldnames)
        rows = list(reader)

    out_fields = [c for c in fieldnames if c not in PIPELINE_ONLY_COLS]
    for col in MANDATORY_META_COLS:
        if col not in out_fields:
            out_fields.append(col)
    if "file_prefix" in fieldnames and "file_prefix" not in out_fields:
        out_fields.append("file_prefix")

    out_rows: list[dict] = []
    for row in rows:
        sample_id = str(row.get("sample_id", "")).strip()
        record = {k: row.get(k, "") for k in out_fields}

        if sample_id in cb_map:
            record["data_path"] = cb_map[sample_id]
        elif args.skip_cellbender or not _as_bool(row.get("run_cellbender", "true"), True):
            path = str(row.get("data_path", "")).strip()
            if args.working_dir and path and not os.path.isabs(path):
                path = os.path.normpath(os.path.join(args.working_dir, path))
            record["data_path"] = path
        else:
            sample_dir = resolve_sample_dir(row, args.working_dir)
            record["data_path"] = cellbender_filtered_output_path(
                sample_dir, sample_id, args.output_suffix
            )

        out_rows.append(record)

    dest = Path(args.output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote harmonize metadata: {dest} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
