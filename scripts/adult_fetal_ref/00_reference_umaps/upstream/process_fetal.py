#!/usr/bin/env python
"""Upstream fetal: Tran QC → assemble → Harmony integrate → sync figures.

One entry point for the fetal Python pipeline (R DevKidCC stays separate).

Stages
------
  sync       copy already-built Harmony figures/tables into 00_reference_umaps/
  tran_qc    GSE124472 MTX → QC (fetal_new_3_datasets/process_tran.py)
  assemble   DevKidCC MTX+CSV → tran.h5ad
  integrate  Stewart+Hochane+Tran → Harmony weeks 9–17
  all        tran_qc → assemble → integrate → sync
             (run DevKidCC R between tran_qc and assemble)

Usage
-----
  python scripts/00_reference_umaps/upstream/process_fetal.py          # default: sync
  python scripts/00_reference_umaps/upstream/process_fetal.py --stage integrate
  # DevKidCC (separate):
  Rscript scripts/00_reference_umaps/upstream/run_devkidcc_tran.R
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib_config import DIRS, WORKSPACE, ensure_dirs  # noqa: E402

FETAL_ROOT = WORKSPACE / "fetal_new_3_datasets"
INTEG = FETAL_ROOT / "integrated_weeks_9_17"
H5AD = INTEG / "fetal_weeks_9_17_harmony_integrated.h5ad"
ATLAS_PY = Path("/data-master/pure-workspace/hmami/conda_envs/atlas/bin/python")

# Original modules live under fetal_new_3_datasets/ — this script is the single driver
STEP_SCRIPTS = {
    "tran_qc": FETAL_ROOT / "process_tran.py",
    "assemble": FETAL_ROOT / "assemble_tran_devkidcc_h5ad.py",
    "integrate": FETAL_ROOT / "harmonize_and_integrate.py",
}
R_DEVKIDCC = FETAL_ROOT / "run_devkidcc_tran.R"

FIGURES = [
    "umap_dataset.png",
    "umap_week.png",
    "umap_harmonized_celltype.png",
    "umap_leiden_res0.6.png",
]
TABLES = [
    "celltype_harmonization_dictionary.csv",
    "harmonized_celltype_counts.csv",
    "harmonized_celltype_by_dataset.csv",
    "harmonized_celltype_by_week.csv",
    "cells_by_week_dataset.csv",
]


def _py() -> str:
    return str(ATLAS_PY if ATLAS_PY.exists() else Path(sys.executable))


def run_step(name: str) -> None:
    script = STEP_SCRIPTS[name]
    if not script.exists():
        raise FileNotFoundError(script)
    print(f"=== Fetal stage: {name} ({script.name}) ===", flush=True)
    subprocess.run([_py(), "-u", str(script)], check=True, cwd=str(FETAL_ROOT))


def sync_outputs() -> None:
    ensure_dirs()
    fig_dir = DIRS["00_fig"] / "fetal_integration"
    tab_dir = DIRS["00_tab"] / "fetal_integration"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    print("=== Fetal stage: sync (figures → 00_reference_umaps) ===", flush=True)
    if not H5AD.exists():
        raise FileNotFoundError(
            f"Missing {H5AD}\n"
            "Rebuild: python scripts/00_reference_umaps/upstream/process_fetal.py --stage integrate"
        )

    for name in FIGURES:
        src = INTEG / "figures" / name
        if src.exists():
            shutil.copy2(src, fig_dir / name)
            print(f"  figure → {fig_dir / name}", flush=True)
    for name in TABLES:
        src = INTEG / name
        if src.exists():
            shutil.copy2(src, tab_dir / name)
            print(f"  table  → {tab_dir / name}", flush=True)
    summary = INTEG / "INTEGRATION_SUMMARY.md"
    if summary.exists():
        shutil.copy2(summary, tab_dir / "INTEGRATION_SUMMARY.md")

    week_tab = tab_dir / "cells_by_week_dataset.csv"
    if week_tab.exists():
        print(pd.read_csv(week_tab).to_string(index=False), flush=True)

    prov = {
        "collection": {
            "Stewart": "prebuilt stewart.h5ad",
            "Hochane": "prebuilt hochane.h5ad",
            "Tran": "GSE124472 MTX → tran_qc → DevKidCC (R) → assemble → tran.h5ad",
        },
        "weeks_kept": "9–17",
        "integration": "PCA → Harmony (batch=dataset) → UMAP + Leiden 0.6",
        "object": str(H5AD),
        "driver": str(Path(__file__).resolve()),
        "devkidcc_R": str(R_DEVKIDCC),
        "figures": str(fig_dir),
    }
    (DIRS["00_tab"] / "fetal_reference_provenance.json").write_text(json.dumps(prov, indent=2))
    print(f"Wrote {DIRS['00_tab'] / 'fetal_reference_provenance.json'}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--stage",
        choices=["sync", "tran_qc", "assemble", "integrate", "all"],
        default="sync",
        help="Pipeline stage (default: sync figures into 00)",
    )
    args = ap.parse_args()

    if args.stage == "sync":
        sync_outputs()
    elif args.stage == "all":
        run_step("tran_qc")
        print(
            f"\n>>> Run DevKidCC next:\n"
            f"    Rscript {R_DEVKIDCC}\n"
            f">>> Then re-run with --stage assemble (or continue after R).\n",
            flush=True,
        )
        run_step("assemble")
        run_step("integrate")
        sync_outputs()
    else:
        run_step(args.stage)
        if args.stage == "integrate":
            sync_outputs()

    print("Fetal upstream done.", flush=True)


if __name__ == "__main__":
    main()
