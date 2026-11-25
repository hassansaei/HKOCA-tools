#!/usr/bin/env python
"""Upstream fetal reference — produced entirely under 00_reference_umaps/.

Pipeline (local only)
---------------------
1. Place source AnnData under:
     organoid_vs_primary_atlas/00_reference_umaps/inputs/fetal/
       stewart.h5ad  hochane.h5ad  tran.h5ad
2. Run this script to filter weeks 9–17, Harmony-integrate, and write
   objects / figures / tables under 00_reference_umaps/.

Optional Tran rebuild (if tran.h5ad is missing):
  - Put GSE124472 10x MTX under 00_reference_umaps/inputs/fetal/tran_raw/
  - Rscript scripts/00_reference_umaps/upstream/run_devkidcc_tran.R
  - Then assemble into inputs/fetal/tran.h5ad (see upstream README)

Usage
-----
  python scripts/00_reference_umaps/upstream/process_fetal.py
  python scripts/00_reference_umaps/upstream/process_fetal.py --stage integrate
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib_config import ROOT  # noqa: E402

UPSTREAM = Path(__file__).resolve().parent
INTEGRATE = UPSTREAM / "integrate_fetal_harmony.py"
INPUT_DIR = ROOT / "00_reference_umaps/inputs/fetal"
OBJ = ROOT / "00_reference_umaps/objects/fetal_weeks_9_17_harmony_integrated.h5ad"


def _require_inputs() -> None:
    missing = [n for n in ("stewart", "hochane", "tran") if not (INPUT_DIR / f"{n}.h5ad").exists()]
    if missing:
        raise FileNotFoundError(
            "Missing local fetal inputs under "
            f"{INPUT_DIR}:\n  "
            + "\n  ".join(f"{m}.h5ad" for m in missing)
            + "\nPlace the three source .h5ad files there, then re-run."
        )


def run_integrate() -> None:
    _require_inputs()
    print(f"=== Fetal integrate (local → {ROOT / '00_reference_umaps'}) ===", flush=True)
    subprocess.run([sys.executable, "-u", str(INTEGRATE)], check=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--stage",
        choices=["integrate"],
        default="integrate",
        help="Only local integrate is supported (default: integrate)",
    )
    args = ap.parse_args()
    if args.stage == "integrate":
        run_integrate()
    if not OBJ.exists():
        raise FileNotFoundError(f"Expected integrated object missing: {OBJ}")
    print(f"Fetal upstream done → {OBJ}", flush=True)


if __name__ == "__main__":
    main()
