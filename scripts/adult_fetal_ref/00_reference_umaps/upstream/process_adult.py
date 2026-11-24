#!/usr/bin/env python
"""Upstream adult: remove diseased cells (+ optional organoid projection rebuild).

Filter
------
  disease == "normal"
  sampled_site_condition ∈ {"healthy", "adjacent"}

Default: audit filter → tables + figure under 00_reference_umaps/
Optional: --rebuild-projection runs the full HCKA healthy+adjacent scANVI /
          organoid projection pipeline (heavy; already done).

Usage
-----
  python scripts/00_reference_umaps/upstream/process_adult.py
  python scripts/00_reference_umaps/upstream/process_adult.py --rebuild-projection
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib_config import DIRS, WORKSPACE, apply_fig_style, ensure_dirs, read_h5_obs_col  # noqa: E402

apply_fig_style()

HCKA_DIR = WORKSPACE / "human_kidney_reference_projection"
HCKA_FULL = HCKA_DIR / "Single_cell_human_kidney_cell_atlas.h5ad"
HCKA_HEALTHY = HCKA_DIR / "01_reference" / "HCKA_healthy_adjacent_reference_scANVI_latent.h5ad"
PROJECTION_PY = HCKA_DIR / "03_hcka_organoid_projection.py"
ATLAS_PY = Path("/data-master/pure-workspace/hmami/conda_envs/atlas/bin/python")


def load_meta(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        disease = read_h5_obs_col(f, "disease")
        site = np.array(
            [str(x).strip().lower() for x in read_h5_obs_col(f, "sampled_site_condition")],
            dtype=object,
        )
        broad = read_h5_obs_col(f, "broad_annotation")
    return pd.DataFrame(
        {
            "disease": np.array([str(x).strip().lower() for x in disease], dtype=object),
            "sampled_site_condition": site,
            "broad_annotation": broad,
        }
    )


def audit_filter() -> None:
    ensure_dirs()
    fig_dir = DIRS["00_fig"] / "adult_processing"
    tab_dir = DIRS["00_tab"]
    fig_dir.mkdir(parents=True, exist_ok=True)

    print("=== Adult: remove diseased cells ===", flush=True)
    meta = load_meta(HCKA_FULL)
    n_all = len(meta)
    keep = (meta["disease"] == "normal") & meta["sampled_site_condition"].isin(
        ["healthy", "adjacent"]
    )
    n_keep, n_drop = int(keep.sum()), n_all - int(keep.sum())

    (
        meta.assign(kept=np.where(keep, "kept", "excluded"))
        .groupby(["disease", "sampled_site_condition", "kept"], as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
        .sort_values("n_cells", ascending=False)
        .to_csv(tab_dir / "adult_filter_by_disease_site.csv", index=False)
    )

    summary = pd.DataFrame(
        [
            {"category": "all_HCKA", "n_cells": n_all},
            {"category": "kept_healthy_adjacent_normal", "n_cells": n_keep},
            {"category": "excluded_diseased_or_other", "n_cells": n_drop},
            {
                "category": "kept_healthy",
                "n_cells": int(
                    ((meta["disease"] == "normal") & (meta["sampled_site_condition"] == "healthy")).sum()
                ),
            },
            {
                "category": "kept_adjacent",
                "n_cells": int(
                    ((meta["disease"] == "normal") & (meta["sampled_site_condition"] == "adjacent")).sum()
                ),
            },
        ]
    )
    summary.to_csv(tab_dir / "adult_filter_summary.csv", index=False)
    print(summary.to_string(index=False), flush=True)

    broad_kept = (
        meta.loc[keep, "broad_annotation"]
        .value_counts()
        .rename_axis("broad_annotation")
        .reset_index(name="n_cells")
    )
    broad_kept["pct"] = 100.0 * broad_kept["n_cells"] / max(n_keep, 1)
    broad_kept.to_csv(tab_dir / "adult_healthy_adjacent_broad_counts.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 6.0))
    ax = axes[0]
    ax.bar(
        [f"Kept\n(healthy+adjacent)\nn={n_keep:,}", f"Excluded\n(diseased/other)\nn={n_drop:,}"],
        [n_keep, n_drop],
        color=["#1565C0", "#BDBDBD"],
        width=0.65,
    )
    ax.set_ylabel("Cells", fontsize=14)
    ax.set_title("Adult HCKA: remove diseased cells", fontsize=16, pad=10)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    site_counts = (
        meta.loc[keep]
        .groupby("sampled_site_condition")
        .size()
        .reindex(["healthy", "adjacent"])
        .fillna(0)
        .astype(int)
    )
    ax.bar(["Healthy", "Adjacent"], site_counts.to_numpy(), color=["#2E7D32", "#FF8F00"], width=0.65)
    for i, v in enumerate(site_counts.to_numpy()):
        ax.text(i, v + max(site_counts.max() * 0.01, 500), f"{v:,}", ha="center", fontsize=12)
    ax.set_ylabel("Cells", fontsize=14)
    ax.set_title("Kept reference by tissue condition", fontsize=16, pad=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out_png = fig_dir / "adult_filter_healthy_adjacent.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}", flush=True)

    prov = {
        "filter": "disease==normal & sampled_site_condition in {healthy, adjacent}",
        "n_all": n_all,
        "n_kept": n_keep,
        "n_excluded": n_drop,
        "healthy_adjacent_object": str(HCKA_HEALTHY),
        "projection_rebuild": str(PROJECTION_PY),
    }
    (tab_dir / "adult_reference_provenance.json").write_text(json.dumps(prov, indent=2))


def rebuild_projection() -> None:
    if not PROJECTION_PY.exists():
        raise FileNotFoundError(PROJECTION_PY)
    py = ATLAS_PY if ATLAS_PY.exists() else Path(sys.executable)
    print(f"=== Adult: rebuild HCKA projection ({PROJECTION_PY.name}) ===", flush=True)
    subprocess.run([str(py), "-u", str(PROJECTION_PY)], check=True, cwd=str(HCKA_DIR))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--rebuild-projection",
        action="store_true",
        help="Also re-run full HCKA healthy+adjacent + organoid projection (slow)",
    )
    args = ap.parse_args()
    audit_filter()
    if args.rebuild_projection:
        rebuild_projection()
    print("Adult upstream done.", flush=True)


if __name__ == "__main__":
    main()
