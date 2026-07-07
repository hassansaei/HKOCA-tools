"""
Figure: Dataset composition after harmonization — cell counts by age,
differentiation protocol, and single-cell protocol, per source file.

For each harmonized .h5ad file, this script:
  1. Reads obs columns (age, differentiation_protocol, sc_protocol) directly
     from the HDF5 store (no full AnnData load needed).
  2. Aggregates cell counts per (file, age, diff_protocol, sc_protocol).
  3. Renders a bubble/swarm plot:
       - y-axis: source file (one horizontal "zone" per file)
       - x-axis: age (days), continuous timeline
       - marker shape: sc protocol
       - marker color: differentiation protocol
       - marker size: number of cells in that (age, diff, sc) group

Output:
    02_harmonization/outputs/fig_harmonization_dataset_age_protocol_composition.png (dpi=300)

Run:
    python fig_harmonization_dataset_age_protocol_composition.py
"""

import re
import random
import colorsys
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Allow importing from common/ regardless of where this script is invoked from
sys.path.append(str(Path(__file__).resolve().parent.parent))
from common.data_paths import HARMONIZED_DATA_DIR, PAPER_ROOT
from common.plotting_style import save_figure, FIGURE_DPI

OUTPUT_PATH = PAPER_ROOT / "02_harmonization" / "outputs" / "fig_harmonization_dataset_age_protocol_composition.png"

# ---------------------------------------------------------
# 0. CONFIGURATION
# ---------------------------------------------------------
selected_meta_cols = {
    "differentiation_protocol": "diff_protocol",
    "sc_protocol": "sc_protocol",
    "age": "Age",
}

CUSTOM_COLORS = {
    "Morizane": "#87CEEB",
    "Takasato": "#FF0000",
    "Freedman": "#275000",
    "Hybrid_Takasato_Morizane": "#800080",
    "Hybrid_Takasato_Uchimura": "#FFA500",
}
MARKER_SEQ = ["^", "s", "D", "o", "v", "p", "*", "X", "h", "P"]
MIN_MARKER_SIZE = 40
MAX_MARKER_SIZE = 160
RANDOM_SEED = 42


# ---------------------------------------------------------
# 1. HDF5 HELPERS
# ---------------------------------------------------------
def clean_age_value(val):
    """Extract numeric digits from mixed age formats (e.g. 'd7' -> 7.0)."""
    if pd.isna(val) or str(val).strip() in ["nan", "None", "Unknown"]:
        return np.nan
    match = re.search(r"[-+]?\d*\.\d+|\d+", str(val).strip())
    return float(match.group()) if match else np.nan


def read_obs_column(obs_group, key):
    """
    Read a full obs column from an AnnData-on-disk HDF5 group, handling
    plain datasets, categorical groups, and nullable-array groups.
    """
    if key not in obs_group:
        return None
    item = obs_group[key]

    if isinstance(item, h5py.Group) and "codes" in item and "categories" in item:
        codes = item["codes"][:]
        categories = item["categories"][:]
        cats = [c.decode("utf-8") if isinstance(c, bytes) else c for c in categories]
        return np.array([cats[c] if 0 <= c < len(cats) else np.nan for c in codes], dtype=object)

    if isinstance(item, h5py.Group) and "data" in item:
        data_item = item["data"]
        arr = data_item[:] if isinstance(data_item, h5py.Dataset) else None
        if arr is None:
            return None
        if arr.dtype.kind in ["S", "O", "U"]:
            arr = np.array([v.decode("utf-8") if isinstance(v, bytes) else v for v in arr], dtype=object)
        if "mask" in item:
            mask = item["mask"][:]
            arr = np.where(mask, np.nan, arr)
        return arr

    if isinstance(item, h5py.Dataset):
        arr = item[:]
        if arr.dtype.kind in ["S", "O", "U"]:
            arr = np.array([v.decode("utf-8") if isinstance(v, bytes) else v for v in arr], dtype=object)
        return arr

    return None


def normalize_value_single(group: str, value: str) -> str:
    v = str(value).strip()
    key = re.sub(r"[\s\-]+", "_", v.lower())
    if group == "differentiation_protocol":
        mapping = {
            "morizane": "Morizane",
            "morizane_modified": "Morizane",
            "takasato": "Takasato",
            "takasato_2016": "Takasato",
            "takasato_modified": "Takasato",
            "hybrid_takasato_morizane": "Hybrid_Takasato_Morizane",
            "hybrid_takasato_uchimura": "Hybrid_Takasato_Uchimura",
        }
        return mapping.get(key, v)
    if group == "sc_protocol":
        mapping = {
            "10x_3_v3": "10x_3_v3.1",
            "10x_v2": "10x_3_v2",
            "drop_seq": "drop seq",
            "dropseq": "drop seq",
        }
        return mapping.get(key, v)
    return v


def generate_distinct_colors(n):
    """Generate n visually distinct colors spread evenly in hue space."""
    golden_ratio = 0.618033988749895
    colors = []
    for i in range(n):
        hue = (i * golden_ratio) % 1.0
        rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
        colors.append(f"#{int(rgb[0]*255):02x}{int(rgb[1]*255):02x}{int(rgb[2]*255):02x}")
    return colors


# ---------------------------------------------------------
# 2. LOAD & AGGREGATE
# ---------------------------------------------------------
def build_plot_dataframe(h5ad_files):
    plot_rows = []
    for path in h5ad_files:
        file_name = path.name
        try:
            with h5py.File(path, "r") as f:
                if "obs" not in f:
                    continue
                obs_g = f["obs"]

                age_arr = read_obs_column(obs_g, selected_meta_cols["age"])
                diff_arr = read_obs_column(obs_g, selected_meta_cols["differentiation_protocol"])
                sc_arr = read_obs_column(obs_g, selected_meta_cols["sc_protocol"])

                length = next((len(a) for a in [age_arr, diff_arr, sc_arr] if a is not None), 0)
                if length == 0:
                    continue

                if age_arr is None:
                    age_arr = np.full(length, np.nan, dtype=object)
                if diff_arr is None:
                    diff_arr = np.full(length, "Unknown", dtype=object)
                if sc_arr is None:
                    sc_arr = np.full(length, "Unknown", dtype=object)

                df_temp = pd.DataFrame({"age": age_arr, "diff": diff_arr, "sc": sc_arr}).fillna(np.nan)
                counts = df_temp.groupby(["age", "diff", "sc"], dropna=False).size().reset_index(name="cell_count")

                for _, row in counts.iterrows():
                    age_val = clean_age_value(row["age"])
                    if pd.isna(age_val):
                        continue
                    diff_val = row["diff"] if not pd.isna(row["diff"]) and str(row["diff"]) not in ["nan", "None"] else "Unknown"
                    sc_val = row["sc"] if not pd.isna(row["sc"]) and str(row["sc"]) not in ["nan", "None"] else "Unknown"

                    plot_rows.append({
                        "file": file_name,
                        "age_x": age_val,
                        "diff_protocol": normalize_value_single("differentiation_protocol", diff_val),
                        "sc_protocol": normalize_value_single("sc_protocol", sc_val),
                        "n_cells": int(row["cell_count"]),
                    })
        except Exception as e:
            print(f"Warning: Skipping file {file_name} due to read error: {e}")

    plot_df = pd.DataFrame(plot_rows)
    if plot_df.empty:
        raise ValueError("No valid data points found. Double check obs column names in selected_meta_cols.")
    return plot_df


# ---------------------------------------------------------
# 3. LAYOUT (per-file zones, jitter, marker sizing)
# ---------------------------------------------------------
def compute_file_zones(plot_df, file_order):
    file_y_info = {}
    current_bottom = 0.0
    for f in file_order:
        df_file = plot_df[plot_df["file"] == f]
        max_overlap = df_file.groupby("age_x").size().max() if not df_file.empty else 1
        height = max(1.0, 1.0 + (max_overlap - 1) * 0.4)
        center = current_bottom + (height / 2.0)
        top = current_bottom + height
        file_y_info[f] = {"center": center, "bottom": current_bottom, "top": top, "height": height}
        current_bottom = top
    return file_y_info, current_bottom


def apply_positions_and_jitter(plot_df, file_y_info):
    plot_df["y_pos"] = plot_df["file"].map(lambda x: file_y_info[x]["center"] if x in file_y_info else 0)
    for (f, age), group in plot_df.groupby(["file", "age_x"]):
        n = len(group)
        if n > 1:
            zone_height = file_y_info[f]["height"]
            spread = zone_height * 0.75
            offsets = np.array([-spread / 4, spread / 4]) if n == 2 else np.linspace(-spread / 2, spread / 2, n)
            plot_df.loc[group.index, "y_pos"] += offsets
    return plot_df


def assign_marker_sizes(plot_df):
    sizes = plot_df["n_cells"].astype(float)
    counts_nonzero = sizes[sizes > 0]
    size_levels = [float(round(s / 5) * 5) for s in np.linspace(MIN_MARKER_SIZE, MAX_MARKER_SIZE, 3)]

    if counts_nonzero.size > 0:
        q = np.quantile(counts_nonzero, [0.0, 1 / 3, 2 / 3, 1.0])
        bins = [float(q[0]), float(q[1]), float(q[2]), float(q[3])]
        for i in range(1, len(bins)):
            if bins[i] < bins[i - 1]:
                bins[i] = bins[i - 1]
        size_cat = pd.Series(np.digitize(sizes, bins[1:-1], right=True), index=plot_df.index)
        plot_df["marker_size"] = size_cat.map({0: size_levels[0], 1: size_levels[1], 2: size_levels[2]})
    else:
        bins = [0.0, 0.0, 0.0, 0.0]
        plot_df["marker_size"] = MIN_MARKER_SIZE
    return plot_df, bins, size_levels


def build_color_map(diff_values):
    random.seed(RANDOM_SEED)
    remaining = [p for p in diff_values if p not in CUSTOM_COLORS]
    colors = dict(CUSTOM_COLORS)
    if remaining:
        for proto, color in zip(remaining, generate_distinct_colors(len(remaining))):
            colors[proto] = color
    return {v: colors.get(v, "#808080") for v in diff_values}


# ---------------------------------------------------------
# 4. PLOT
# ---------------------------------------------------------
def make_figure(plot_df, file_order, file_y_info, current_bottom, bins, size_levels):
    sc_values = sorted(plot_df["sc_protocol"].astype(str).unique())
    diff_values = sorted(plot_df["diff_protocol"].astype(str).unique())
    marker_map = {v: MARKER_SEQ[i % len(MARKER_SEQ)] for i, v in enumerate(sc_values)}
    color_map = build_color_map(diff_values)

    numeric_vals = sorted(plot_df["age_x"].unique())
    tickvals = numeric_vals
    ticktext = [f"d{int(v)}" if float(v).is_integer() else f"d{v}" for v in numeric_vals]

    fig, ax = plt.subplots(figsize=(18, 10))
    for sc in sc_values:
        for diff in diff_values:
            sub = plot_df[(plot_df["sc_protocol"] == sc) & (plot_df["diff_protocol"] == diff)]
            if sub.empty:
                continue
            ax.scatter(
                sub["age_x"], sub["y_pos"], s=sub["marker_size"],
                marker=marker_map[sc], color=color_map[diff],
                edgecolor="black", linewidth=0.3, alpha=0.7, zorder=3,
            )

    legend_elements = [Line2D([0], [0], linestyle="None", color="none", label="sc protocol (shape):")]
    for sc in sc_values:
        legend_elements.append(Line2D([0], [0], marker=marker_map[sc], linestyle="None",
                                       color="black", markeredgecolor="black", markersize=8, label=f"  {sc}"))
    legend_elements.append(Line2D([0], [0], linestyle="None", color="none", label="diff protocol (color):"))
    for diff in diff_values:
        legend_elements.append(Line2D([0], [0], marker="o", linestyle="None",
                                       color=color_map[diff], markeredgecolor="black", markersize=8, label=f"  {diff}"))
    legend_elements.append(Line2D([0], [0], linestyle="None", color="none", label="cells (size):"))
    for i in range(3):
        label = f"{int(bins[i]):,}-{int(bins[i+1]):,}"
        legend_elements.append(Line2D([0], [0], marker="o", linestyle="None", color="gray",
                                       markeredgecolor="black", markersize=max(4, np.sqrt(size_levels[i])), label=f"  {label}"))

    y_ticks = [file_y_info[f]["center"] for f in file_order if f in file_y_info]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f for f in file_order if f in file_y_info])

    for i, f in enumerate(file_order):
        if f not in file_y_info:
            continue
        info = file_y_info[f]
        zone_color = "#ffffff" if i % 2 == 0 else "#f5f5f5"
        ax.axhspan(info["bottom"], info["top"], color=zone_color, zorder=0)
        ax.axhline(info["bottom"], color="#cccccc", linewidth=1, zorder=1)
    ax.axhline(current_bottom, color="#cccccc", linewidth=1, zorder=1)

    ax.set_xlabel("Age")
    ax.set_ylabel("Dataset (file)")
    ax.set_xticks(tickvals)
    ax.set_xticklabels(ticktext, rotation=45, ha="right")
    ax.set_xlim(min(tickvals) - 1, max(tickvals) + 1)
    ax.set_ylim(0, current_bottom)
    ax.margins(y=0)
    ax.grid(axis="x", linestyle="--", alpha=0.5, zorder=2)
    ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=True, fontsize=9)

    fig.tight_layout(rect=[0, 0, 0.8, 1])
    return fig


# ---------------------------------------------------------
# 5. MAIN
# ---------------------------------------------------------
def main():
    h5ad_files = sorted(HARMONIZED_DATA_DIR.glob("*.h5ad"))
    if not h5ad_files:
        raise FileNotFoundError(f"No .h5ad files found in {HARMONIZED_DATA_DIR}")

    print(f"Scanning {len(h5ad_files)} harmonized files...")
    plot_df = build_plot_dataframe(h5ad_files)

    file_order = sorted(plot_df["file"].unique())
    file_y_info, current_bottom = compute_file_zones(plot_df, file_order)
    plot_df = apply_positions_and_jitter(plot_df, file_y_info)
    plot_df, bins, size_levels = assign_marker_sizes(plot_df)

    fig = make_figure(plot_df, file_order, file_y_info, current_bottom, bins, size_levels)
    save_figure(fig, OUTPUT_PATH, dpi=FIGURE_DPI)
    plt.close(fig)


if __name__ == "__main__":
    main()
