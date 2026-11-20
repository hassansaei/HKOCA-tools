#!/usr/bin/env python3
"""Aggregate iteration metrics and write scorecard figures.

Reads subsample_iter_*/scib_raw_metrics_full.csv and writes mean/std tables
plus PNG scorecards under results/ (and checkpoints/ for continuity).
"""
from __future__ import annotations

import os

import matplotlib
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCH_DIR = os.path.dirname(SCRIPT_DIR)
with open(os.path.join(BENCH_DIR, "config.yaml"), encoding="utf-8") as _fh:
    _CFG = yaml.safe_load(_fh)

CHECKPOINTS = os.path.join(BENCH_DIR, "checkpoints")
RESULTS = os.path.join(BENCH_DIR, "results")
N_ITERS = int(_CFG["benchmark"]["n_iters"])
BIO_KEYS = ["ARI", "NMI", "ASW_celltype", "iso_label_F1", "iso_label_ASW", "cLISI"]
BATCH_KEYS = ["ASW_batch", "graph_connectivity", "iLISI", "kBET", "PCR"]



def save_both(fig, name: str, dpi: int = 150) -> None:
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(CHECKPOINTS, exist_ok=True)
    for outdir in (CHECKPOINTS, RESULTS):
        fig.savefig(os.path.join(outdir, name), dpi=dpi, bbox_inches="tight")


print("Aggregating Level_2 benchmark results.")

all_records = []
for i in range(1, N_ITERS + 1):
    csv_path = os.path.join(CHECKPOINTS, f"subsample_iter_{i}", "scib_raw_metrics_full.csv")
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        all_records.append(df)
        print(f"  loaded iter {i}: {len(df)} methods")
    else:
        print(f"  MISSING iter {i}: {csv_path}")

if not all_records:
    raise SystemExit("No scIB CSV files found — abort viz")

iter_df = pd.concat(all_records, ignore_index=True)
mean_df = iter_df.groupby("method").mean(numeric_only=True)
std_df = iter_df.groupby("method").std(numeric_only=True)

os.makedirs(RESULTS, exist_ok=True)
for outdir in (CHECKPOINTS, RESULTS):
    mean_df.to_csv(os.path.join(outdir, "scib_mean_metrics.csv"))
    std_df.to_csv(os.path.join(outdir, "scib_std_metrics.csv"))

print("\nWilcoxon P-values...")
winner_method = mean_df["overall"].idxmax()
sorted_methods = mean_df["overall"].sort_values(ascending=False).index.tolist()
competitors = [m for m in sorted_methods if m not in [winner_method, "pca"]]
runner_up = competitors[0] if competitors else None

comparisons = [(winner_method, "pca")]
if runner_up:
    comparisons.append((winner_method, runner_up))

for method_a, method_b in comparisons:
    df_a = iter_df[iter_df["method"] == method_a].set_index("seed")["overall"]
    df_b = iter_df[iter_df["method"] == method_b].set_index("seed")["overall"]
    common = df_a.index.intersection(df_b.index)
    if len(common) > 0:
        W, p = stats.wilcoxon(df_a.loc[common], df_b.loc[common], alternative="greater")
        print(f"  {method_a} > {method_b}: p={p:.4f} (W={W})")

print("\nGenerating figures...")
all_cols = BIO_KEYS + BATCH_KEYS

fig, axes = plt.subplots(1, 2, figsize=(24, 7))
sns.heatmap(mean_df[all_cols], ax=axes[0], annot=True, fmt=".3f", cmap="YlOrRd", vmin=0, vmax=1, linewidths=0.4)
axes[0].set_title(f"Mean scIB scores ({N_ITERS} iterations)", fontsize=25, fontweight="bold")
axes[0].tick_params(axis="y", labelsize=20)
axes[0].tick_params(axis="x", labelsize=20)
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha="right", fontsize=20)
axes[0].set_yticklabels(axes[0].get_yticklabels(), rotation=0, fontsize=20)
for text in axes[0].texts:
    text.set_fontsize(20)

sns.heatmap(std_df[all_cols], ax=axes[1], annot=True, fmt=".3f", cmap="Blues", linewidths=0.4)
axes[1].set_title("SD (uncertainty)", fontsize=25, fontweight="bold")
axes[1].tick_params(axis="y", labelsize=25)
axes[1].tick_params(axis="x", labelsize=25)
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha="right", fontsize=20)
axes[1].set_yticklabels(axes[1].get_yticklabels(), rotation=0, fontsize=20)
for text in axes[1].texts:
    text.set_fontsize(20)

plt.tight_layout()
save_both(fig, "scib_mean_sd_heatmap_size.png", dpi=150)
plt.close()
print("  Saved: scib_mean_sd_heatmap_size.png")

fig, ax = plt.subplots(figsize=(14, 10))
ax.bar(
    sorted_methods,
    mean_df.loc[sorted_methods, "overall"],
    yerr=std_df.loc[sorted_methods, "overall"],
    capsize=5,
    color="steelblue",
    edgecolor="white",
    width=0.6,
)
ax.set_ylabel("scIB score  (mean ± SD)", fontsize=35, fontweight="bold")
ax.set_title("Integration benchmarking : overall score (Level_2)", fontsize=28, fontweight="bold")
ax.set_ylim(0, 1)
ax.tick_params(axis="y", labelsize=35)
plt.xticks(rotation=35, ha="right", fontsize=35)
plt.tight_layout()
save_both(fig, "scib_overall_errorbars_size.png", dpi=300)
plt.close()
print("  Saved: scib_overall_errorbars_size.png")

plot_df = mean_df.copy()
plot_df["Bio Score (60%)"] = mean_df[BIO_KEYS].mean(axis=1)
plot_df["Batch Score (40%)"] = mean_df[BATCH_KEYS].mean(axis=1)
plot_df = plot_df.rename(columns={"overall": "Overall Score"}).sort_values("Overall Score", ascending=False)

plot_norm = plot_df.copy()
for col in BIO_KEYS + BATCH_KEYS + ["Bio Score (60%)", "Batch Score (40%)", "Overall Score"]:
    cmin, cmax = plot_df[col].min(), plot_df[col].max()
    plot_norm[col] = (plot_df[col] - cmin) / (cmax - cmin) if cmax - cmin > 1e-9 else 1.0

fig, axes = plt.subplots(1, 4, figsize=(26, 9), gridspec_kw={"width_ratios": [6, 5, 2, 1.2]})
plt.subplots_adjust(wspace=0.08)
kws = {"fmt": ".3f", "linewidths": 0.5, "annot_kws": {"size": 20}}

sns.heatmap(data=plot_norm[BIO_KEYS], annot=plot_df[BIO_KEYS], cmap="Blues", ax=axes[0], cbar=False, **kws)
axes[0].set_title("Biological Conservation", fontsize=25, fontweight="bold", pad=15)
axes[0].set_ylabel("Integration Method", fontsize=25, fontweight="bold")
axes[0].tick_params(axis="y", rotation=0, labelsize=20)
axes[0].tick_params(axis="x", labelsize=20)
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha="right", fontsize=20)

sns.heatmap(data=plot_norm[BATCH_KEYS], annot=plot_df[BATCH_KEYS], cmap="Oranges", ax=axes[1], cbar=False, **kws)
axes[1].set_title("Batch Correction", fontsize=25, fontweight="bold", pad=15)
axes[1].set_yticks([])
axes[1].set_ylabel("")
axes[1].tick_params(axis="x", labelsize=20)
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha="right", fontsize=20)

sns.heatmap(
    data=plot_norm[["Bio Score (60%)", "Batch Score (40%)"]],
    annot=plot_df[["Bio Score (60%)", "Batch Score (40%)"]],
    cmap="Greens",
    ax=axes[2],
    cbar=False,
    **kws,
)
axes[2].set_title("Category Scores", fontsize=25, fontweight="bold", pad=15)
axes[2].set_yticks([])
axes[2].set_ylabel("")
axes[2].tick_params(axis="x", labelsize=20)
axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, ha="right", fontsize=20)

sns.heatmap(
    data=plot_norm[["Overall Score"]],
    annot=plot_df[["Overall Score"]],
    cmap="Reds",
    ax=axes[3],
    cbar=True,
    **kws,
)
axes[3].set_title("Overall", fontsize=25, fontweight="bold", pad=15)
axes[3].set_yticks([])
axes[3].set_ylabel("")
axes[3].tick_params(axis="x", labelsize=20)
axes[3].set_xticklabels(axes[3].get_xticklabels(), rotation=45, ha="right", fontsize=20)

for ax in axes:
    for text in ax.texts:
        text.set_fontsize(20)

plt.suptitle("HKOCA Integration Benchmark — Level_2 Master Scorecard", fontsize=32, fontweight="bold", y=1.03)
save_both(fig, "MASTER_SCORECARD_size.png", dpi=300)
plt.close()
print("  Saved: MASTER_SCORECARD_size.png")
print("\nAll visualizations generated successfully.")
