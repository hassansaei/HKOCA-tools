import os
import sys
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats

try:
    import yaml
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyyaml", "-q"])
    import yaml


def load_config():
    argv = sys.argv[1:]
    if "--config" in argv:
        config_file = Path(argv[argv.index("--config") + 1]).resolve()
    else:
        config_file = Path(__file__).resolve().parent / "benchmark_config.yaml"
    with open(config_file, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


CFG = load_config()

VIZ_CFG = CFG["visualization"]
METRICS_CFG = CFG["metrics"]
CHECKPOINT_BASE = CFG["paths"]["checkpoints_dir"]
N_ITERATIONS = VIZ_CFG["n_iterations"]
BASELINE_METHOD = VIZ_CFG["baseline_method"]
BIO_KEYS = METRICS_CFG["bio_keys"]
BATCH_KEYS = METRICS_CFG["batch_keys"]
BIO_WEIGHT = METRICS_CFG["bio_weight"]
BATCH_WEIGHT = METRICS_CFG["batch_weight"]

print("Aggregating benchmark results.")

all_records = []
for i in range(1, N_ITERATIONS + 1):
    csv_path = os.path.join(
        CHECKPOINT_BASE,
        CFG["iteration"]["checkpoint_subdir_pattern"].format(iter_id=i),
        METRICS_CFG["output_file"],
    )
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        all_records.append(df)

if not all_records:
    print("No scIB CSV files found.")
    exit()

iter_df = pd.concat(all_records, ignore_index=True)
mean_df = iter_df.groupby("method").mean()
std_df = iter_df.groupby("method").std()

mean_df.to_csv(os.path.join(CHECKPOINT_BASE, VIZ_CFG["mean_metrics_file"]))
std_df.to_csv(os.path.join(CHECKPOINT_BASE, VIZ_CFG["std_metrics_file"]))

print("\nCalculating Wilcoxon P-values...")
winner_method = mean_df["overall"].idxmax()
sorted_methods = mean_df["overall"].sort_values(ascending=False).index.tolist()
competitors = [m for m in sorted_methods if m not in [winner_method, BASELINE_METHOD]]
runner_up = competitors[0] if competitors else None

COMPARISONS = [(winner_method, BASELINE_METHOD)]
if runner_up: COMPARISONS.append((winner_method, runner_up))

for method_a, method_b in COMPARISONS:
    df_a = iter_df[iter_df["method"] == method_a].set_index("seed")["overall"]
    df_b = iter_df[iter_df["method"] == method_b].set_index("seed")["overall"]
    common_seeds = df_a.index.intersection(df_b.index)
    if len(common_seeds) > 0:
        W, p = stats.wilcoxon(df_a.loc[common_seeds], df_b.loc[common_seeds], alternative="greater")
        print(f"  {method_a} > {method_b}: p-value = {p:.4f} (W={W})")

print("\nGenerating visualizations...")
bio_metrics = BIO_KEYS
batch_metrics = BATCH_KEYS
ALL_COLS = bio_metrics + batch_metrics

heatmap_figsize = tuple(VIZ_CFG["heatmap_figsize"])
fig, axes = plt.subplots(1, 2, figsize=heatmap_figsize)
sns.heatmap(mean_df[ALL_COLS], ax=axes[0], annot=True, fmt=".3f", cmap="YlOrRd", vmin=0, vmax=1, linewidths=0.4)
axes[0].set_title(f"Mean scIB scores ({N_ITERATIONS} iterations)", fontsize=13, fontweight="bold")
sns.heatmap(std_df[ALL_COLS], ax=axes[1], annot=True, fmt=".3f", cmap="Blues", linewidths=0.4)
axes[1].set_title("SD (uncertainty)", fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(os.path.join(CHECKPOINT_BASE, VIZ_CFG["heatmap_file"]), dpi=VIZ_CFG["heatmap_dpi"], bbox_inches="tight")
plt.close()
print(f"  Saved: {VIZ_CFG['heatmap_file']}")

errorbars_figsize = tuple(VIZ_CFG["errorbars_figsize"])
fig, ax = plt.subplots(figsize=errorbars_figsize)
ax.bar(sorted_methods, mean_df.loc[sorted_methods, "overall"], yerr=std_df.loc[sorted_methods, "overall"], capsize=5, color="steelblue", edgecolor="white", width=0.6)
ax.set_ylabel("Overall scIB score  (mean ± SD)", fontsize=11)
ax.set_title("Integration benchmarking — overall score with uncertainty", fontsize=13)
ax.set_ylim(0, 1)
plt.xticks(rotation=35, ha="right", fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(CHECKPOINT_BASE, VIZ_CFG["errorbars_file"]), dpi=VIZ_CFG["errorbars_dpi"], bbox_inches="tight")
plt.close()
print(f"  Saved: {VIZ_CFG['errorbars_file']}")

plot_df = mean_df.copy()
plot_df[f"Bio Score ({int(BIO_WEIGHT * 100)}%)"] = mean_df[bio_metrics].mean(axis=1)
plot_df[f"Batch Score ({int(BATCH_WEIGHT * 100)}%)"] = mean_df[batch_metrics].mean(axis=1)
plot_df = plot_df.rename(columns={"overall": "Overall Score"}).sort_values("Overall Score", ascending=False)

bio_col = f"Bio Score ({int(BIO_WEIGHT * 100)}%)"
batch_col = f"Batch Score ({int(BATCH_WEIGHT * 100)}%)"
score_cols = [bio_col, batch_col, "Overall Score"]

plot_norm = plot_df.copy()
for col in bio_metrics + batch_metrics + score_cols:
    cmin, cmax = plot_df[col].min(), plot_df[col].max()
    plot_norm[col] = (plot_df[col] - cmin) / (cmax - cmin) if cmax - cmin > 1e-9 else 1.0

scorecard_figsize = tuple(VIZ_CFG["scorecard_figsize"])
fig, axes = plt.subplots(1, 4, figsize=scorecard_figsize, gridspec_kw={'width_ratios': VIZ_CFG["scorecard_width_ratios"]})
plt.subplots_adjust(wspace=0.08)
kws = {"fmt": ".3f", "linewidths": 0.5, "annot_kws": {"size": 11}}

sns.heatmap(data=plot_norm[bio_metrics], annot=plot_df[bio_metrics], cmap="Blues", ax=axes[0], cbar=False, **kws)
axes[0].set_title("Biological Conservation", fontsize=15, fontweight="bold", pad=15)
axes[0].set_ylabel("Integration Method", fontsize=14, fontweight="bold")
axes[0].tick_params(axis='y', rotation=0, labelsize=12)
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha='right')

sns.heatmap(data=plot_norm[batch_metrics], annot=plot_df[batch_metrics], cmap="Oranges", ax=axes[1], cbar=False, **kws)
axes[1].set_title("Batch Correction", fontsize=15, fontweight="bold", pad=15)
axes[1].set_yticks([]); axes[1].set_ylabel("")
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha='right')

sns.heatmap(data=plot_norm[[bio_col, batch_col]], annot=plot_df[[bio_col, batch_col]], cmap="Greens", ax=axes[2], cbar=False, **kws)
axes[2].set_title("Category Scores", fontsize=15, fontweight="bold", pad=15)
axes[2].set_yticks([]); axes[2].set_ylabel("")
axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, ha='right')

sns.heatmap(data=plot_norm[["Overall Score"]], annot=plot_df[["Overall Score"]], cmap="Reds", ax=axes[3], cbar=True, **kws)
axes[3].set_title("Overall", fontsize=15, fontweight="bold", pad=15)
axes[3].set_yticks([]); axes[3].set_ylabel("")

plt.suptitle("HKOCA Integration Benchmark Master Scorecard", fontsize=22, fontweight="bold", y=1.03)
fig.savefig(os.path.join(CHECKPOINT_BASE, VIZ_CFG["scorecard_file"]), dpi=VIZ_CFG["scorecard_dpi"], bbox_inches="tight")
plt.close()
print(f"  Saved: {VIZ_CFG['scorecard_file']}")

print("\nAll visualizations generated successfully.")
