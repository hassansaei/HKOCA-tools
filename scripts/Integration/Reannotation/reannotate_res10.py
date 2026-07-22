# Level 3 reannotation at Leiden resolution 10
# Configuration: reannotation_config.yaml

import os
import sys
import subprocess
from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
        config_file = Path(__file__).resolve().parent / "reannotation_config.yaml"
    with open(config_file, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


CFG = load_config()
PATHS = CFG["paths"]
CLUSTER_CFG = CFG["clustering"]
ANNOT_CFG = CFG["annotation"]
MARKER_CFG = CFG["marker_refinement"]
PLOT_CFG = CFG["plotting"]

INPUT_H5AD = PATHS["input_h5ad"]
OUTPUT_H5AD = PATHS["output_h5ad"]
FIGURE_DIR = PATHS["figure_dir"]
os.makedirs(FIGURE_DIR, exist_ok=True)

print("Loading full integrated atlas into RAM...")
adata = sc.read_h5ad(INPUT_H5AD)

print("\nComputing neighborhood graph on X_emb...")
sc.pp.neighbors(
    adata,
    use_rep=CLUSTER_CFG["use_rep"],
    n_neighbors=CLUSTER_CFG["n_neighbors"],
    random_state=CLUSTER_CFG["random_state"],
)

print(f"Running high-resolution Leiden clustering (resolution = {CLUSTER_CFG['resolution']})...")
sc.tl.leiden(
    adata,
    resolution=CLUSTER_CFG["resolution"],
    key_added=CLUSTER_CFG["leiden_key"],
    random_state=CLUSTER_CFG["random_state"],
)
num_clusters = adata.obs[CLUSTER_CFG["leiden_key"]].nunique()
print(f"Generated {num_clusters} high-resolution micro-clusters.")

print("\nExecuting dominant cell type assignment...")
l3_to_l2_map = adata.obs.groupby(ANNOT_CFG["source_level3"])[ANNOT_CFG["source_level2"]].agg(
    lambda x: pd.Series.mode(x)[0]
).to_dict()
l3_to_l1_map = adata.obs.groupby(ANNOT_CFG["source_level3"])[ANNOT_CFG["source_level1"]].agg(
    lambda x: pd.Series.mode(x)[0]
).to_dict()

cluster_to_majority_l3 = {}
for cluster in adata.obs[CLUSTER_CFG["leiden_key"]].cat.categories:
    cluster_cells = adata.obs[adata.obs[CLUSTER_CFG["leiden_key"]] == cluster]
    majority_vote = cluster_cells[ANNOT_CFG["source_level3"]].mode()[0]
    cluster_to_majority_l3[cluster] = majority_vote

adata.obs[ANNOT_CFG["output_level3"]] = (
    adata.obs[CLUSTER_CFG["leiden_key"]].map(cluster_to_majority_l3).astype("category")
)
adata.obs[ANNOT_CFG["output_level2"]] = (
    adata.obs[ANNOT_CFG["output_level3"]].map(l3_to_l2_map).astype("category")
)
adata.obs[ANNOT_CFG["output_level1"]] = (
    adata.obs[ANNOT_CFG["output_level3"]].map(l3_to_l1_map).astype("category")
)

changed_cells = (
    adata.obs[ANNOT_CFG["source_level3"]].astype(str)
    != adata.obs[ANNOT_CFG["output_level3"]].astype(str)
).sum()
print(
    f"Reannotation complete. {changed_cells:,} cells "
    f"({(changed_cells / adata.n_obs) * 100:.1f}%) were relabeled."
)

print("\n" + "=" * 60)
print("Starting marker gene refinement (100k subset)")
print("=" * 60)

print(f"Extracting {MARKER_CFG['n_subset']:,} random cells...")
adata_sub = sc.pp.subsample(
    adata,
    n_obs=MARKER_CFG["n_subset"],
    copy=True,
    random_state=MARKER_CFG["random_state"],
)

if adata_sub.raw is not None:
    adata_sub.X = adata_sub.raw.X

print(f"Running Wilcoxon rank-sum test for DEGs across {MARKER_CFG['groupby']}...")
sc.tl.rank_genes_groups(
    adata_sub,
    groupby=MARKER_CFG["groupby"],
    method=MARKER_CFG["method"],
    use_raw=MARKER_CFG["use_raw"],
    pts=MARKER_CFG["pts"],
)

top_10_genes_dict = {}
for cluster in adata_sub.obs[MARKER_CFG["groupby"]].cat.categories:
    genes = sc.get.rank_genes_groups_df(adata_sub, group=cluster)["names"].head(
        MARKER_CFG["n_top_genes"]
    ).tolist()
    top_10_genes_dict[cluster] = genes

print("Extracted top 10 marker genes for each cell type.")

print("Computing hierarchical clustering (dendrogram)...")
sc.tl.dendrogram(adata_sub, groupby=MARKER_CFG["groupby"])

print("Rendering hierarchical marker matrix plot...")
sc.set_figure_params(dpi=PLOT_CFG["dpi"], facecolor=PLOT_CFG["facecolor"])
sc.pl.rank_genes_groups_matrixplot(
    adata_sub,
    n_genes=MARKER_CFG["n_top_genes"],
    use_raw=MARKER_CFG["use_raw"],
    dendrogram=True,
    cmap=PLOT_CFG["cmap"],
    standard_scale=PLOT_CFG["standard_scale"],
    title=PLOT_CFG["title"],
    show=False,
)
plot_path = os.path.join(FIGURE_DIR, PLOT_CFG["matrixplot_filename"])
plt.savefig(plot_path, dpi=PLOT_CFG["dpi"], bbox_inches="tight")
plt.close()

print(f"Marker gene plot saved to: {plot_path}")

print("\nSaving final reannotated atlas to disk...")
adata.write_h5ad(OUTPUT_H5AD)
print("Pipeline complete.")
