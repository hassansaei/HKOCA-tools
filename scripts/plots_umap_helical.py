import anndata as ad
import pickle
import numpy as np
import umap
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from pathlib import Path
Path("/out/plots/umap/").mkdir(parents=True, exist_ok=True)


def plot_umap(adata, mapper, legend, model):
    plot_df = pd.DataFrame(mapper.embedding_, columns=['px', 'py'])
    labels = adata.obs[legend]
    plot_df['label'] = labels.values

    # Set up a figure with specific size
    fig, ax = plt.subplots(figsize=(10, 10))  # Square figure

    # Create the scatter plot without automatic legend
    sns.scatterplot(
        data=plot_df,
        x='px', y='py',
        hue='label',
        s=0.5,
        palette="pastel",
        legend=False,
        ax=ax
    )

    # Create a manual legend
    unique_labels = plot_df['label'].unique()
    palette = sns.color_palette("pastel", n_colors=len(unique_labels))
    label_to_color = dict(zip(unique_labels, palette))
    handles = [mpatches.Patch(color=label_to_color[label], label=label) for label in unique_labels]

    # Add the legend outside the plot
    ax.legend(
        handles=handles,
        title=f"{legend}",
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        borderaxespad=0.,
        frameon=False
    )

    # Make sure the plot area is perfectly square
    ax.set_aspect('equal', adjustable='datalim')

    # Tight layout to avoid cutting off labels
    plt.tight_layout()

    # Save the figure
    plt.savefig(f"/out/plots/umap/umap_{model}_{legend}.png", bbox_inches='tight')

if __name__ == "__main__":
    # Load the data
    file_path = "/mnt/data/processed/HKOCA_concatenated_harmonized_processed_filtered_raw.h5ad"
    adata = ad.read_h5ad(file_path)
    
    for model in ["geneformer_v2", "scGPT", "UCE"]:

        with open(f"/out/data/processed/{model}_umap.pkl", 'rb') as file:
            mapper = pickle.load(file)

        for legend in ['diff_protocol', 'study', 'sequencing', 'sc_protocol', 'Age']:
            plot_umap(adata, mapper, legend, model)