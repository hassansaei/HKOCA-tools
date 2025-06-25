import scanpy as sc
import anndata as ad

import os
import pandas as pd
import numpy as np

import scvi
import torch
torch.set_float32_matmul_precision("medium")
import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt

def combine_annotations(row):
    if pd.isna(row['level_3']):
        return row['level_2']
    else:
        return row['level_3']


output_dir = "/out/plots/integration_comparison"
os.makedirs(output_dir, exist_ok=True)


file_path = "/out/data/HKOCA_concatenated_harmonized_processed_filtered_raw_annotated_snapseed-indiv.h5ad"
adata = ad.read_h5ad(file_path)
adata.obs['celltype'] = adata.obs.apply(combine_annotations, axis=1)

meta = ['celltype', 'diff_protocol', 'study', 'sequencing', 'sc_protocol', 'Age']


#load scvi latent representation
scvi_latent = ad.read_h5ad("/out/data/scvi/scvi_latent.h5ad")
adata.obsm["X_scvi"] = scvi_latent.X
sc.pp.neighbors(adata, use_rep="X_scvi")
sc.tl.umap(adata, min_dist = 0.3)
for col in meta:
    sc.pl.umap(adata, show=False, color=col)
    plt.savefig(f'{output_dir}/umap_scvi_{col}.png', bbox_inches='tight')

#load scanvi latent representation
scanvi_latent = ad.read_h5ad("/out/data/scvi/scanvi_latent.h5ad")
adata.obsm["X_scanvi"] = scanvi_latent.X
sc.pp.neighbors(adata, use_rep="X_scanvi")
sc.tl.umap(adata, min_dist = 0.3)
for col in meta:
    sc.pl.umap(adata, show=False, color=col)
    plt.savefig(f'{output_dir}/umap_scanvi_{col}.png', bbox_inches='tight')

#load scPoli latent representation
scanvi_latent = ad.read_h5ad("/out/data/scpoli/scpoli_latent.h5ad")
adata.obsm["X_scpoli"] = scanvi_latent.X
sc.pp.neighbors(adata, use_rep="X_scpoli")
sc.tl.umap(adata, min_dist = 0.3)
for col in meta:
    sc.pl.umap(adata, show=False, color=col)
    plt.savefig(f'{output_dir}/umap_scpoli_{col}.png', bbox_inches='tight')

#load geneformer embedding
geneformer_latent = np.load("/out/data/processed/geneformer_v2.npy")
adata.obsm["X_geneformer"] = geneformer_latent
sc.pp.neighbors(adata, use_rep="X_geneformer")
sc.tl.umap(adata, min_dist = 0.3)
for col in meta:
    sc.pl.umap(adata, show=False, color=col)
    plt.savefig(f'{output_dir}/umap_geneformer_{col}.png', bbox_inches='tight')

#load scGPT embedding
scgpt_latent = np.load("/out/data/processed/scGPT.npy")
adata.obsm["X_scGPT"] = scgpt_latent
sc.pp.neighbors(adata, use_rep="X_scGPT")
sc.tl.umap(adata, min_dist = 0.3)
for col in meta:
    sc.pl.umap(adata, show=False, color=col)
    plt.savefig(f'{output_dir}/umap_scGPT_{col}.png', bbox_inches='tight')


adata.write_h5ad("/out/data/HKOCA_concatenated_harmonized_processed_filtered_raw_annotated_snapseed-indiv_embeddings.h5ad")
