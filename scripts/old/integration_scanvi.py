import argparse
import os
import pandas as pd

import scanpy as sc
import anndata as ad

import scvi
import torch
torch.set_float32_matmul_precision("medium")
import warnings
warnings.filterwarnings('ignore')    

output_dir = "/out/data/scvi"

def combine_annotations(row):
    if pd.isna(row['level_3']):
        return row['level_2']
    else:
        return row['level_3']

file_path = "/out/data/HKOCA_concatenated_harmonized_processed_filtered_raw_annotated_snapseed-indiv.h5ad"
adata = ad.read_h5ad(file_path)
adata.obs['celltype'] = adata.obs.apply(combine_annotations, axis=1)
adata.obs['celltype'] = adata.obs['celltype'].astype(str)
adata.X = adata.layers["raw"].copy()

vae = scvi.model.SCVI.load(f"{output_dir}/scvi_model.pt", adata=adata)
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="celltype",
    unlabeled_category="unknown",
)
lvae.train()

scanvi_latent = lvae.get_latent_representation(adata)

out_ad = ad.AnnData(
    X=scanvi_latent
)
out_ad.obs.index = adata.obs.index

out_ad.write_h5ad(f"{output_dir}/scanvi_latent.h5ad")
