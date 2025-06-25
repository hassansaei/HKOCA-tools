import argparse
import scanpy as sc
import anndata as ad

import scvi
import torch
torch.set_float32_matmul_precision("medium")
import warnings
warnings.filterwarnings('ignore')
import os
import pandas as pd


def combine_annotations(row):
    if pd.isna(row['level_3']):
        return row['level_2']
    else:
        return row['level_3']

output_dir = "/out/data/scvi"
os.makedirs(output_dir, exist_ok=True)


file_path = "/out/data/HKOCA_concatenated_harmonized_processed_filtered_raw_annotated_snapseed-indiv.h5ad"
adata = ad.read_h5ad(file_path)
adata.obs['celltype'] = adata.obs.apply(combine_annotations, axis=1)
adata.X = adata.layers["raw"].copy()


scvi.model.SCVI.setup_anndata(adata, batch_key="sample_id")

scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
scvi_model.train()
scvi_model.save(f"{output_dir}/scvi_model.pt")

scvi_latent = scvi_model.get_latent_representation()

out_ad = ad.AnnData(
    X=scvi_latent
)

out_ad.obs.index = adata.obs.index

out_ad.write_h5ad(f"{output_dir}/scvi_latent.h5ad")

