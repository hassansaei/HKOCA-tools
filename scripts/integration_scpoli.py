import os
os.environ["SCIPY_ARRAY_API"]="1"


import scanpy as sc
from snapseed.utils import read_yaml
import snapseed as snap
import anndata as ad
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scarches.models.scpoli import scPoli


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


scpoli_model = scPoli(
    adata=adata,
    condition_keys="sample_id",
    cell_type_keys="celltype",
    embedding_dims=5,
    hidden_layer_sizes=[1024]
)

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

scpoli_model.train(
    unlabeled_prototype_training=False,
    n_epochs=7,
    pretraining_epochs=5,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=10,
    alpha_epoch_anneal=100
)



output_dir = "/out/data/scpoli"
os.makedirs(output_dir, exist_ok=True)

scpoli_model.save(f"{output_dir}/scpoli_model.pkl")

scpoli_model.get_conditional_embeddings().write_h5ad(
    f"{output_dir}/hkoca-tool_conditional_embeddings.h5ad")

adata.X = adata.X.to_array()

scpoli_latent = scpoli_model.get_latent(
    adata,
    mean = True
)

out_ad = ad.AnnData(
    X=scpoli_latent
)

out_ad.obs.index = adata.obs.index

out_ad.write_h5ad(f"{output_dir}/scpoli_latent.h5ad")
