import scanpy as sc
from snapseed.utils import read_yaml
import snapseed as snap
import anndata as ad
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import os
os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/'

file_path = "/mnt/data/processed/HKOCA_concatenated_harmonized_processed_filtered_raw.h5ad"
adata = ad.read_h5ad(file_path)
adata.layers['raw'] = adata.X.copy()
adata.X = adata.layers['X_norm']

marker_genes = read_yaml("/out/data/snapseed_markers.yaml")


# Get unique sample IDs
sample_ids = adata.obs['sample_id'].unique()

obs_df = pd.DataFrame()

# Loop over each sample ID
for sample_id in sample_ids:
    # Subset the data for the current sample ID
    sample_adata = adata[adata.obs['sample_id'] == sample_id].copy()
    # Log1p transformation
    sc.pp.log1p(sample_adata)
    # Find highly variable genes
    sc.pp.highly_variable_genes(sample_adata, n_top_genes=2000)
    # Scale the data
    sc.pp.scale(sample_adata)
    # PCA dimension reduction
    sc.tl.pca(sample_adata, n_comps=50)
    # Compute the neighborhood graph
    sc.pp.neighbors(sample_adata)
    print(f'leiden clustering of {sample_id}')
    # Leiden clustering
    sc.tl.leiden(sample_adata)
    #Annotate using snapseed
    annot = snap.annotate_hierarchy(
	    sample_adata,
	    marker_genes,
	    group_name="leiden"
	)

    # annotate individual cells
    sample_adata.obs = sample_adata.obs.merge(annot['assignments'], how='left', left_on='leiden', right_index=True)
    # transfer annotation to the main object
    #adata.obs.loc[sample_adata.obs.index, ["leiden","level_1", "level_2", "level_3"]] = sample_adata.obs[["leiden","level_1", "level_2", "level_3"]]
    obs_df = pd.concat([obs_df, sample_adata.obs], axis=0)

    print(f'{sample_id} is processed')

# Save the annotated data
adata.obs = adata.obs.merge(obs_df[["level_1", "level_2", "level_3"]], how='left', left_index=True, right_index=True)
#save adata object
adata.write_h5ad("/mnt/data/processed/HKOCA_concatenated_harmonized_processed_filtered_raw_snapseed.h5ad")
