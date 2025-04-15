import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
sc._settings.ScanpyConfig.n_jobs = -1
import pandas as pd
import numpy as np
import os
import argparse
# input parameters from command line

def parse_args():
    """
    Parse command line arguments.
    
    Returns:
    - args: Parsed arguments
    """
    parser = argparse.ArgumentParser(description='Filter cells based on QC metrics.')
    parser.add_argument("h5ad", type=str,
                    help="path to h5ad file")
    parser.add_argument("params", type=str,
                    help="path to params file")
    args = parser.parse_args()
    return args


def calculate_stat(adata):
    """
    calucalte qc metrics if were not calculated before

    Parameters:
    - adata: AnnData object
    """
    if "total_counts" not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata)

def parse_filter_params(path, adata, dataset):
    """
    parses file with filter parameters based on the dataset
    Verifies that the parameters are in the right format and present for all datasets
    Parameters:
    - path: path to the filter parameters file
    - adata: AnnData object
    - dataset: dataset type (e.g., 'sample', 'study')
    - dataset_name: name of the dataset
    """
    df = pd.read_csv(path, sep=",", header=0, index_col=0)
    # check that the columns are in the right format
    columns = ["n_genes_by_count_min", "n_genes_by_count_max","total_counts_min","total_counts_max"]
    if not all(col in df.columns for col in columns):
        raise ValueError(f"The filter parameters file must contain the columns: {columns}")
    # check that the index contains all the datasets from adata
    if not all(dts in adata.obs[dataset].unique() for dts in df.index):
        raise ValueError("The filter parameters file must contain all the datasets from adata")
    df = df[columns]
    return df

    

def filter_cells(adata, dataset: str, dataset_name: str, min_genes=200, max_genes=2500, min_counts=500, max_counts=10000):
    """
    Filter cells based on the number of genes and counts.
    
    Parameters:
    - adata: AnnData object
    - min_genes: Minimum number of genes per cell
    - max_genes: Maximum number of genes per cell
    - min_counts: Minimum number of counts per cell
    - max_counts: Maximum number of counts per cell
    - dataset: Dataset type (e.g., 'sample', 'study')
    - dataset_name: Name of the dataset
    
    Returns:
    - filtered_adata: Filtered AnnData object
    """
    
    # Filter cells based on the specified criteria
    filtered_adata = adata[(adata.obs[dataset] == dataset_name) &
                           (adata.obs['n_genes_by_counts'] >= min_genes) & 
                           (adata.obs['n_genes_by_counts'] <= max_genes) & 
                           (adata.obs['total_counts'] >= min_counts) & 
                           (adata.obs['total_counts'] <= max_counts)]
    
    return filtered_adata

def plot_qc_metrics(adata, dataset: str, dataset_name: str, filtered = False, prefix = None):
    """
    Plot QC metrics for the dataset.
    
    Parameters:
    - adata: AnnData object
    - dataset: Dataset type (e.g., 'sample', 'study')
    - dataset_name: Name of the dataset
    - filtered: Boolean indicating if the data is filtered
    - prefix: Prefix for the output file names (optional)
    """
    # check if the dataset is in the adata object
    if dataset not in adata.obs.columns:
        raise ValueError(f"The dataset '{dataset}' is not present in the AnnData object.")
    # check that output directory exists
    if not os.path.exists('/out/qc/'):
        os.makedirs('/out/qc/')

    #select the dataset
    adata_sel = adata[adata.obs[dataset] == dataset_name]
    
    # Create a violin plot for the number of genes and counts
    
    plt.figure(figsize=(12, 6))
    sns.displot(adata_sel.obs["total_counts"], bins=100, kde=False)
    if filtered:
        plt.title(f'Filtered Read Count for {dataset_name}')
        plt.savefig(f'/out/qc/{dataset_name}_{prefix}_filtered_read_count.png')
    else:
        plt.title(f'Unfiltered Read Count for {dataset_name}')
        plt.savefig(f'/out/qc/{dataset_name}_{prefix}_unfiltered_read_count.png')
    
    plt.close()

    plt.figure(figsize=(12, 6))
    sns.displot(adata_sel.obs["n_genes_by_counts"], bins=100, kde=False)
    if filtered:
        plt.title(f'Filtered Gene Count for {dataset_name}')
        plt.savefig(f'/out/qc/{dataset_name}_{prefix}_filtered_gene_count.png')
    else:
        plt.title(f'Unfiltered Gene Count for {dataset_name}')
        plt.savefig(f'/out/qc/{dataset_name}_{prefix}_unfiltered_gene_count.png')
    plt.close()

    plt.figure(figsize=(12, 6))
    sc.pl.umap(adata_sel, color = "n_genes_by_counts")
    if filtered:
        plt.title(f'Filtered UMAP Gene Count for {dataset_name}')
        plt.savefig(f'/out/qc/{dataset_name}_{prefix}_filtered_umap_gene_count.png')
    else:
        plt.title(f'Unfiltered UMAP Gene Count for {dataset_name}')
        plt.savefig(f'/out/qc/{dataset_name}_{prefix}_unfiltered_umap_gene_count.png')
    plt.close()

#Function that will select in loop the datasets and filter them and outputs the filtered adata file which conatins all the datasets
def filter_adata(adata, dataset: str, filter_params: pd.DataFrame):
    """
    Filter cells based on the number of genes and counts for all datasets.
    
    Parameters:
    - adata: AnnData object
    - dataset: Dataset type (e.g.,'sample', 'study')
    - filter_params: DataFrame containing filter parameters for each dataset
    """
    filtered_adata_list = []
    for dataset_name in adata.obs[dataset].unique():
        print(f"Saving dataset: {dataset_name}")
        # Get filter parameters for the current dataset
        min_genes = filter_params.loc[dataset_name, 'n_genes_by_count_min']
        max_genes = filter_params.loc[dataset_name, 'n_genes_by_count_max']
        min_counts = filter_params.loc[dataset_name, 'total_counts_min']
        max_counts = filter_params.loc[dataset_name, 'total_counts_max']

        # Filter cells
        filtered_adata = filter_cells(adata, dataset, dataset_name, min_genes, max_genes, min_counts, max_counts)
        filtered_adata_list.append(filtered_adata)
    # Concatenate all filtered datasets
    adata_updated = sc.concat(filtered_adata_list, join="outer", label=dataset)   
    return adata_updated 
        
    


if __name__ == "__main__":
    args = parse_args()
    print("loading data")
    adata = sc.read_h5ad(args.h5ad)
    # Calculate QC metrics
    calculate_stat(adata)
    dataset = 'study'

    # Parse filter parameters
    filter_params = parse_filter_params(args.params, adata, 'study')

    for dataset_name in adata.obs[dataset].unique():
        print(f"Filtering dataset: {dataset_name}")
        # Get filter parameters for the current dataset
        min_genes = filter_params.loc[dataset_name, 'n_genes_by_count_min']
        max_genes = filter_params.loc[dataset_name, 'n_genes_by_count_max']
        min_counts = filter_params.loc[dataset_name, 'total_counts_min']
        max_counts = filter_params.loc[dataset_name, 'total_counts_max']

        # Filter cells
        filtered_adata = filter_cells(adata, dataset, dataset_name, min_genes, max_genes, min_counts, max_counts)
        unfiltered_adata = adata[adata.obs[dataset] == dataset_name]
        
        # Plot QC metrics
        plot_qc_metrics(filtered_adata, dataset, dataset_name, filtered=True, prefix='processed')
        plot_qc_metrics(unfiltered_adata, dataset, dataset_name, filtered=False, prefix='processed')
        # create a filtered adata object
    
    adata_filtered = filter_adata(adata, dataset, filter_params)

    # Save the filtered AnnData object
    print("Saving filtered data")
    adata_filtered.write('/mnt/data/processed/HKOCA_concatenated_harmonized_processed_filtered.h5ad')
