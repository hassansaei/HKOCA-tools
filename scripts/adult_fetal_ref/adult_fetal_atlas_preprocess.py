import argparse
import logging
import os
import sys
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def setup_logging(log_file_path=None):
    """Initializes logging configuration."""
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    handlers = [logging.StreamHandler(sys.stdout)]
    
    if log_file_path:
        handlers.append(logging.FileHandler(log_file_path))
        
    logging.basicConfig(level=logging.INFO, format=log_format, handlers=handlers)

def load_config(config_path):
    """Loads YAML configuration file."""
    try:
        with open(config_path, 'r') as file:
            return yaml.safe_load(file)
    except Exception as e:
        logging.error(f"Failed to load configuration file at {config_path}: {e}")
        sys.exit(1)

def setup_directories(config):
    """Creates all required output directories based on the configuration."""
    base_dir = Path(config["paths"]["base_dir"])
    dirs = {
        "out_ref": base_dir / config["paths"]["out_ref"],
        "out_adult": base_dir / config["paths"]["out_adult"],
        "out_fetal": base_dir / config["paths"]["out_fetal"],
        "out_req": base_dir / config["paths"]["out_req"]
    }
    
    for name, dir_path in dirs.items():
        try:
            dir_path.mkdir(parents=True, exist_ok=True)
            logging.info(f"Directory ready: {dir_path}")
        except Exception as e:
            logging.error(f"Failed to create directory {dir_path}: {e}")
            sys.exit(1)
            
    return dirs, base_dir

def checkpoint_overview(adata, label):
    """Logs dataset overview and cell type distributions."""
    cts = adata.obs["cell_type"].value_counts()
    logging.info("-" * 55)
    logging.info(f"  {label}")
    logging.info("-" * 55)
    logging.info(f"  Total cells   : {adata.n_obs:,}")
    logging.info(f"  Total genes   : {adata.n_vars:,}")
    logging.info(f"  Cell types    : {cts.shape[0]}")
    logging.info(f"\n  {'Cell type':<40} {'Count':>6}  {'%':>6}")
    logging.info("-" * 56)
    
    for ct, n in cts.items():
        logging.info(f"  {ct:<40} {n:>6,}  {100*n/adata.n_obs:>5.1f}%")
    logging.info("-" * 55)
    return cts

def apply_cell_type_consolidation(adata, mapping_dict):
    """Remaps cell types according to the provided mapping dictionary."""
    try:
        adata.obs["cell_type"] = adata.obs["cell_type"].replace(mapping_dict)
        logging.info("Cell type consolidation applied successfully.")
        return adata
    except KeyError as e:
        logging.error(f"Failed to apply cell type consolidation. Missing column: {e}")
        raise

def extract_markers(adata, label, groupby="cell_type", n_genes=10000):
    """Runs Wilcoxon rank-genes using raw counts for percentage computation."""
    try:
        if adata.raw is None:
            adata.raw = adata.copy()
        
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method="wilcoxon",
            n_genes=n_genes,
            pts=True,              
            use_raw=True,          
        )
        result = sc.get.rank_genes_groups_df(adata, group=None)
        result.columns = ["group", "gene", "logFC", "score", "pval", "pval_adj", "pct_in", "pct_out"]
        
        logging.info(f"Markers extracted - {label}")
        logging.info(f"  Groups : {result['group'].nunique()}")
        logging.info(f"  Rows   : {len(result):,}")
        logging.info(f"  pct_in  - min: {result['pct_in'].min():.3f}  max: {result['pct_in'].max():.3f}  mean: {result['pct_in'].mean():.3f}")
        logging.info(f"  pct_out - min: {result['pct_out'].min():.3f}  max: {result['pct_out'].max():.3f}  mean: {result['pct_out'].mean():.3f}")
        return result
    except Exception as e:
        logging.error(f"Marker extraction failed for {label}: {e}")
        raise

def check_threshold_elbow(marker_df, cell_type, tolerance, params):
    """Diagnostic tool to test stringent score and elbow method limits."""
    try:
        df = marker_df[marker_df['group'] == cell_type].copy()
        df = df[~df['gene'].str.startswith(params['exclude_gene_prefix'])]
        df = df[(df['pct_in'] >= params['min_pct_in']) & 
                (df['logFC'] >= params['min_logfc']) & 
                (df['pval_adj'] < params['max_pval_adj'])]
        
        if df.empty:
            logging.warning(f"No genes passed baseline filters for diagnostic test on '{cell_type}'.")
            return None
        
        df['stringent_score'] = df['logFC'] * (df['pct_in'] ** 2) / (df['pct_out'] + tolerance)
        df_sorted = df.sort_values(by='stringent_score', ascending=False).head(params['max_n']).copy()
        
        y = df_sorted['stringent_score'].values
        x = np.arange(len(y))
        
        if len(y) > 2:
            x_norm = x / x.max()
            y_norm = (y - y.min()) / (y.max() - y.min() + 1e-9)
            distances = 1.0 - (x_norm + y_norm)
            elbow_idx = np.argmax(distances) 
        else:
            elbow_idx = len(y) - 1 

        cutoff_rank = max(elbow_idx + 1, params['min_kept'])
        
        df_sorted['rank'] = x + 1
        df_sorted['Action'] = np.where(df_sorted['rank'] <= cutoff_rank, 'Kept', 'Dropped')
        df_sorted['stringent_score'] = df_sorted['stringent_score'].round(2)
        df_sorted['pct_in'] = df_sorted['pct_in'].round(3)
        df_sorted['pct_out'] = df_sorted['pct_out'].round(3)
        df_sorted['logFC'] = df_sorted['logFC'].round(2)
        
        logging.info(f"--- Diagnostic: {cell_type} ---")
        logging.info(f"Tolerance Dial     : {tolerance}")
        logging.info(f"Elbow detected at rank : {elbow_idx + 1}")
        logging.info(f"Final Cutoff Applied : Rank {cutoff_rank}")
        
        cols = ['rank', 'Action', 'gene', 'stringent_score', 'pct_in', 'pct_out', 'logFC']
        return df_sorted[cols].reset_index(drop=True)
    except Exception as e:
        logging.error(f"Diagnostic test failed for {cell_type}: {e}")
        return None

def normalize_name(name):
    """Normalizes labels for consistency."""
    return str(name).strip().replace(" ", "_").replace("/", "_")

def generate_canonical_maps(raw_df, raw_map):
    """Dynamically builds final canonical mappings based on unique group values."""
    final_map = {}
    for raw in raw_df['group'].unique():
        norm = normalize_name(raw)
        final_map[raw] = raw_map.get(norm, norm)
    return final_map

def export_final_markers_to_tsv(marker_df, output_folder, tolerance, dataset_label, canonical_map, params):
    """Filters, calculates stringent scores, applies elbow method, and exports TSVs."""
    try:
        df_all = marker_df.copy()
        df_all["canonical"] = df_all["group"].map(canonical_map).fillna(dataset_label + "_" + df_all["group"].astype(str))
        canonical_types = df_all["canonical"].unique()
        logging.info(f"Processing {len(canonical_types)} canonical cell types for {dataset_label} into '{output_folder}/'")

        for canonical_name in sorted(canonical_types):
            df = df_all[df_all["canonical"] == canonical_name].copy()
            df = df[~df["gene"].str.startswith(params['exclude_gene_prefix'])]
            df = df[(df["pct_in"] >= params['min_pct_in']) & 
                    (df["logFC"] >= params['min_logfc']) & 
                    (df["pval_adj"] < params['max_pval_adj'])]

            if df.empty:
                continue

            df["stringent_score"] = df["logFC"] * (df["pct_in"] ** 2) / (df["pct_out"] + tolerance)
            df = df.sort_values("stringent_score", ascending=False).drop_duplicates(subset="gene", keep="first")
            df_sorted = df.head(params['max_n']).copy()

            y = df_sorted["stringent_score"].values
            x = np.arange(len(y))

            if len(y) > 2:
                x_norm = x / x.max()
                y_norm = (y - y.min()) / (y.max() - y.min() + 1e-9)
                distances = 1.0 - (x_norm + y_norm)
                elbow_idx = np.argmax(distances)
            else:
                elbow_idx = len(y) - 1

            cutoff_rank = max(elbow_idx + 1, params['min_kept'])
            final_genes = df_sorted.iloc[:cutoff_rank].copy()

            final_genes["stringent_score"] = final_genes["stringent_score"].round(2)
            final_genes["pct_in"]          = final_genes["pct_in"].round(3)
            final_genes["pct_out"]         = final_genes["pct_out"].round(3)
            final_genes["logFC"]           = final_genes["logFC"].round(2)
            final_genes["group"]           = canonical_name   

            cols = ["group", "gene", "stringent_score", "pct_in", "pct_out", "logFC", "pval_adj"]
            cols += [c for c in final_genes.columns if c not in cols]
            final_genes = final_genes[cols]

            safe_name = canonical_name.replace("/", "_").replace(" ", "_")
            file_path = output_folder / f"{dataset_label}_{safe_name}.tsv"
            final_genes.to_csv(file_path, sep="\t", index=False)
            logging.info(f"  Saved {len(final_genes)} genes to {file_path.name}")
    except Exception as e:
        logging.error(f"Failed during TSV export for {dataset_label}: {e}")
        raise

def extract_top_genes_for_deliverables(df_raw, canonical_map, tolerance, params, is_adult=False):
    """Extracts required metrics for deliverables."""
    df = df_raw.copy()
    df["canonical"] = df["group"].map(canonical_map)
    
    df = df[~df["gene"].str.startswith(params['exclude_gene_prefix'])]
    df = df[(df["pct_in"] >= params['min_pct_in']) & 
            (df["logFC"] >= params['min_logfc']) & 
            (df["pval_adj"] < params['max_pval_adj'])]
    df["stringent_score"] = df["logFC"] * (df["pct_in"] ** 2) / (df["pct_out"] + tolerance)
    
    top_overlap_dict, top_export_list = {}, []
    
    for canon in df["canonical"].unique():
        sub = df[df["canonical"] == canon].copy()
        if is_adult and canon == "endothelial_cell":
            sub["norm_score"] = sub.groupby("group")["stringent_score"].transform(lambda x: x / (x.max() + 1e-9))
            agg = sub.groupby("gene").agg(freq=("group", "count"), mean_norm=("norm_score", "mean"))
            agg["score"] = agg["freq"] * agg["mean_norm"]
            ranked = agg.sort_values(["freq", "score"], ascending=[False, False]).reset_index()
        else:
            ranked = sub.groupby("gene")["stringent_score"].max().sort_values(ascending=False).reset_index()
            ranked.rename(columns={"stringent_score": "score"}, inplace=True)
            
        ranked["canonical_cell_type"] = canon
        top_overlap_dict[canon] = set(ranked.head(params['overlap_n'])["gene"]) 
        
        top_export = ranked.head(params['export_n']).copy()
        top_export["rank"] = range(1, len(top_export) + 1)
        top_export_list.append(top_export[["canonical_cell_type", "gene", "score", "rank"]])

    return pd.concat(top_export_list, ignore_index=True), top_overlap_dict, df

def qc_filter_low_abundance(adata, atlas_name, min_cells):
    """Removes cell types with insufficient cells to prevent spurious marker calls."""
    if "cell_type" not in adata.obs:
        return adata
        
    counts = adata.obs["cell_type"].value_counts()
    valid_types = counts[counts >= min_cells].index
    dropped = counts[counts < min_cells]
    
    if not dropped.empty:
        logging.warning(f"[{atlas_name}] QC ALERT: Dropping cell types with < {min_cells} cells: {dropped.to_dict()}")
        
    adata_filtered = adata[adata.obs["cell_type"].isin(valid_types)].copy()
    adata_filtered.obs["cell_type"] = adata_filtered.obs["cell_type"].astype("category")
    adata_filtered.obs["cell_type"] = adata_filtered.obs["cell_type"].cat.remove_unused_categories()
    return adata_filtered


def qc_post_mapping_checks(adult_canons, fetal_canons, corr_matrix, expected_lineages, min_corr):
    """Validates the presence of expected major lineages and checks cross-atlas overlap correlation."""
    # 1. Expected Major Lineages Check
    missing_adult = [l for l in expected_lineages if l not in adult_canons]
    missing_fetal = [l for l in expected_lineages if l not in fetal_canons]
    
    if missing_adult:
        logging.warning(f"[QC FAIL] Adult atlas is missing expected core lineages: {missing_adult}")
    if missing_fetal:
        logging.warning(f"[QC FAIL] Fetal atlas is missing expected core lineages: {missing_fetal}")

    # 2. Minimum Overlap Score Threshold Check
    shared_types = set(adult_canons).intersection(set(fetal_canons))
    
    for ct in shared_types:
        if ct in corr_matrix.index and ct in corr_matrix.columns:
            direct_correlation = corr_matrix.loc[ct, ct]
            
            if pd.isna(direct_correlation) or direct_correlation < min_corr:
                logging.warning(f"[QC WARNING] Low cross-atlas overlap for '{ct}'. "
                                f"Pearson r = {direct_correlation:.2f}. Potential mapping failure.")


def merge_pan_endothelial(out_adult_dir, endo_files, top_n, merged_group_name, merged_file_name):
    """Merges specified endothelial subsets into a unified pan-endothelial marker list."""
    try:
        dfs = []
        for file_name in endo_files:
            file_path = out_adult_dir / file_name
            if file_path.exists():
                df = pd.read_csv(file_path, sep="\t")
                if not df.empty:
                    max_score = df["stringent_score"].max()
                    df["norm_score"] = df["stringent_score"] / max_score
                    dfs.append(df)
            else:
                logging.warning(f"Could not find {file_name} for pan-endothelial merge.")

        if dfs:
            logging.info(f"Merging {len(dfs)} endothelial subsets into unified marker list.")
            merged_df = pd.concat(dfs, ignore_index=True)

            agg_df = merged_df.groupby("gene").agg(
                freq=("group", "count"),                      
                mean_norm_score=("norm_score", "mean"),       
                pct_in=("pct_in", "mean"),                    
                pct_out=("pct_out", "mean"),                  
                logFC=("logFC", "mean"),                      
                pval_adj=("pval_adj", "min")                  
            ).reset_index()

            agg_df["pan_endothelial_score"] = agg_df["freq"] * agg_df["mean_norm_score"]
            top_endo = agg_df.sort_values(["freq", "pan_endothelial_score"], ascending=[False, False]).head(top_n).copy()

            top_endo["group"] = merged_group_name
            top_endo["stringent_score"] = top_endo["pan_endothelial_score"].round(2)
            top_endo["pct_in"] = top_endo["pct_in"].round(3)
            top_endo["pct_out"] = top_endo["pct_out"].round(3)
            top_endo["logFC"] = top_endo["logFC"].round(2)

            final_cols = ["group", "gene", "stringent_score", "pct_in", "pct_out", "logFC", "pval_adj"]
            top_endo = top_endo[final_cols]

            output_path = out_adult_dir / merged_file_name
            top_endo.to_csv(output_path, sep="\t", index=False)
            logging.info(f"Generated '{merged_file_name}' with {len(top_endo)} markers.")
        else:
            logging.warning("No endothelial TSV files found to merge.")
    except Exception as e:
        logging.error(f"Failed during pan-endothelial merge: {e}")
        raise
                
def generate_deliverables(adult_adata, fetal_adata, adult_markers_raw, fetal_markers_raw, adult_canonical_map, fetal_canonical_map, dirs, config):
    """Compiles and exports the comprehensive reporting outputs required."""
    try:
        out_req = dirs["out_req"]
        params = config["deliverables"]
        logging.info(f"Processing outputs to: {out_req}/")

        # 1. Label Mapping
        mapping_data = [{"Original_Label": k, "Canonical_Label": v, "Atlas": "Adult"} for k, v in adult_canonical_map.items()]
        mapping_data += [{"Original_Label": k, "Canonical_Label": v, "Atlas": "Fetal"} for k, v in fetal_canonical_map.items()]
        pd.DataFrame(mapping_data).to_csv(out_req / "cell_type_label_mapping.tsv", sep="\t", index=False)

        # 2 & 3. Top N TSVs
        adult_top_df, adult_top_genes, adult_filtered = extract_top_genes_for_deliverables(
            adult_markers_raw, adult_canonical_map, config["elbow_params"]["adult_tolerance"], params, is_adult=True)
        fetal_top_df, fetal_top_genes, fetal_filtered = extract_top_genes_for_deliverables(
            fetal_markers_raw, fetal_canonical_map, config["elbow_params"]["fetal_tolerance"], params, is_adult=False)
        
        adult_top_df.to_csv(out_req / "adult_kidney_markers_top20.tsv", sep="\t", index=False)
        fetal_top_df.to_csv(out_req / "fetal_kidney_markers_top20.tsv", sep="\t", index=False)

        # 4. Cell Type Comparison
        adult_canons = set(adult_canonical_map.values())
        fetal_canons = set(fetal_canonical_map.values())
        all_canons = sorted(list(adult_canons | fetal_canons))

        comparison_data = [
            {"Canonical_Cell_Type": ct, "In_Adult_Atlas": "Yes" if ct in adult_canons else "No", "In_Fetal_Atlas": "Yes" if ct in fetal_canons else "No"}
            for ct in all_canons
        ]
        pd.DataFrame(comparison_data).to_csv(out_req / "cell_type_comparison.tsv", sep="\t", index=False)

        # 5. Overlap Heatmap Matrices & Shared Output TSV
        overlap_records = []
        count_matrix = pd.DataFrame(index=sorted(list(adult_canons)), columns=sorted(list(fetal_canons)), dtype=int)

        for a_canon in sorted(list(adult_canons)):
            a_genes = adult_top_genes.get(a_canon, set())
            for f_canon in sorted(list(fetal_canons)):
                f_genes = fetal_top_genes.get(f_canon, set())
                
                intersect = a_genes.intersection(f_genes)
                count_matrix.loc[a_canon, f_canon] = int(len(intersect))
                
                if a_canon == f_canon:
                    for g in intersect:
                        a_score = adult_filtered[(adult_filtered["canonical"] == a_canon) & (adult_filtered["gene"] == g)]["stringent_score"].max()
                        f_score = fetal_filtered[(fetal_filtered["canonical"] == f_canon) & (fetal_filtered["gene"] == g)]["stringent_score"].max()
                        overlap_records.append({
                            "Canonical_Cell_Type": a_canon, "Gene": g,
                            "Adult_Stringent_Score": round(a_score, 2), "Fetal_Stringent_Score": round(f_score, 2)
                        })

        overlap_df = pd.DataFrame(overlap_records)
        if not overlap_df.empty:
            overlap_df["Combined_Score"] = overlap_df["Adult_Stringent_Score"] + overlap_df["Fetal_Stringent_Score"]
            overlap_df = overlap_df.sort_values(by=["Canonical_Cell_Type", "Combined_Score"], ascending=[True, False]).drop(columns=["Combined_Score"])
        overlap_df.to_csv(out_req / "marker_overlap_adult_vs_fetal.tsv", sep="\t", index=False)

        # 6. Marker Overlap Figure (Pseudobulk Correlation + Counts)
        if config["plotting"]["save_figs"]:
            logging.info("Generating pseudobulk correlation figure...")
            
            # Add canonical labels directly to the AnnData copies locally
            adult_adata.obs["canonical_label"] = adult_adata.obs["cell_type"].map(adult_canonical_map)
            fetal_adata.obs["canonical_label"] = fetal_adata.obs["cell_type"].map(fetal_canonical_map)
            
            a_sub = adult_adata[adult_adata.obs["canonical_label"].notna()]
            f_sub = fetal_adata[fetal_adata.obs["canonical_label"].notna()]
            
            # Combine all top marker genes
            adult_genes_all = set([g for subset in adult_top_genes.values() for g in subset])
            fetal_genes_all = set([g for subset in fetal_top_genes.values() for g in subset])
            all_genes = list(adult_genes_all | fetal_genes_all)
            
            # Only keep genes that exist in both datasets
            shared_genes = [g for g in all_genes if g in a_sub.var_names and g in f_sub.var_names]
            
            # Create pseudobulk profiles
            a_df = sc.get.obs_df(a_sub, keys=shared_genes + ["canonical_label"])
            f_df = sc.get.obs_df(f_sub, keys=shared_genes + ["canonical_label"])
            
            a_pb = a_df.groupby("canonical_label").mean()
            f_pb = f_df.groupby("canonical_label").mean()
            
            # Align the count matrix (from step 5) with the pseudobulk index
            common_adult = [t for t in a_pb.index if t in count_matrix.index]
            common_fetal = [t for t in f_pb.index if t in count_matrix.columns]
            
            a_pb = a_pb.loc[common_adult]
            f_pb = f_pb.loc[common_fetal]
            c_mat = count_matrix.loc[common_adult, common_fetal].fillna(0).astype(int)
            
            # Compute Pearson correlation matrix
            corr_matrix = pd.DataFrame(index=a_pb.index, columns=f_pb.index, dtype=float)
            for a_type in a_pb.index:
                for f_type in f_pb.index:
                    corr_matrix.loc[a_type, f_type] = a_pb.loc[a_type].corr(f_pb.loc[f_type], method="pearson")
            
            # Change this line inside generate_deliverables function
            qc_post_mapping_checks(
                adult_canons, fetal_canons, corr_matrix, 
                expected_lineages=config["qc"]["expected_lineages"], 
                min_corr=config["qc"]["min_corr"]
            )

            # Plot
            plt.figure(figsize=(14, 12))
            sns.heatmap(corr_matrix, cmap="RdBu_r", center=0, vmin=-1, vmax=1, 
                        annot=c_mat, fmt="d", annot_kws={"size": 10}, 
                        cbar_kws={'label': 'Pearson Correlation (Color)'})
            
            plt.title(f"Cross-Atlas Similarity Across Top {params['overlap_n']} Genes\n(Colors: Transcriptomic Correlation | Text: Shared Marker Genes)", 
                      fontsize=16, pad=20)
            plt.xlabel("Fetal Cell Types (Canonical)", fontsize=14)
            plt.ylabel("Adult Cell Types (Canonical)", fontsize=14)
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout()
            plt.savefig(out_req / "marker_overlap_figure.png", dpi=300, bbox_inches="tight")
            plt.close()

        # 7. Marker Method Params YAML
        export_params = {
            "baseline_filters": {
                "min_pct_in": params["min_pct_in"], 
                "min_logFC": params["min_logfc"], 
                "max_pval_adj": params["max_pval_adj"], 
                "exclude_gene_prefix": params["exclude_gene_prefix"]
            },
            "stringent_score_formula": "logFC * (pct_in ^ 2) / (pct_out + tolerance)",
            "tolerances": {"adult": config["elbow_params"]["adult_tolerance"], "fetal": config["elbow_params"]["fetal_tolerance"]},
            "exports": {"top_n_individual_tsv": params["export_n"], "top_n_overlap_heatmap": params["overlap_n"]}
        }
        with open(out_req / "marker_method_params.yaml", "w") as f:
            yaml.dump(export_params, f, sort_keys=False)
            
        logging.info("Deliverables exported successfully.")
    except Exception as e:
        logging.error(f"Failed to generate deliverables: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Reference Atlas Gene Marker Extraction Pipeline")
    parser.add_argument("--config", type=str, required=True, help="Path to the configuration YAML file.")
    args = parser.parse_args()

    setup_logging()
    logging.info("Starting Marker Extraction Pipeline...")
    
    config = load_config(args.config)
    dirs, base_dir = setup_directories(config)
    
    # Configure Scanpy/Matplotlib
    sc.settings.verbosity = config["scanpy_settings"]["verbosity"]
    sc.settings.set_figure_params(
        dpi=config["scanpy_settings"]["dpi"], 
        facecolor=config["scanpy_settings"]["facecolor"], 
        frameon=config["scanpy_settings"]["frameon"]
    )
    matplotlib.rcParams["figure.dpi"] = config["scanpy_settings"]["dpi"]

    # Load datasets
    try:
        adult_h5ad = base_dir / config["paths"]["adult_h5ad"]
        fetal_h5ad = base_dir / config["paths"]["fetal_h5ad"]
        logging.info(f"Loading adult dataset from {adult_h5ad}")
        adult = sc.read_h5ad(adult_h5ad)
        logging.info(f"Loading fetal dataset from {fetal_h5ad}")
        fetal = sc.read_h5ad(fetal_h5ad)
    except Exception as e:
        logging.error(f"Failed to load input h5ad files: {e}")
        sys.exit(1)

    # Standardize column mappings
    adult.obs["cell_type"] = adult.obs[config["columns"]["adult_celltype_col"]].astype(str)
    fetal.obs["cell_type"] = fetal.obs[config["columns"]["fetal_celltype_col"]].astype(str)

    # Checkpoints
    checkpoint_overview(adult, "ADULT REFERENCE")
    checkpoint_overview(fetal, "FETAL REFERENCE")

    # Cell type consolidation
    adult = apply_cell_type_consolidation(adult, config["mappings"]["endothelial_consolidation"])
    fetal = apply_cell_type_consolidation(fetal, config["mappings"]["endothelial_consolidation"])
    
    # Run Pre-Marker QC
    adult = qc_filter_low_abundance(adult, atlas_name="Adult", min_cells=config["qc"]["min_cells"])
    fetal = qc_filter_low_abundance(fetal, atlas_name="Fetal", min_cells=config["qc"]["min_cells"])

    # Marker Extraction
    adult_markers_raw = extract_markers(adult, "Adult reference", n_genes=config["wilcoxon"]["n_genes"])
    fetal_markers_raw = extract_markers(fetal, "Fetal reference", n_genes=config["wilcoxon"]["n_genes"])

    # Optional Diagnostic Tests
    if config.get("run_diagnostics", False):
        check_threshold_elbow(
            adult_markers_raw, 
            cell_type=config["diagnostics"]["adult_test_cell_type"], 
            tolerance=config["elbow_params"]["adult_tolerance"], 
            params=config["deliverables"]
        )
        check_threshold_elbow(
            fetal_markers_raw, 
            cell_type=config["diagnostics"]["fetal_test_cell_type"], 
            tolerance=config["elbow_params"]["fetal_tolerance"], 
            params=config["deliverables"]
        )

    # Canonical Mappings
    adult_canonical_map = generate_canonical_maps(adult_markers_raw, config["mappings"]["adult_canonical_raw"])
    fetal_canonical_map = generate_canonical_maps(fetal_markers_raw, config["mappings"]["fetal_canonical_raw"])

    # Pipeline TSV Exports
    export_final_markers_to_tsv(
        adult_markers_raw, dirs["out_adult"], config["elbow_params"]["adult_tolerance"], 
        "adult", adult_canonical_map, config["deliverables"]
    )
    export_final_markers_to_tsv(
        fetal_markers_raw, dirs["out_fetal"], config["elbow_params"]["fetal_tolerance"], 
        "fetal", fetal_canonical_map, config["deliverables"]
    )

    # Post-Processing: Pan-Endothelial Merge
    merge_pan_endothelial(
        dirs["out_adult"], 
        config["pan_endothelial"]["files_to_merge"], 
        config["pan_endothelial"]["top_n"],
        config["pan_endothelial"]["merged_group_name"],
        config["pan_endothelial"]["merged_file_name"]
    )

    generate_deliverables(
        adult, fetal,
        adult_markers_raw, fetal_markers_raw, 
        adult_canonical_map, fetal_canonical_map, 
        dirs, config
    )

    logging.info("Pipeline execution completed successfully.")

if __name__ == "__main__":
    main()