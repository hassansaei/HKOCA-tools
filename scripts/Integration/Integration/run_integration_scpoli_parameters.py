# scPoli full-atlas integration
# Configuration: scpoli_integration_config.yaml

import os
import gc
import sys
import logging
import warnings
import subprocess
import types
from pathlib import Path

import numpy as np
import scipy.sparse as sp

import anndata as ad

_boot_logger = logging.getLogger(__name__)

if not hasattr(ad, 'read'):
    ad.read = ad.read_h5ad

if not hasattr(ad, 'io'):
    io_module = types.ModuleType('anndata.io')
    _names = ['read_csv', 'read_loom', 'read_text', 'read_excel',
              'read_hdf', 'read_mtx', 'read_zarr', 'read_h5ad',
              'read_umi_tools', 'read_elem', 'write_elem']

    _experimental = None
    try:
        import anndata.experimental as _experimental
    except ImportError:
        pass

    for _name in _names:
        if hasattr(ad, _name):
            setattr(io_module, _name, getattr(ad, _name))
        elif _experimental is not None and hasattr(_experimental, _name):
            setattr(io_module, _name, getattr(_experimental, _name))
        else:
            _boot_logger.warning(f"anndata.{_name} not found anywhere; may cause issues downstream.")

    ad.io = io_module
    sys.modules['anndata.io'] = io_module

if not hasattr(ad, 'abc'):
    from anndata._core.sparse_dataset import CSRDataset, CSCDataset
    abc_module = types.ModuleType('anndata.abc')
    abc_module.CSRDataset = CSRDataset
    abc_module.CSCDataset = CSCDataset
    ad.abc = abc_module
    sys.modules['anndata.abc'] = abc_module

sys.modules['anndata'] = ad

import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from scarches.models import scPoli

try:
    import yaml
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyyaml", "-q"])
    import yaml

warnings.filterwarnings('ignore')


def load_config():
    argv = sys.argv[1:]
    if "--config" in argv:
        config_file = Path(argv[argv.index("--config") + 1]).resolve()
    else:
        config_file = Path(__file__).resolve().parent / "scpoli_integration_config.yaml"
    with open(config_file, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


CFG = load_config()

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

OUTPUT_DIR = CFG["paths"]["output_dir"]
os.makedirs(OUTPUT_DIR, exist_ok=True)

FINAL_H5AD_OUT = os.path.join(OUTPUT_DIR, CFG["paths"]["output_h5ad"])
BASELINE_H5AD = CFG["paths"]["baseline_h5ad"]

FEAT_CFG = CFG["feature_selection"]
PCA_CFG = CFG["pca"]
SCPOLI_CFG = CFG["scpoli"]
TRAIN_CFG = CFG["training"]
PLOT_CFG = CFG["plotting"]

sc.settings.figdir = OUTPUT_DIR
mpl.rcParams['figure.figsize'] = tuple(PLOT_CFG["figsize"])
mpl.rcParams['figure.dpi'] = PLOT_CFG["dpi"]
mpl.rcParams['savefig.dpi'] = PLOT_CFG["dpi"]
mpl.rcParams['savefig.format'] = 'png'
sc.set_figure_params(dpi=PLOT_CFG["dpi"], facecolor='white', frameon=False)

logger.info("Loading dataset...")
adata = sc.read_h5ad(BASELINE_H5AD)
logger.info(f"Initial data loaded: {adata.n_obs} cells.")

logger.info("Dropping unused layers and raw slot...")
if adata.raw is not None:
    del adata.raw
for layer in list(adata.layers.keys()):
    del adata.layers[layer]
gc.collect()

if not sp.issparse(adata.X):
    logger.info("Compressing dense matrix into CSR sparse format...")
    adata.X = sp.csr_matrix(adata.X)
gc.collect()

target_keys = CFG["metadata"]["sanitize_keys"]
for key in target_keys:
    if key in adata.obs.columns:
        adata.obs[key] = adata.obs[key].astype(str).replace(['nan', 'NaN', 'None', ''], np.nan)
        nan_count = adata.obs[key].isna().sum()
        if nan_count > 0:
            logger.warning(f"Found {nan_count} missing values in '{key}'. Purging...")
            valid_idx = adata.obs[key].notna()
            adata = adata[valid_idx]
            gc.collect()

if 'sample_id' in adata.obs.columns: adata.obs['sample_id'] = adata.obs['sample_id'].astype('category')
if 'study' in adata.obs.columns: adata.obs['study'] = adata.obs['study'].astype('category')
if 'Level_2' in adata.obs.columns: adata.obs['Level_2'] = adata.obs['Level_2'].astype('category')
if 'Level_1' in adata.obs.columns: adata.obs['Level_1'] = adata.obs['Level_1'].astype('category')

logger.info(f"Sanitization complete. Cleaned matrix has {adata.n_obs} cells remaining.")

logger.info(f"Selecting highly variable genes ({FEAT_CFG['n_top_genes']})...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=FEAT_CFG["n_top_genes"],
    flavor=FEAT_CFG["hvg_flavor"],
    batch_key=FEAT_CFG["batch_key"]
)
adata = adata[:, adata.var['highly_variable']].copy()
gc.collect()

logger.info("Processing baseline representation (PCA)...")
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=PCA_CFG["normalize_target_sum"])
sc.pp.log1p(adata)

if 'X_pca' not in adata.obsm:
    sc.pp.pca(adata, n_comps=PCA_CFG["n_comps"])

sc.pp.neighbors(adata, use_rep='X_pca', key_added='neighbors_pca', n_neighbors=PCA_CFG["n_neighbors"], random_state=PCA_CFG["random_state"])
sc.tl.umap(adata, neighbors_key='neighbors_pca', random_state=PCA_CFG["random_state"])
adata.obsm['X_umap_pca'] = adata.obsm['X_umap'].copy()

logger.info("Initializing scPoli framework targeting Level_2 biology...")

adata.X = adata.layers['counts'].copy()
del adata.layers['counts']
gc.collect()

model = scPoli(
    adata=adata,
    condition_keys=SCPOLI_CFG["condition_keys"],
    cell_type_keys=SCPOLI_CFG["cell_type_key"],
    hidden_layer_sizes=SCPOLI_CFG["hidden_layer_sizes"],
    latent_dim=SCPOLI_CFG["latent_dim"],
    embedding_dims=SCPOLI_CFG["embedding_dims"],
    recon_loss=SCPOLI_CFG["recon_loss"]
)

logger.info("Training scPoli model...")
model.train(
    n_epochs=TRAIN_CFG["n_epochs"],
    pretraining_epochs=TRAIN_CFG["pretraining_epochs"],
    eta=TRAIN_CFG["eta"],
    alpha_epoch_anneal=TRAIN_CFG["alpha_epoch_anneal"],
    batch_size=TRAIN_CFG["batch_size"],
    early_stopping_kwargs={
    "early_stopping_metric": TRAIN_CFG["early_stopping_metric"],
    "patience": TRAIN_CFG["patience"],
    "threshold": TRAIN_CFG["threshold"],
    "reduce_lr": TRAIN_CFG["reduce_lr"],
    "lr_patience": TRAIN_CFG["lr_patience"],
    "lr_factor": TRAIN_CFG["lr_factor"],
    }
)

logger.info("Generating integrated embeddings...")
adata.obsm['X_emb'] = model.get_latent(adata, mean=True)

logger.info("Processing integrated representation (scPoli)...")
sc.pp.neighbors(adata, use_rep='X_emb', key_added='neighbors_scpoli', n_neighbors=PCA_CFG["n_neighbors"], random_state=PCA_CFG["random_state"])
sc.tl.umap(adata, neighbors_key='neighbors_scpoli', random_state=PCA_CFG["random_state"])
adata.obsm['X_umap_scpoli'] = adata.obsm['X_umap'].copy()

logger.info("Saving final sample-integrated reference object...")
adata.write_h5ad(FINAL_H5AD_OUT, compression='gzip')
logger.info(f"H5AD written to: {FINAL_H5AD_OUT}")

logger.info("Rendering comparison graphics...")

distinct_palette = sc.pl.palettes.default_102
level_1_colors = PLOT_CFG["level_1_colors"]

def generate_clean_comparison_plot(color_key, filename, title_suffix, palette=None, legend_loc='right margin'):
    fig, axes = plt.subplots(1, 2, figsize=tuple(PLOT_CFG["comparison_figsize"]))

    sc.pl.embedding(
        adata, basis='X_umap_pca', color=color_key, ax=axes[0], legend_loc=None,
        show=False, frameon=False, palette=palette, title=f"BEFORE: PCA ({title_suffix})"
    )
    sc.pl.embedding(
        adata, basis='X_umap_scpoli', color=color_key, ax=axes[1], legend_loc=legend_loc,
        show=False, frameon=False, palette=palette, title=f"AFTER: scPoli ({title_suffix})"
    )

    plt.savefig(os.path.join(OUTPUT_DIR, filename), dpi=PLOT_CFG["dpi"], bbox_inches='tight')
    plt.close()

generate_clean_comparison_plot('sample_id', 'Comparison_Sample_Alignment.png', 'Sample ID', legend_loc=None)
generate_clean_comparison_plot('study', 'Comparison_Study_Alignment.png', 'Study', palette=distinct_palette)
generate_clean_comparison_plot('Level_2', 'Comparison_Level2_Biology.png', 'Level 2 Cell Types', palette=distinct_palette)
generate_clean_comparison_plot('Level_1', 'Comparison_Level1_Biology.png', 'Level 1 Broad Categories', palette=level_1_colors)

logger.info("Master integration pipeline complete!")
