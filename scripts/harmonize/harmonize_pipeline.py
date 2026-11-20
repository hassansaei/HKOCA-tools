#!/usr/bin/env python3
"""Data loaders + gene harmonisation for 21 atlas datasets (scRNA-seq).

Config-driven port of data_loader_harmonized_21datasets.ipynb.

Usage
-----
    python data_loader_harmonized_21datasets.py
    python data_loader_harmonized_21datasets.py --config path/to.config
    python data_loader_harmonized_21datasets.py --studies STUDY1,STUDY2
"""
from __future__ import annotations

import os
import re
import glob
import gzip
import sys
import argparse
import configparser
import platform
import warnings
import urllib.request

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import anndata
import scanpy as sc

try:
    import sfaira  # noqa: F401
except ImportError:
    sfaira = None

warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)
anndata.settings.allow_write_nullable_strings = True

def geo_gsm_url(gsm, filename):
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{filename}"

def geo_gse_url(gse, filename):
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/suppl/{filename}"

def redownload_urls(urls, target_dir):
    os.makedirs(target_dir, exist_ok=True)
    for url in urls:
        fname = os.path.basename(url)
        fpath = os.path.join(target_dir, fname)
        if os.path.exists(fpath):
            os.replace(fpath, fpath + '.corrupt')
        try:
            print(f"  Re-downloading {fname}...")
            urllib.request.urlretrieve(url, fpath)
        except Exception as e:
            print(f"  Failed: {e}")

def read_10x_mtx_with_redownload(path, prefix, urls):
    try:
        return sc.read_10x_mtx(path, var_names='gene_symbols', prefix=prefix)
    except Exception as e:
        print(f"  Read failed ({e}), re-downloading...")
        redownload_urls(urls, path)
        return sc.read_10x_mtx(path, var_names='gene_symbols', prefix=prefix)

def read_table_with_redownload(path, filename, sep, compression, index_col, urls):
    fpath = os.path.join(path, filename)
    try:
        return pd.read_csv(fpath, sep=sep, index_col=index_col, compression=compression)
    except Exception as e:
        print(f"  Read failed ({e}), re-downloading...")
        redownload_urls(urls, path)
        return pd.read_csv(fpath, sep=sep, index_col=index_col, compression=compression)

def _decode(arr):
    return [x.decode('utf-8') if isinstance(x, (bytes, bytearray)) else str(x) for x in arr]

def read_10x_h5_safe(path):
    try:
        return sc.read_10x_h5(path)
    except Exception as err:
        import h5py
        from scipy.sparse import csr_matrix
        with h5py.File(path, 'r') as f:
            g = f['matrix'] if 'matrix' in f else f[list(f.keys())[0]]
            barcodes   = _decode(g['barcodes'][()])
            gene_names = None
            if 'features' in g:
                feat = g['features']
                if 'name'        in feat: gene_names = _decode(feat['name'][()])
                elif 'gene_names' in feat: gene_names = _decode(feat['gene_names'][()])
                elif 'id'        in feat: gene_names = _decode(feat['id'][()])
            if gene_names is None and 'gene_names' in g: gene_names = _decode(g['gene_names'][()])
            if gene_names is None and 'genes'      in g: gene_names = _decode(g['genes'][()])
            if gene_names is None:
                raise RuntimeError(f'No gene names in {path}')
            X = csr_matrix((g['data'][()], g['indices'][()], g['indptr'][()]),
                           shape=g['shape'][()]).T
        return sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes),
                          var=pd.DataFrame(index=gene_names))

def read_10x_h5_with_redownload(path, gsm):
    try:
        return read_10x_h5_safe(path)
    except Exception:
        url = geo_gsm_url(gsm, os.path.basename(path))
        if os.path.exists(path): os.replace(path, path + '.corrupt')
        urllib.request.urlretrieve(url, path)
        return read_10x_h5_safe(path)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG = os.path.join(SCRIPT_DIR, "data_loader_harmonized_21datasets.config")

CFG = {}
DATA_ROOT = None
OUTPUT_ROOT = None
HARMONIZED_DIR = None
H5AD_DIR = None
GTF_FILE = None
ALLOWED_GENES = None
STANDARD_META_COLS = None
META_FILL = "unknown"
REFUSE_OVERWRITE = True


def _split_csv(value: str) -> list[str]:
    return [x.strip() for x in value.replace("\n", ",").split(",") if x.strip()]


def load_config(path: str) -> dict:
    """Load INI config into CFG and module-level path globals."""
    global CFG, DATA_ROOT, OUTPUT_ROOT, HARMONIZED_DIR, H5AD_DIR, GTF_FILE
    global STANDARD_META_COLS, META_FILL, REFUSE_OVERWRITE

    cp = configparser.ConfigParser()
    cp.optionxform = str  # keep key case (Age, etc.)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Config not found: {path}")
    cp.read(path)

    CFG = {s: dict(cp.items(s)) for s in cp.sections()}
    DATA_ROOT = cp.get("paths", "data_root")
    OUTPUT_ROOT = cp.get("paths", "output_root")
    HARMONIZED_DIR = os.path.join(OUTPUT_ROOT, cp.get("paths", "harmonized_subdir"))
    H5AD_DIR = os.path.join(OUTPUT_ROOT, cp.get("paths", "h5ad_subdir"))
    GTF_FILE = cp.get("paths", "gtf_file")
    STANDARD_META_COLS = _split_csv(cp.get("metadata", "standard_cols"))
    META_FILL = cp.get("metadata", "fill_value", fallback="unknown")
    REFUSE_OVERWRITE = cp.getboolean("io", "refuse_overwrite", fallback=True)

    os.makedirs(HARMONIZED_DIR, exist_ok=True)
    os.makedirs(H5AD_DIR, exist_ok=True)
    return CFG


def build_allowed_genes(gtf_file: str | None = None) -> set[str]:
    """Build reference gene set from GTF biotypes listed in config."""
    global ALLOWED_GENES
    gtf_file = gtf_file or GTF_FILE
    biotypes = set(_split_csv(CFG["reference"]["gene_biotypes"]))
    genes = []
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            info = {
                kv.split(" ")[0]: kv.split(" ")[1].strip('";')
                for kv in cols[8].split("; ")
                if " " in kv
            }
            if info.get("gene_biotype") in biotypes:
                genes.append(info["gene_name"])
    ALLOWED_GENES = set(genes)
    print(f"Reference genes ({', '.join(sorted(biotypes))}): {len(ALLOWED_GENES)}")
    return ALLOWED_GENES


def enabled_studies(cli_studies: str | None = None) -> list[str]:
    atlas = _split_csv(CFG["studies"]["list"])
    if cli_studies:
        wanted = _split_csv(cli_studies)
    else:
        enabled = CFG["studies"].get("enabled", "all").strip()
        wanted = atlas if enabled.lower() == "all" else _split_csv(enabled)
    unknown = [s for s in wanted if s not in atlas]
    if unknown:
        raise ValueError(f"Unknown study keys (not in config [studies] list): {unknown}")
    return wanted


def maybe_write_h5ad(adata, path: str, label: str = "h5ad") -> None:
    if os.path.exists(path):
        if REFUSE_OVERWRITE:
            print(f"  {label} exists, not overwriting: {path}")
            return
        raise FileExistsError(f"Refusing to overwrite ({label}): {path}")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    adata.write(path)
    print(f"  Saved {label}: {path} {getattr(adata, 'shape', '')}")


def fill_missing_meta(adata, fill_value=None):
    """Ensure all STANDARD_META_COLS exist in adata.obs."""
    fill_value = META_FILL if fill_value is None else fill_value
    for col in STANDARD_META_COLS:
        if col not in adata.obs.columns:
            adata.obs[col] = fill_value
    for col in STANDARD_META_COLS:
        adata.obs[col] = adata.obs[col].astype(str)
    return adata


def memory_efficient_harmonize_sparse(input_file, output_file):
    """Harmonise gene space to ALLOWED_GENES; preserve obs; no QC."""
    if os.path.exists(output_file) and REFUSE_OVERWRITE:
        print(f"  Harmonised exists, not overwriting: {output_file}")
        return
    if os.path.exists(output_file):
        raise FileExistsError(f"Refusing to overwrite: {output_file}")
    print(f"  Harmonising {os.path.basename(input_file)} [sparse mode]...")
    adata = sc.read_h5ad(input_file)
    print(f"    Input  shape : {adata.shape}")
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata_new = harmonize_matrix_sparse(adata, ALLOWED_GENES)
    print(f"    Output shape : {adata_new.shape}")
    adata_new.write(output_file)
    print(f"    Saved -> {output_file}")

def harmonize_matrix(df):
    """
    Align a genes×cells DataFrame to the reference ALLOWED_GENES set.
    Common genes are retained; missing genes are added as sparse zero rows.
    Returns cells×genes.
    """
    dataset_genes  = set(df.index)
    common_genes   = sorted(ALLOWED_GENES.intersection(dataset_genes))
    missing_genes  = sorted(ALLOWED_GENES - dataset_genes)
    df = df.loc[common_genes, :]
    zero_df = pd.DataFrame.sparse.from_spmatrix(
        sp.csr_matrix((len(missing_genes), df.shape[1])),
        index=missing_genes, columns=df.columns)
    return pd.concat([df, zero_df]).loc[sorted(ALLOWED_GENES), :].T

def memory_efficient_harmonize(input_file, output_file, harmonize_matrix):
    """
    Load AnnData, harmonise gene space via pandas sparse, save.
    All obs metadata is preserved. No QC filtering.
    """
    adata = sc.read_h5ad(input_file)
    if sp.issparse(adata.X):
        df = pd.DataFrame.sparse.from_spmatrix(
            adata.X, index=adata.obs_names, columns=adata.var_names)
    else:
        df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    df = harmonize_matrix(df.T)
    adata_new = sc.AnnData(X=df, obs=adata.obs.copy(), var=pd.DataFrame(index=df.columns))
    adata_new.write(output_file)
    print(f"  Saved harmonised → {output_file}  shape={adata_new.shape}")

def harmonize_matrix_sparse(adata, ALLOWED_GENES):
    """
    Harmonise an AnnData using scipy sparse operations only - no dense conversion.
    Recommended for large datasets (> 200 K cells).
    """
    target_genes  = sorted(list(ALLOWED_GENES))
    dataset_genes = adata.var_names.tolist()
    common_genes  = sorted(list(set(dataset_genes) & set(target_genes)))
    if dataset_genes == target_genes:
        return adata
    adata_c = adata[:, common_genes].copy()
    g2i = {g: i for i, g in enumerate(target_genes)}
    coo = sp.coo_matrix(adata_c.X if sp.issparse(adata_c.X) else adata_c.X)
    l2g = np.array([g2i[g] for g in common_genes])
    new_X = sp.coo_matrix(
        (coo.data, (coo.row, l2g[coo.col])),
        shape=(adata.shape[0], len(target_genes))
    ).tocsr()
    return sc.AnnData(X=new_X, obs=adata.obs.copy(), var=pd.DataFrame(index=target_genes))


def load_d10_1016_j_cell_2020_04_004():
    """Load and harmonise d10_1016_j_cell_2020_04_004."""
    # 1 - GSE147863 | d10_1016_j_cell_2020_04_004
    study   = 'd10_1016_j_cell_2020_04_004'
    h5_path = os.path.join(DATA_ROOT, study, 'GSE147863_RAW',
                           'GSM4447249_KidneyOrganoid_FilteredGeneBCMatrices.h5')

    adata = read_10x_h5_with_redownload(h5_path, 'GSM4447249')
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    adata.obs['sample_id']     = 'GSM4447249'
    adata.obs['study']         = study
    adata.obs['source']        = 'hiPSC'
    adata.obs['diff_protocol'] = 'Morizane_modified'
    adata.obs['sc_protocol']   = '10x_3_v3'
    adata.obs['sequencing']    = 'Illumina_NextSeq_500'
    adata.obs['genome_build']  = 'hg38'
    adata.obs['Age']           = 'd16'
    adata.obs['type']          = 'scRNAseq'
    adata.obs['disease']       = 'Control'
    adata.obs['condition']     = 'WT'
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1016_j_celrep_2020_108514():
    """Load and harmonise d10_1016_j_celrep_2020_108514."""
    # 2 - GSE131086 | d10_1016_j_celrep_2020_108514
    study     = 'd10_1016_j_celrep_2020_108514'
    file_path = os.path.join(DATA_ROOT, study, 'GSE131086_RAW')

    protocol_mapping = {
        'GSM3763146': 'Takasato',
        'GSM3763147': 'Uchimura'
    }

    adatas = []
    for names in os.listdir(file_path):
        df    = pd.read_csv(os.path.join(file_path, names), sep='\t', index_col=0, compression='gzip')
        adata = sc.AnnData(df.T)

        gsm_id = names.split('_')[0]

        adata.obs['sample_id'] = gsm_id
        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        adata.obs['study']         = study
        adata.obs['source']        = 'hiPSC'
        adata.obs['diff_protocol'] = protocol_mapping.get(gsm_id, 'unknown')
        adata.obs['sc_protocol']   = '10X_5'
        adata.obs['sequencing']    = 'Illumina_NovaSeq_6000'
        adata.obs['genome_build']  = 'hg38'
        adata.obs['Age']           = 'd26'
        adata.obs['type']          = 'scRNAseq'
        adata.obs['disease']       = 'Control'
        adatas.append(adata)

    adata = sc.concat(adatas, join='outer')
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1016_j_cmet_2022_04_009():
    """Load and harmonise d10_1016_j_cmet_2022_04_009."""
    # 3 - GSE181002 | d10_1016_j_cmet_2022_04_009
    # CONTROL ONLY: 5mM_Mock (excluded 5mM_SARS - COVID-19 treatment)
    study   = 'd10_1016_j_cmet_2022_04_009'
    h5_path = os.path.join(DATA_ROOT, study, 'GSE181002_RAW',
                           'GSM5481909_5mM_Mock_filtered_feature_bc_matrix.h5')

    adata = sc.read_10x_h5(h5_path)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    adata.obs['sample_id']     = 'GSM5481909'
    adata.obs['study']         = study
    adata.obs['source']        = 'hiPSC'
    adata.obs['diff_protocol'] = 'Morizane_modified'
    adata.obs['sc_protocol']   = '10x_3_v3.1'
    adata.obs['sequencing']    = 'Illumina_NovaSeq_6000'
    adata.obs['genome_build']  = 'hg38'
    adata.obs['Age']           = 'd16'
    adata.obs['type']          = 'scRNAseq'
    adata.obs['disease']       = 'Control'
    adata.obs['condition']     = '5mM_Mock'
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1016_j_devcel_2019_06_001():
    """Load and harmonise d10_1016_j_devcel_2019_06_001."""
    # 4 - GSE124472 | d10_1016_j_devcel_2019_06_001
    import tarfile

    study    = 'd10_1016_j_devcel_2019_06_001'
    basepath = os.path.join(DATA_ROOT, study, 'GSE124472_RAW')

    adatas = []
    for archive_name in os.listdir(basepath):
        if not archive_name.endswith('.tar.gz'):
            continue

        gsm_id = archive_name.split('_')[0]
        extract_dir = os.path.join(basepath, gsm_id)

        if not os.path.exists(extract_dir):
            print(f'  Extracting {archive_name}...')
            os.makedirs(extract_dir, exist_ok=True)
            with tarfile.open(os.path.join(basepath, archive_name), 'r:gz') as tar:
                tar.extractall(path=extract_dir)

        mtx_dir = extract_dir
        for root, dirs, files in os.walk(extract_dir):
            if any('matrix.mtx' in f for f in files):
                mtx_dir = root
                break

        adata = sc.read_10x_mtx(mtx_dir)
        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        adata.obs['sample_id']     = gsm_id
        adata.obs['study']         = study
        adata.obs['source']        = 'hiPSC'
        adata.obs['diff_protocol'] = 'Morizane_modified' # Maintained from your prompt
        adata.obs['sc_protocol']   = '10x_3_v2'
        adata.obs['sequencing']    = 'Illumina_NextSeq_500'
        adata.obs['genome_build']  = 'hg38'
        adata.obs['type']          = 'scRNAseq'
        adata.obs['disease']       = 'Control'
        adatas.append(adata)
        print(f'  Loaded {gsm_id}')

    adata = sc.concat(adatas, join='outer')
    mapping = {'GSM3534660': 'd16', 'GSM3534661': 'd16',
               'GSM3534662': 'd28', 'GSM3534663': 'd28'}
    adata.obs['Age'] = adata.obs['sample_id'].map(mapping)
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1016_j_stem_2018_04_022():
    """Load and harmonise d10_1016_j_stem_2018_04_022."""
    # 5 - GSE109718 | d10_1016_j_stem_2018_04_022
    study     = 'd10_1016_j_stem_2018_04_022'
    file_path = os.path.join(DATA_ROOT, study, 'GSE109718_RAW')

    adatas = []
    for names in os.listdir(file_path):
        if not names.endswith('.gz'):
            continue

        df    = pd.read_csv(os.path.join(file_path, names), sep='\t', index_col=0, compression='gzip')
        adata = sc.AnnData(df.T)

        gsm_id = names.split('_')[0]
        adata.obs['sample_id'] = gsm_id

        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        adata.obs['study']         = study
        adata.obs['source']        = 'hiPSC-ESC'
        adata.obs['diff_protocol'] = 'Freedman'
        adata.obs['sc_protocol']   = 'Drop-seq'
        adata.obs['sequencing']    = 'Illumina_HiSeq_2000'
        adata.obs['genome_build']  = 'hg38'
        adata.obs['Age']           = 'd21'
        adata.obs['type']          = 'scRNAseq'
        adata.obs['disease']       = 'Control'

        if 'NoVegf' in names:
            adata.obs['condition'] = 'NoVegf'
        else:
            adata.obs['condition'] = 'WT'

        adatas.append(adata)

    adata = sc.concat(adatas, join='outer')
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1016_j_stem_2018_10_010():
    """Load and harmonise d10_1016_j_stem_2018_10_010."""
    # 6 - GSE118184 | d10_1016_j_stem_2018_10_010
    import os
    import pandas as pd
    import scipy.sparse as sp
    import scanpy as sc
    import anndata as ad

    study       = 'd10_1016_j_stem_2018_10_010'
    folder_path = os.path.join(DATA_ROOT, study, 'GSE118184')

    print(f"Loading files for {study}...")

    adata_list = []
    for fname in os.listdir(folder_path):
        if not fname.endswith('.gz'):
            continue

        # Exclude the BDNF inhibitor treated sample
        if 'BDNF.inhibitor' in fname:
            continue

        print(f"  -> Processing: {fname}")

        df  = pd.read_csv(os.path.join(folder_path, fname), sep='\t', index_col=0, compression='gzip').T
        X   = sp.csr_matrix(df.values)

        tmp = sc.AnnData(X=X)
        tmp.obs_names = df.index
        tmp.var_names = df.columns

        tmp.obs_names_make_unique()
        tmp.var_names_make_unique()

        tmp.obs['sample_id'] = fname
        adata_list.append(tmp)

    adata = ad.concat(adata_list, join='outer', index_unique='_')

    adata.obs['sc_protocol']   = 'Drop-seq'
    adata.obs['sequencing']    = 'Illumina_HiSeq_2500'
    adata.obs['genome_build']  = 'hg38'
    adata.obs['type']          = 'scRNAseq'
    adata.obs['disease']       = 'Control'
    adata.obs['condition']     = 'Control'
    adata.obs['study']         = study

    age_map = {
        'GSE118184_Takasato.iPS.timecourse.txt.gz':   ('Unknown', 'Takasato', 'iPSC'),

        'GSE118184_Takasato_ES.dge.txt.gz':           ('d26', 'Takasato', 'ESC'),
        'GSE118184_Takasato_iPS_Batch1_2.dge.txt.gz': ('d26', 'Takasato', 'iPSC'),
        'GSE118184_Takasato_iPS_Batch3.dge.txt.gz':   ('d26', 'Takasato', 'iPSC'),
        'GSE118184_Takasato_iPS_day34.dge.txt.gz':    ('d34', 'Takasato', 'iPSC'),

        'GSE118184_Morizane_ES.dge.txt.gz':           ('d26', 'Morizane', 'ESC'),
        'GSE118184_Morizane_iPS_Batch1_2.dge.txt.gz': ('d26', 'Morizane', 'iPSC'),
        'GSE118184_Morizane_iPS_Batch3.dge.txt.gz':   ('d26', 'Morizane', 'iPSC'),
    }

    adata.obs['Age']           = adata.obs['sample_id'].map(lambda x: age_map.get(x, ('Unknown', 'Unknown', 'Unknown'))[0])
    adata.obs['diff_protocol'] = adata.obs['sample_id'].map(lambda x: age_map.get(x, ('Unknown', 'Unknown', 'Unknown'))[1])
    adata.obs['source']        = adata.obs['sample_id'].map(lambda x: age_map.get(x, ('Unknown', 'Unknown', 'Unknown'))[2])

    adata.obs['sample_id'] = adata.obs['sample_id'].str.replace('.dge.txt.gz', '', regex=False).str.replace('.txt.gz', '', regex=False)

    extracted_days = adata.obs.index.to_series().str.extract(r'(Day\d+)')[0]

    mask = extracted_days.notna()
    if mask.any():
        adata.obs.loc[mask, 'Age'] = extracted_days[mask]

    adata.obs['Age'] = adata.obs['Age'].astype(str).str.replace(r'Day', 'd', regex=True)

    for col in ['sample_id', 'sc_protocol', 'sequencing', 'genome_build', 'type', 'disease', 'condition', 'study', 'Age', 'diff_protocol', 'source']:
        adata.obs[col] = adata.obs[col].astype('category')

    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out, compression='gzip')
    print(f'✅ Saved {study}: {adata.n_obs} cells x {adata.n_vars} genes')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1016_j_stem_2019_06_009():
    """Load and harmonise d10_1016_j_stem_2019_06_009."""
    # 7 - GSE132023 | d10_1016_j_stem_2019_06_009
    study      = 'd10_1016_j_stem_2019_06_009'
    count_path = os.path.join(
        DATA_ROOT, study, 'raw', 'GSE132023_Day10_12_14_count.txt',
        'day10_12_14_GEO_trawdata_barcode_v2_wHEADER_sort-trans_edit.txt'
    )

    df   = pd.read_csv(count_path, sep='\t', index_col=0, low_memory=False)
    meta = df.iloc[0]
    df   = df.iloc[1:].apply(pd.to_numeric, errors='coerce').fillna(0)
    X    = sp.csr_matrix(df.to_numpy(dtype=np.float32).T)
    adata = sc.AnnData(X=X, obs=pd.DataFrame(index=df.columns),
                       var=pd.DataFrame(index=df.index))
    adata.obs['sample_id'] = meta.values
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    mapping1 = {'Day10rep1': 'GSM3834523', 'Day10rep2': 'GSM3834524',
                'Day12':     'GSM3834525', 'Day14rep1': 'GSM3834526',
                'Day14rep2': 'GSM3834527'}
    mapping2 = {'Day10rep1': 'd10', 'Day10rep2': 'd10', 'Day12': 'd12',
                'Day14rep1': 'd14', 'Day14rep2': 'd14'}
    adata.obs['sample']        = adata.obs['sample_id'].map(mapping1).fillna('unknown')
    adata.obs['Age']           = adata.obs['sample_id'].map(mapping2).fillna('unknown')
    adata.obs['study']         = study
    adata.obs['source']        = 'hPSCs'
    adata.obs['diff_protocol'] = 'Low'
    adata.obs['sc_protocol']   = '10x_3_v1'
    adata.obs['sequencing']    = 'Illumina_HiSeq_4000'
    adata.obs['genome_build']  = 'hg38'
    adata.obs['type']          = 'scRNAseq'
    adata.obs['disease']       = 'Control'
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1016_j_stem_2021_11_001():
    """Load and harmonise d10_1016_j_stem_2021_11_001."""
    # 8 - GSE181595 | d10_1016_j_stem_2021_11_001
    # 7 samples: d21×2, d28×3, d35×2 - all controls
    study    = 'd10_1016_j_stem_2021_11_001'
    basepath = os.path.join(DATA_ROOT, study, 'GSE181595_RAW')

    h5_files = {
        'GSM5507784_3D_day21_1_filtered_gene_bc_matrices_h5.h5': {'sample_id': 'GSM5507784', 'Age': 'd21'},
        'GSM5507785_3D_day21_2_filtered_gene_bc_matrices_h5.h5': {'sample_id': 'GSM5507785', 'Age': 'd21'},
        'GSM5507786_3D_day28_1_filtered_gene_bc_matrices_h5.h5': {'sample_id': 'GSM5507786', 'Age': 'd28'},
        'GSM5507787_3D_day28_2_filtered_gene_bc_matrices_h5.h5': {'sample_id': 'GSM5507787', 'Age': 'd28'},
        'GSM5507788_3D_day28_3_filtered_gene_bc_matrices_h5.h5': {'sample_id': 'GSM5507788', 'Age': 'd28'},
        'GSM5507789_3D_day35_1_filtered_gene_bc_matrices_h5.h5': {'sample_id': 'GSM5507789', 'Age': 'd35'},
        'GSM5507790_3D_day35_2_filtered_gene_bc_matrices_h5.h5': {'sample_id': 'GSM5507790', 'Age': 'd35'},
    }

    adatas = []
    for fname, meta in h5_files.items():
        adata_sample = read_10x_h5_with_redownload(os.path.join(basepath, fname), meta['sample_id'])
        adata_sample.obs_names_make_unique()
        adata_sample.var_names_make_unique()

        adata_sample.obs['sample_id']     = meta['sample_id']
        adata_sample.obs['study']         = study
        adata_sample.obs['source']        = 'hiPSC'
        adata_sample.obs['diff_protocol'] = 'Takasato_modified'
        adata_sample.obs['sc_protocol']   = '10x_3_v3.1'
        adata_sample.obs['sequencing']    = 'Illumina_HiSeq_2500'
        adata_sample.obs['genome_build']  = 'hg38'
        adata_sample.obs['Age']           = meta['Age']  # Exclusively mapped to 'Age'
        adata_sample.obs['type']          = 'scRNAseq'
        adata_sample.obs['disease']       = 'Control'

        adatas.append(adata_sample)

    adata = sc.concat(adatas, join='outer')
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1016_j_stem_2022_06_005():
    """Load and harmonise d10_1016_j_stem_2022_06_005."""
    # 9 - GSE164564 | d10_1016_j_stem_2022_06_005
    import tarfile
    import gzip
    import shutil
    import os
    import pandas as pd
    import scanpy as sc

    study    = 'd10_1016_j_stem_2022_06_005'
    basepath = os.path.join(DATA_ROOT, study, 'GSE164564_RAW')

    adatas = []
    for archive_name in os.listdir(basepath):
        if not archive_name.endswith('.tar.gz'):
            continue

        gsm_id = archive_name.split('_')[0]
        extract_dir = os.path.join(basepath, gsm_id)
        os.makedirs(extract_dir, exist_ok=True)

        print(f'  Processing {gsm_id}...')
        file_path = os.path.join(basepath, archive_name)
        current_target = file_path

        layer = 0
        while True:
            with open(current_target, 'rb') as f:
                magic = f.read(2)
            if magic == b'\x1f\x8b':
                layer += 1
                print(f'    Decompressing gzip layer {layer}...')
                temp_out = os.path.join(extract_dir, f'unzipped_L{layer}.tmp')
                with gzip.open(current_target, 'rb') as f_in:
                    with open(temp_out, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                current_target = temp_out
            else:
                break

        if tarfile.is_tarfile(current_target):
            print('    Payload is a tar archive. Extracting...')
            with tarfile.open(current_target, 'r') as tar:
                tar.extractall(path=extract_dir)
        else:
            print('    Payload is not a tar file. Assuming it is a direct counts matrix.')
            final_csv = os.path.join(extract_dir, 'counts_matrix.csv')
            shutil.copy(current_target, final_csv)

        for f in os.listdir(extract_dir):
            if f.endswith('.tmp'):
                os.remove(os.path.join(extract_dir, f))

        try:
            adata_sample = None
            mtx_dir = None

            for root, dirs, files in os.walk(extract_dir):
                if any(f.endswith('.mtx') or f.endswith('.mtx.gz') for f in files):
                    mtx_dir = root
                    for f in files:
                        old_path = os.path.join(root, f)
                        fl = f.lower()
                        new_name = None
                        if fl.endswith('.mtx.gz'): new_name = 'matrix.mtx.gz'
                        elif fl.endswith('.mtx'): new_name = 'matrix.mtx'
                        elif 'barcode' in fl: new_name = 'barcodes.tsv.gz' if fl.endswith('.gz') else 'barcodes.tsv'
                        elif 'feature' in fl or 'gene' in fl: new_name = 'features.tsv.gz' if fl.endswith('.gz') else 'features.tsv'

                        if new_name and new_name != f:
                            os.rename(old_path, os.path.join(root, new_name))
                    break

            if mtx_dir:
                print('    Loading via scanpy MTX reader...')
                adata_sample = sc.read_10x_mtx(mtx_dir, var_names='gene_symbols')

            elif os.path.exists(os.path.join(extract_dir, 'counts_matrix.csv')):
                print('    Loading as a dense tabular matrix...')
                csv_path = os.path.join(extract_dir, 'counts_matrix.csv')

                with open(csv_path, 'r') as f:
                    first_line = f.readline()
                sep = '\t' if '\t' in first_line else ','

                df = pd.read_csv(csv_path, sep=sep, index_col=0)

                if df.shape[0] > df.shape[1] or df.shape[0] > 15000:
                    df = df.T
                adata_sample = sc.AnnData(df)

            if adata_sample is not None:
                adata_sample.obs_names_make_unique()
                adata_sample.var_names_make_unique()

                adata_sample.obs['sample_id']     = gsm_id
                adata_sample.obs['study']         = study
                adata_sample.obs['source']        = 'hiPSC'
                adata_sample.obs['diff_protocol'] = 'Morizane_modified'
                adata_sample.obs['sc_protocol']   = '10x_3_v3'
                adata_sample.obs['sequencing']    = 'Illumina_NextSeq_500'
                adata_sample.obs['genome_build']  = 'hg38'
                adata_sample.obs['type']          = 'scRNAseq'
                adata_sample.obs['disease']       = 'Control'

                adatas.append(adata_sample)
                print(f'    -> Successfully loaded {gsm_id}: {adata_sample.shape}')
            else:
                print(f"  [❌] No valid MTX or tabular data found in {gsm_id}")

        except Exception as e:
            print(f"  [❌] Failed to load matrix for {gsm_id}: {e}")

    if adatas:
        adata = sc.concat(adatas, join='outer')
        mapping = {'GSM5014155': 'd8', 'GSM5014156': 'd10', 'GSM5014157': 'd14',
                   'GSM5014158': 'd16', 'GSM5014159': 'd28'}
        adata.obs['Age'] = adata.obs['sample_id'].map(mapping)
        fill_missing_meta(adata)

        raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
        os.makedirs(os.path.dirname(raw_out), exist_ok=True)
        adata.write(raw_out)
        print(f'Saved {study}: {adata.shape}')

        memory_efficient_harmonize_sparse(
            input_file=raw_out,
            output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
        )
    else:
        print(f"No valid datasets were loaded for {study}.")

def load_d10_1038_s41467_019_13382_0():
    """Load and harmonise d10_1038_s41467_019_13382_0."""
    # 10 - GSE136314 | d10_1038_s41467_019_13382_0
    # Uses CellRanger output layout; 26 folders mapped to 9 GSM samples
    study          = 'd10_1038_s41467_019_13382_0'
    CELLRANGER_DIR = os.path.join(DATA_ROOT, study, 'cellranger')

    raw_gse = os.path.join(DATA_ROOT, study, 'GSE136314_RAW')
    if os.path.isdir(raw_gse):
        for fpath in glob.glob(os.path.join(raw_gse, '*_features.tsv.gz')):
            df_feat = pd.read_csv(fpath, sep='\t', header=None)
            if df_feat.shape[1] < 3:
                df_feat['_dup'] = df_feat.iloc[:, 0]
                df_feat['Gene Expression'] = 'Gene Expression'
                with gzip.open(fpath, 'wt', encoding='utf-8') as fh:
                    df_feat.to_csv(fh, sep='\t', header=False, index=False, lineterminator='\n')

    gsm_mapping = {
        'Org_Expt1_Day29_1_AS': 'GSM4044536', 'Org_Expt1_Day29_3_AS': 'GSM4044536',
        'Org_Expt2_Day29_1_AS': 'GSM4044536', 'Org_Expt2_Day29_2_AS': 'GSM4044536',
        'Org_Expt2_Day29_3_AS': 'GSM4044536',
        'Org_Expt1_Day15_1_AS_1': 'GSM4044537', 'Org_Expt1_Day15_1_AS_2': 'GSM4044537',
        'Org_Expt1_Day15_1_AS_3': 'GSM4044537', 'Org_Expt2_Day15_1_AS_1': 'GSM4044537',
        'Org_Expt2_Day15_1_AS_2': 'GSM4044537', 'Org_Expt2_Day15_1_AS_3': 'GSM4044537',
        'Org_Expt1_Day7_1_AS':  'GSM4044538', 'Org_Expt2_Day7_2_AS': 'GSM4044538',
        'Org_Expt2_Day29_1_ThF': 'GSM4044540', 'Org_Expt2_Day29_2_ThF': 'GSM4044540',
        'Org_Expt2_Day29_3_ThF': 'GSM4044540',
        'Org_Expt2_Day15_1_ThF_1': 'GSM4044541', 'Org_Expt2_Day15_1_ThF_2': 'GSM4044541',
        'Org_Expt2_Day15_1_ThF_3': 'GSM4044541',
        'Org_Expt1_Day7_1_ThF': 'GSM4044542', 'Org_Expt2_Day7_1_ThF': 'GSM4044542',
        'ThF_ipsc': 'GSM4044543',
        'Org_Day32_ThF_1': 'GSM4044544', 'Org_Day32_ThF_2': 'GSM4044544',
        'Org_Day32_ThF_3': 'GSM4044544', 'Org_Day32_ThF_4': 'GSM4044544',
    }
    age_mapping = {
        'GSM4044536': 'd29', 'GSM4044537': 'd15', 'GSM4044538': 'd7',
        'GSM4044539': 'd0',  'GSM4044540': 'd29', 'GSM4044541': 'd15',
        'GSM4044542': 'd7',  'GSM4044543': 'd0',  'GSM4044544': 'd32',
    }

    triplet_paths = glob.glob(os.path.join(CELLRANGER_DIR, '*/filtered_feature_bc_matrix'))
    print(f'Found {len(triplet_paths)} CellRanger folders')
    adatas = []
    for path in triplet_paths:
        folder_name = path.split('/')[-2]
        tmp = sc.read_10x_mtx(path, var_names='gene_symbols')
        tmp.obs_names = f'{folder_name}-' + tmp.obs_names
        gsm_id  = gsm_mapping.get(folder_name, 'unknown')
        age_val = age_mapping.get(gsm_id, 'unknown')
        sc_prot = '10x_3_v3.1' if gsm_id == 'GSM4044544' else '10x_3_v2'
        tmp.obs['sample_id']        = gsm_id
        tmp.obs['folder_replicate'] = folder_name
        tmp.obs['Age']              = age_val
        tmp.obs['study']            = study
        tmp.obs['source']           = 'hiPSC'
        tmp.obs['diff_protocol']    = 'Takasato'
        tmp.obs['sc_protocol']      = sc_prot
        tmp.obs['sequencing']       = 'HiSeq_X_Ten'
        tmp.obs['genome_build']     = 'hg38'
        tmp.obs['type']             = 'scRNAseq'
        tmp.obs['disease']          = 'Control'
        adatas.append(tmp)

    adata = sc.concat(adatas, join='outer', index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1038_s41467_022_33623_z():
    """Load and harmonise d10_1038_s41467_022_33623_z."""
    # 11 - GSE184928 | d10_1038_s41467_022_33623_z
    study    = 'd10_1038_s41467_022_33623_z'
    basepath = os.path.join(DATA_ROOT, study, 'GSE184928_RAW')

    adatas = []
    for folder_name in os.listdir(basepath):
        folder_path = os.path.join(basepath, folder_name)
        if not os.path.isdir(folder_path):
            continue

        prefix = ""
        for f in os.listdir(folder_path):
            if f.endswith('matrix.mtx.gz'):
                prefix = f.replace('matrix.mtx.gz', '')
                break

        if not prefix:
            print(f"  Warning: No matrix.mtx.gz found in {folder_name}")
            continue

        adata_sample = sc.read_10x_mtx(folder_path, prefix=prefix, var_names='gene_symbols')
        adata_sample.obs_names_make_unique()
        adata_sample.var_names_make_unique()

        adata_sample.obs['sample_id']     = folder_name
        adata_sample.obs['study']         = study
        adata_sample.obs['source']        = 'hiPSC'
        adata_sample.obs['diff_protocol'] = 'Vanslambrouck_2022'
        adata_sample.obs['sc_protocol']   = '10x_3_v3.1'
        adata_sample.obs['sequencing']    = 'Illumina_NovaSeq_6000'
        adata_sample.obs['genome_build']  = 'hg38'
        adata_sample.obs['type']          = 'scRNAseq'
        adata_sample.obs['disease']       = 'Control'
        adatas.append(adata_sample)
        print(f'  Loaded {folder_name} (prefix: {prefix})')

    adata = sc.concat(adatas, join='outer')

    # FIXED MAPPING: Replaced the incorrect GSM7821878 with GSM5600484
    mapping = {'GSM5600482': 'd13', 'GSM5600483': 'd27', 'GSM5600484': 'd25'}
    adata.obs['Age'] = adata.obs['sample_id'].map(mapping)

    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1038_s41551_025_01542_1():
    """Load and harmonise d10_1038_s41551_025_01542_1."""
    # 12 - GSE229433 | d10_1038_s41551_025_01542_1
    study    = 'd10_1038_s41551_025_01542_1'
    basepath = os.path.join(DATA_ROOT, study, 'GSE229433_RAW')

    h5_files = {
        'GSM7163470_filtered_feature_bc_matrix_05K.h5' : {'sample_id': 'GSM7163470', 'condition': '05K'},
        'GSM7163471_filtered_feature_bc_matrix_250K.h5': {'sample_id': 'GSM7163471', 'condition': '250K'},
    }

    adatas = []
    for fname, meta in h5_files.items():
        adata_sample = read_10x_h5_with_redownload(os.path.join(basepath, fname), meta['sample_id'])
        adata_sample.obs_names_make_unique()
        adata_sample.var_names_make_unique()

        adata_sample.obs['sample_id']     = meta['sample_id']
        adata_sample.obs['study']         = study
        adata_sample.obs['source']        = 'hPSC'
        adata_sample.obs['diff_protocol'] = 'Garreta_suspension'
        adata_sample.obs['sc_protocol']   = '10x_3_v3'
        adata_sample.obs['sequencing']    = 'Illumina_NextSeq_500'
        adata_sample.obs['genome_build']  = 'hg38'
        adata_sample.obs['Age']           = 'd16'
        adata_sample.obs['type']          = 'scRNAseq'
        adata_sample.obs['disease']       = 'Control'
        adata_sample.obs['condition']     = meta['condition']

        adatas.append(adata_sample)
        print(f"  Loaded {meta['sample_id']} (Condition: {meta['condition']})")

    adata = sc.concat(adatas, join='outer')
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1038_s41563_020_00853_9():
    """Load and harmonise d10_1038_s41563_020_00853_9."""
    # 13 - GSE152014 | d10_1038_s41563_020_00853_9
    import tarfile

    study    = 'd10_1038_s41563_020_00853_9'
    basepath = os.path.join(DATA_ROOT, study, 'GSE152014_RAW')

    sample_meta = {
        'man': {'sample_id': 'GSM4598695', 'condition': 'Manual_230k'},
        'R0':  {'sample_id': 'GSM4598696', 'condition': 'Bioprinted_R0_110k'},
        'R40': {'sample_id': 'GSM4598697', 'condition': 'Bioprinted_R40_110k'}
    }

    adatas = []
    for archive_name in os.listdir(basepath):
        # Skip non-archives and HTO (Hash-Tag Oligo) cell hashing matrices
        if not archive_name.endswith('.tar.gz') or 'HTO' in archive_name:
            continue

        sample_key = 'unknown'
        for key in sample_meta.keys():
            if f"_{key}_" in archive_name:
                sample_key = key
                break

        if sample_key == 'unknown':
            continue

        meta = sample_meta[sample_key]
        extract_dir = os.path.join(basepath, meta['sample_id'])

        if not os.path.exists(extract_dir):
            print(f'  Extracting {archive_name}...')
            os.makedirs(extract_dir, exist_ok=True)
            try:
                with tarfile.open(os.path.join(basepath, archive_name), 'r:*') as tar:
                    tar.extractall(path=extract_dir)
            except tarfile.ReadError:
                print(f"  [❌] FAILED to extract {archive_name}. File may be corrupted.")
                continue

        mtx_dir = extract_dir
        for root, dirs, files in os.walk(extract_dir):
            if any('matrix.mtx' in f for f in files):
                mtx_dir = root
                break

        try:
            adata = sc.read_10x_mtx(mtx_dir)
            adata.obs_names_make_unique()
            adata.var_names_make_unique()

            adata.obs['sample_id']     = meta['sample_id']
            adata.obs['study']         = study
            adata.obs['source']        = 'hiPSC'
            adata.obs['diff_protocol'] = 'Takasato' 
            adata.obs['sc_protocol']   = '10x_3_v3.1'
            adata.obs['sequencing']    = 'Illumina_NovaSeq_6000'
            adata.obs['genome_build']  = 'hg38'
            adata.obs['Age']           = 'd12'  # From GEO: dissociated after 12 days
            adata.obs['type']          = 'scRNAseq'
            adata.obs['disease']       = 'Control'
            adata.obs['condition']     = meta['condition']

            adatas.append(adata)
            print(f"  Loaded {meta['sample_id']} ({meta['condition']})")
        except Exception as e:
            print(f"  [❌] Failed to load matrix for {archive_name}: {e}")

    if adatas:
        adata = sc.concat(adatas, join='outer')
        fill_missing_meta(adata)

        raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
        os.makedirs(os.path.dirname(raw_out), exist_ok=True)
        adata.write(raw_out)
        print(f'Saved {study}: {adata.shape}')

        memory_efficient_harmonize_sparse(
            input_file=raw_out,
            output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
        )
    else:
        print(f"No valid datasets were loaded for {study}.")

def load_d10_1038_s41592_018_0253_2():
    """Load and harmonise d10_1038_s41592_018_0253_2."""
    # 14 - GSE114802 | d10_1038_s41592_018_0253_2
    import gzip
    from scipy.io import mmread

    study    = 'd10_1038_s41592_018_0253_2'

    basepath = os.path.join(DATA_ROOT, study, 'GSE114802_RAW')
    if not os.path.exists(basepath):
        basepath = os.path.join(DATA_ROOT, study)

    files  = set(['_'.join(x.split('_')[:2]) + '_' for x in os.listdir(basepath) if 'matrix.mtx' in x])

    adatas = []
    for prefix in files:
        matrix_path   = os.path.join(basepath, f'{prefix}matrix.mtx.gz')
        features_path = os.path.join(basepath, f'{prefix}genes.tsv.gz')
        barcodes_path = os.path.join(basepath, f'{prefix}barcodes.tsv.gz')

        if not os.path.exists(features_path):
            features_path = os.path.join(basepath, f'{prefix}features.tsv.gz')

        print(f'  Loading {prefix} manually via scipy mmread...')

        with gzip.open(matrix_path, 'rb') as fh:
            X = mmread(fh).tocsr().T

        features = pd.read_csv(features_path, sep='\t', header=None, compression='gzip')
        barcodes = pd.read_csv(barcodes_path, sep='\t', header=None, compression='gzip')

        gene_ids     = features[0].values
        gene_symbols = features[1].values if features.shape[1] >= 2 else features[0].values

        adata = sc.AnnData(X=X,
                           obs=pd.DataFrame(index=barcodes[0].values),
                           var=pd.DataFrame({'gene_ids': gene_ids}, index=gene_symbols))

        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        clean_sample_id = prefix.rstrip('_')

        adata.obs['sample_id']     = clean_sample_id
        adata.obs['study']         = study
        adata.obs['source']        = 'hiPSC'
        adata.obs['diff_protocol'] = 'Takasato'
        adata.obs['sc_protocol']   = '10x_3_v2'
        adata.obs['sequencing']    = 'Illumina_HiSeq_2500'
        adata.obs['genome_build']  = 'hg38'
        adata.obs['Age']           = 'd25'
        adata.obs['type']          = 'scRNAseq'
        adata.obs['disease']       = 'Control'
        adatas.append(adata)
        print(f'    -> Loaded {clean_sample_id}: {adata.shape}')

    adata = sc.concat(adatas, join='outer')
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1038_s42003_024_07069_6():
    """Load and harmonise d10_1038_s42003_024_07069_6."""
    # 15 - GSE263556 | d10_1038_s42003_024_07069_6
    # PPARα signalling / PT maturation - all 5 control timepoints
    # EXCLUDED: d23 control/treatment mixed matrix
    study    = 'd10_1038_s42003_024_07069_6'
    raw_base = os.path.join(DATA_ROOT, study, 'GSE263556_RAW')

    prefixes = {
        'd15': 'GSM8194405_d15_kidney_organoids_',
        'd19': 'GSM8194406_d19_kidney_organoids_',
        'd23': 'GSM8194407_d23_kidney_organoids_',
        'd27': 'GSM8194408_d27_kidney_organoids_',
        'd31': 'GSM8194409_d31_kidney_organoids_',
    }

    adatas = []
    for tp, prefix in prefixes.items():
        mtx_files = [f'{prefix}matrix.mtx.gz',
                     f'{prefix}barcodes.tsv.gz',
                     f'{prefix}features.tsv.gz']

        urls  = [geo_gse_url('GSE263556', fn) for fn in mtx_files]
        adata_sample = read_10x_mtx_with_redownload(raw_base, prefix, urls)

        adata_sample.obs_names_make_unique()
        adata_sample.var_names_make_unique()

        gsm_id = prefix.split('_')[0]

        adata_sample.obs['sample_id']     = gsm_id
        adata_sample.obs['Age']           = tp
        adata_sample.obs['study']         = study
        adata_sample.obs['source']        = 'hiPSC'
        adata_sample.obs['diff_protocol'] = 'Takasato_modified'
        adata_sample.obs['sc_protocol']   = '10x_3_v3.1'
        adata_sample.obs['sequencing']    = 'Illumina_NovaSeq_6000'
        adata_sample.obs['genome_build']  = 'hg38'
        adata_sample.obs['type']          = 'scRNAseq'
        adata_sample.obs['disease']       = 'Control'

        adatas.append(adata_sample)
        print(f'  Loaded {gsm_id} ({tp}): {adata_sample.shape}')

    import anndata as ad
    adata = ad.concat(adatas, join='outer', index_unique='-')
    adata.obs_names_make_unique()

    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1101_505396():
    """Load and harmonise d10_1101_505396."""
    # 16 - GSE138835 | d10_1101_505396
    # Uses genes.tsv.gz (v2 layout) instead of features.tsv.gz
    study    = 'd10_1101_505396'
    basepath = os.path.join(DATA_ROOT, study, 'GSE138835_RAW')
    prefix   = 'GSM4120180_'

    mtx_files_17 = [f'{prefix}matrix.mtx.gz',
                    f'{prefix}barcodes.tsv.gz',
                    f'{prefix}genes.tsv.gz']
    urls_17 = [geo_gsm_url('GSM4120180', fn) for fn in mtx_files_17]

    try:
        adata = read_10x_mtx_with_redownload(basepath, prefix, urls_17)
    except Exception:
        print('Standard MTX read failed - loading manually via scipy mmread')
        matrix_path   = os.path.join(basepath, f'{prefix}matrix.mtx.gz')
        features_path = os.path.join(basepath, f'{prefix}genes.tsv.gz')
        barcodes_path = os.path.join(basepath, f'{prefix}barcodes.tsv.gz')
        if not all(os.path.exists(p) for p in [matrix_path, features_path, barcodes_path]):
            redownload_urls(urls_17, basepath)
        with gzip.open(matrix_path, 'rb') as fh:
            X = mmread(fh).tocsr().T
        features    = pd.read_csv(features_path, sep='\t', header=None, compression='gzip')
        barcodes_df = pd.read_csv(barcodes_path, sep='\t', header=None, compression='gzip')
        gene_ids     = features[0].values
        gene_symbols = features[1].values if features.shape[1] >= 2 else features[0].values
        adata = sc.AnnData(X=X,
                           obs=pd.DataFrame(index=barcodes_df[0].values),
                           var=pd.DataFrame({'gene_ids': gene_ids}, index=gene_symbols))

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    adata.obs['sample_id']     = 'GSM4120180'
    adata.obs['study']         = study
    adata.obs['source']        = 'hiPSC'
    adata.obs['diff_protocol'] = 'Takasato'
    adata.obs['sc_protocol']   = '10x_3_v2'
    adata.obs['sequencing']    = 'Illumina_HiSeq_2500'
    adata.obs['genome_build']  = 'hg38'
    adata.obs['Age']           = 'd12'
    adata.obs['type']          = 'scRNAseq'
    adata.obs['disease']       = 'Control'
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1126_sciadv_abj5633():
    """Load and harmonise d10_1126_sciadv_abj5633."""
    # 17 - GSE165104 | d10_1126_sciadv_abj5633
    # 24 organoid samples: 8 timepoints × 3 replicates
    study    = 'd10_1126_sciadv_abj5633'
    basepath = os.path.join(DATA_ROOT, study, 'GSE165104', 'GSE165104_RAW')

    organoid_samples = {
        'GSM5026212_L417_83_':  {'sample_id': 'GSM5026212', 'Age': 'd0',  'rep': 1},
        'GSM5026213_L417_84_':  {'sample_id': 'GSM5026213', 'Age': 'd0',  'rep': 2},
        'GSM5026214_L417_85_':  {'sample_id': 'GSM5026214', 'Age': 'd0',  'rep': 3},
        'GSM5026215_L417_86_':  {'sample_id': 'GSM5026215', 'Age': 'd7',  'rep': 1},
        'GSM5026216_L417_87_':  {'sample_id': 'GSM5026216', 'Age': 'd7',  'rep': 2},
        'GSM5026217_L417_88_':  {'sample_id': 'GSM5026217', 'Age': 'd7',  'rep': 3},
        'GSM5026218_L417_91_':  {'sample_id': 'GSM5026218', 'Age': 'd10', 'rep': 1},
        'GSM5026219_L417_92_':  {'sample_id': 'GSM5026219', 'Age': 'd10', 'rep': 2},
        'GSM5026220_L417_93_':  {'sample_id': 'GSM5026220', 'Age': 'd10', 'rep': 3},
        'GSM5026221_L417_94_':  {'sample_id': 'GSM5026221', 'Age': 'd12', 'rep': 1},
        'GSM5026222_L417_95_':  {'sample_id': 'GSM5026222', 'Age': 'd12', 'rep': 2},
        'GSM5026223_L417_96_':  {'sample_id': 'GSM5026223', 'Age': 'd12', 'rep': 3},
        'GSM5026224_L417_97_':  {'sample_id': 'GSM5026224', 'Age': 'd14', 'rep': 1},
        'GSM5026225_L417_98_':  {'sample_id': 'GSM5026225', 'Age': 'd14', 'rep': 2},
        'GSM5026226_L417_99_':  {'sample_id': 'GSM5026226', 'Age': 'd14', 'rep': 3},
        'GSM5026227_L417_100_': {'sample_id': 'GSM5026227', 'Age': 'd24', 'rep': 1},
        'GSM5026228_L417_101_': {'sample_id': 'GSM5026228', 'Age': 'd24', 'rep': 2},
        'GSM5026229_L417_102_': {'sample_id': 'GSM5026229', 'Age': 'd24', 'rep': 3},
        'GSM5026233_L417_106_': {'sample_id': 'GSM5026233', 'Age': 'd26', 'rep': 1},
        'GSM5026234_L417_107_': {'sample_id': 'GSM5026234', 'Age': 'd26', 'rep': 2},
        'GSM5026235_L417_108_': {'sample_id': 'GSM5026235', 'Age': 'd26', 'rep': 3},
        'GSM5026239_L417_112_': {'sample_id': 'GSM5026239', 'Age': 'd28', 'rep': 1},
        'GSM5026240_L417_113_': {'sample_id': 'GSM5026240', 'Age': 'd28', 'rep': 2},
        'GSM5026241_L417_114_': {'sample_id': 'GSM5026241', 'Age': 'd28', 'rep': 3},
    }

    adatas = []
    for prefix, meta in organoid_samples.items():
        try:
            mtx_files = [f'{prefix}matrix.mtx.gz', f'{prefix}barcodes.tsv.gz',
                         f'{prefix}features.tsv.gz']
            urls  = [geo_gsm_url(meta['sample_id'], fn) for fn in mtx_files]
            tmp   = read_10x_mtx_with_redownload(basepath, prefix, urls)
            tmp.obs_names_make_unique()
            tmp.var_names_make_unique()
            tmp.obs['sample_id'] = meta['sample_id']
            tmp.obs['Age']       = meta['Age']
            adatas.append(tmp)
            print(f"  Loaded {meta['sample_id']} ({meta['Age']} rep{meta['rep']}): {tmp.shape}")
        except Exception as e:
            print(f"  SKIP {meta['sample_id']}: {e}")

    adata = sc.concat(adatas, join='outer')
    adata.obs['study']         = study
    adata.obs['source']        = 'iPSC'
    adata.obs['diff_protocol'] = 'Hybrid_Takasato_Morizane'
    adata.obs['sc_protocol']   = '10x_3_v3'
    adata.obs['sequencing']    = 'Illumina_HiSeq_2500'
    adata.obs['genome_build']  = 'hg38'
    adata.obs['type']          = 'scRNAseq'
    adata.obs['disease']       = 'Control'
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, 'GSE165104', f'{study}_GSE165104.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}_GSE165104: {adata.shape} ({len(adatas)} samples)')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_GSE165104_harmonized.h5ad')
    )

def load_d10_1172_jci_insight_122697():
    """Load and harmonise d10_1172_jci_insight_122697."""
    # 18 - GSE115986 | d10_1172_jci_insight_122697
    study     = 'd10_1172_jci_insight_122697'
    file_path = os.path.join(DATA_ROOT, study, 'GSE115986_RAW')

    adatas = []
    for names in os.listdir(file_path):
        if not names.endswith('.gz'):
            continue

        df    = pd.read_csv(os.path.join(file_path, names), sep='\t', index_col=0, compression='gzip')
        adata = sc.AnnData(df.T)

        gsm_id = names.split('_')[0]
        adata.obs['sample_id'] = gsm_id

        adata.obs_names_make_unique()
        adata.var_names_make_unique()

        adata.obs['study']         = study
        adata.obs['source']        = 'ESC'
        adata.obs['diff_protocol'] = 'Freedman'
        adata.obs['sc_protocol']   = 'Drop-seq'
        adata.obs['sequencing']    = 'Illumina_HiSeq_2000'
        adata.obs['genome_build']  = 'hg38'
        adata.obs['type']          = 'scRNAseq'
        adata.obs['disease']       = 'Control'
        adatas.append(adata)

    adata = sc.concat(adatas, join='outer')

    age_mapping = {
        'GSM3204315': 'd21',  # NLW33
        'GSM3204316': 'd19',  # NLW34a
        'GSM3204317': 'd19',  # NLW34b
        'GSM3204318': 'd19',  # NLW35a
        'GSM3204319': 'd19',  # NLW35b
        'GSM3204320': 'd19',  # NLW36a
        'GSM3204321': 'd19',  # NLW36b
        'GSM3204322': 'd20',  # NLW39
    }

    adata.obs['Age'] = adata.obs['sample_id'].map(age_mapping)

    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_d10_1242_dev_200198():
    """Load and harmonise d10_1242_dev_200198."""
    # 19 - GSE181954 | d10_1242_dev_200198
    study    = 'd10_1242_dev_200198'
    basepath = os.path.join(DATA_ROOT, study, 'GSE181954_RAW')

    prefix = 'GSM5514769_'

    adata = sc.read_10x_mtx(basepath, var_names='gene_symbols', prefix=prefix)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    adata.obs['sample_id']     = 'GSM5514769'  # Using the specific GSM ID for consistency
    adata.obs['study']         = study
    adata.obs['source']        = 'hiPSC'
    adata.obs['diff_protocol'] = 'Hybrid_Takasato_Uchimura'
    adata.obs['sc_protocol']   = '10x_3_v3.1'
    adata.obs['sequencing']    = 'Illumina_NovaSeq_6000'
    adata.obs['genome_build']  = 'hg38'
    adata.obs['Age']           = 'd24'
    adata.obs['type']          = 'scRNAseq'
    adata.obs['disease']       = 'Control'
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_dno_doi_kidney_organoid_GSE312980():
    """Load and harmonise dno_doi_kidney_organoid_GSE312980."""
    # 20 - GSE312980 | dno_doi_kidney_organoid_GSE312980
    study    = 'dno_doi_kidney_organoid_GSE312980'
    basepath = os.path.join(DATA_ROOT, study, 'GSE312980_RAW')

    samples = {
        'GSM9357919_Illumina-Kidney-CloneR2-3D_': {'sample_id': 'GSM9357919', 'condition': '3D'},
    }

    adatas = []
    for prefix, meta in samples.items():
        mtx_files = [f'{prefix}matrix.mtx.gz', f'{prefix}barcodes.tsv.gz',
                     f'{prefix}features.tsv.gz']
        urls  = [geo_gsm_url(meta['sample_id'], fn) for fn in mtx_files]
        adata = read_10x_mtx_with_redownload(basepath, prefix, urls)
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        adata.obs['sample_id'] = meta['sample_id']
        adata.obs['condition'] = meta['condition']
        adatas.append(adata)

    adata = sc.concat(adatas, join='outer')
    adata.obs['study']         = study
    adata.obs['source']        = 'hiPSC'
    adata.obs['diff_protocol'] = 'Freedman'
    adata.obs['sc_protocol']   = 'PIPseq_V_T10_3'
    adata.obs['sequencing']    = 'Illumina_NextSeq_2000'
    adata.obs['genome_build']  = 'hg38'
    adata.obs['Age']           = 'd18'
    adata.obs['type']          = 'scRNAseq'
    adata.obs['disease']       = 'Control'
    fill_missing_meta(adata)

    raw_out = os.path.join(DATA_ROOT, study, f'{study}.h5ad')
    os.makedirs(os.path.dirname(raw_out), exist_ok=True)
    adata.write(raw_out)
    print(f'Saved {study}: {adata.shape}')

    memory_efficient_harmonize_sparse(
        input_file=raw_out,
        output_file=os.path.join(HARMONIZED_DIR, f'{study}_harmonized.h5ad')
    )

def load_dno_doi_kidney_organoid_LMRH():
    """Load and harmonise dno_doi_kidney_organoid_LMRH."""
    lmrh_cfg = CFG["lmrh"]
    study = lmrh_cfg["study_id"]
    lmrh_dir = os.path.join(DATA_ROOT, study, lmrh_cfg["input_subdir"])

    adatas_lmrh = []
    for filename in sorted(os.listdir(lmrh_dir)):
        if not filename.endswith(".h5ad"):
            continue
        file_path = os.path.join(lmrh_dir, filename)
        print(f"  Loading LMRH file {filename}...")
        adata_sample = sc.read_h5ad(file_path)
        adata_sample.obs_names_make_unique()
        adata_sample.var_names_make_unique()
        if "sample_id" not in adata_sample.obs.columns:
            prefix = filename.replace(".h5ad", "")
            adata_sample.obs["sample_id"] = prefix.split("_")[0]
        adata_sample.obs["study"] = study
        adata_sample.obs["source"] = lmrh_cfg["source"]
        adata_sample.obs["diff_protocol"] = lmrh_cfg["diff_protocol"]
        adata_sample.obs["sc_protocol"] = lmrh_cfg["sc_protocol"]
        adata_sample.obs["sequencing"] = lmrh_cfg["sequencing"]
        adata_sample.obs["genome_build"] = CFG["metadata"]["default_genome_build"]
        adata_sample.obs["type"] = CFG["metadata"]["default_type"]
        adata_sample.obs["disease"] = CFG["metadata"]["default_disease"]
        adatas_lmrh.append(adata_sample)
        print(f"    -> {adata_sample.shape}")

    if not adatas_lmrh:
        raise RuntimeError(f"No LMRH h5ad files found in {lmrh_dir}")

    adata_lmrh = sc.concat(adatas_lmrh, join="outer")
    fill_missing_meta(adata_lmrh)

    raw_lmrh = os.path.join(DATA_ROOT, study, f"{study}.h5ad")
    maybe_write_h5ad(adata_lmrh, raw_lmrh, label="raw LMRH")

    harmonized_lmrh_only = os.path.join(
        HARMONIZED_DIR, f"{study}{lmrh_cfg['only_harmonized_suffix']}"
    )
    memory_efficient_harmonize_sparse(input_file=raw_lmrh, output_file=harmonized_lmrh_only)

    harmonized_clones = []
    if CFG["lmrh"].get("include_clone_controls", "true").lower() in {"1", "true", "yes"}:
        for section in ("clone_g7", "clone_ttc21b"):
            spec = CFG[section]
            mtx_path = os.path.join(DATA_ROOT, spec["study_id"], spec["mtx_subdir"])
            print(f"Loading {spec['study_id']} from {mtx_path}...")
            adata = sc.read_10x_mtx(mtx_path, var_names="gene_symbols", cache=False)
            adata.obs_names_make_unique()
            adata.var_names_make_unique()
            adata.obs["sample_id"] = spec["sample_id"]
            adata.obs["study"] = spec["study_id"]
            adata.obs["source"] = spec["source"]
            adata.obs["diff_protocol"] = spec["diff_protocol"]
            adata.obs["sc_protocol"] = spec["sc_protocol"]
            adata.obs["sequencing"] = spec["sequencing"]
            adata.obs["genome_build"] = CFG["metadata"]["default_genome_build"]
            adata.obs["Age"] = spec["Age"]
            adata.obs["type"] = CFG["metadata"]["default_type"]
            adata.obs["disease"] = CFG["metadata"]["default_disease"]
            adata.obs["condition"] = spec["condition"]
            fill_missing_meta(adata)

            raw_out = os.path.join(DATA_ROOT, spec["study_id"], f"{spec['study_id']}.h5ad")
            maybe_write_h5ad(adata, raw_out, label=f"raw {spec['study_id']}")
            harm_out = os.path.join(HARMONIZED_DIR, f"{spec['study_id']}_harmonized.h5ad")
            memory_efficient_harmonize_sparse(input_file=raw_out, output_file=harm_out)
            harmonized_clones.append(harm_out)

    print("Building freeze LMRH object (LMRH + G7 + TTC21B)...")
    adata_lmrh_h = sc.read_h5ad(harmonized_lmrh_only)
    parts = [adata_lmrh_h] + [sc.read_h5ad(p) for p in harmonized_clones]
    for a in parts:
        a.obs_names_make_unique()
    adata_freeze = sc.concat(parts, join="outer", merge="same", index_unique="-")
    print(f"Freeze LMRH shape: {adata_freeze.shape}")
    print(adata_freeze.obs["study"].value_counts())
    print(adata_freeze.obs["sample_id"].value_counts())

    freeze_out = os.path.join(HARMONIZED_DIR, f"{study}{lmrh_cfg['freeze_harmonized_suffix']}")
    maybe_write_h5ad(adata_freeze, freeze_out, label="freeze LMRH")

STUDY_LOADERS = {
    'd10_1016_j_cell_2020_04_004': load_d10_1016_j_cell_2020_04_004,
    'd10_1016_j_celrep_2020_108514': load_d10_1016_j_celrep_2020_108514,
    'd10_1016_j_cmet_2022_04_009': load_d10_1016_j_cmet_2022_04_009,
    'd10_1016_j_devcel_2019_06_001': load_d10_1016_j_devcel_2019_06_001,
    'd10_1016_j_stem_2018_04_022': load_d10_1016_j_stem_2018_04_022,
    'd10_1016_j_stem_2018_10_010': load_d10_1016_j_stem_2018_10_010,
    'd10_1016_j_stem_2019_06_009': load_d10_1016_j_stem_2019_06_009,
    'd10_1016_j_stem_2021_11_001': load_d10_1016_j_stem_2021_11_001,
    'd10_1016_j_stem_2022_06_005': load_d10_1016_j_stem_2022_06_005,
    'd10_1038_s41467_019_13382_0': load_d10_1038_s41467_019_13382_0,
    'd10_1038_s41467_022_33623_z': load_d10_1038_s41467_022_33623_z,
    'd10_1038_s41551_025_01542_1': load_d10_1038_s41551_025_01542_1,
    'd10_1038_s41563_020_00853_9': load_d10_1038_s41563_020_00853_9,
    'd10_1038_s41592_018_0253_2': load_d10_1038_s41592_018_0253_2,
    'd10_1038_s42003_024_07069_6': load_d10_1038_s42003_024_07069_6,
    'd10_1101_505396': load_d10_1101_505396,
    'd10_1126_sciadv_abj5633': load_d10_1126_sciadv_abj5633,
    'd10_1172_jci_insight_122697': load_d10_1172_jci_insight_122697,
    'd10_1242_dev_200198': load_d10_1242_dev_200198,
    'dno_doi_kidney_organoid_GSE312980': load_dno_doi_kidney_organoid_GSE312980,
    'dno_doi_kidney_organoid_LMRH': load_dno_doi_kidney_organoid_LMRH,
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Harmonise 21 atlas datasets (config-driven loaders)."
    )
    p.add_argument(
        "--config",
        default=DEFAULT_CONFIG,
        help=f"Path to INI config (default: {DEFAULT_CONFIG})",
    )
    p.add_argument(
        "--studies",
        default=None,
        help="Comma-separated study keys to run (default: [studies] enabled/all)",
    )
    p.add_argument(
        "--list-studies",
        action="store_true",
        help="Print configured study keys and exit",
    )
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    load_config(args.config)

    print("Node:", platform.node())
    print("CPUs:", os.cpu_count())
    print(f"DATA_ROOT      : {DATA_ROOT}")
    print(f"OUTPUT_ROOT    : {OUTPUT_ROOT}")
    print(f"HARMONIZED_DIR : {HARMONIZED_DIR}")
    print(f"GTF_FILE       : {GTF_FILE}")

    if args.list_studies:
        for s in _split_csv(CFG["studies"]["list"]):
            print(s)
        return 0

    build_allowed_genes()
    studies = enabled_studies(args.studies)
    print(f"Running {len(studies)} studies...")
    for study in studies:
        print("\n" + "=" * 70)
        print(study)
        print("=" * 70)
        STUDY_LOADERS[study]()
    print("\nDone. Harmonised objects are under:", HARMONIZED_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
