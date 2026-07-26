#!/usr/bin/env python3
"""
scRNA-seq Atlas Harmonization Script
=======================================
Reads a metadata CSV, loads each sample (H5 / H5AD / MTX),
concatenates samples per study, then harmonizes gene space to
the 24,100 protein-coding + lncRNA genes from GRCh38.104.

Usage
-----
    hkoca qc-filter [-w DIR] [OPTIONS]
    python -m hkoca.qc_filter [-w DIR] [OPTIONS]

Config files: harmonize.config and qc_config.dcf (CWD or package).
Precedence: CLI flags > config files > environment variables.

All arguments are optional when the corresponding paths are set in
harmonize.config or as environment variables.

Options
-------
    -w, --working-dir DIR
                     Base directory used to resolve all relative paths in
                     the config file and the metadata CSV.
                     Defaults to the current working directory.
    --config PATH    Path to a harmonize.config INI file.
                     Auto-discovered as harmonize.config next to the script,
                     or in the current working directory.
    --csv PATH       Path to the metadata CSV describing all samples.
                     Overrides [paths] metadata_csv and METADATA_CSV env var.
    --gtf PATH       Path to the Ensembl GRCh38.104 GTF (plain or .gz).
                     Overrides [paths] gtf_file and GTF_FILE env var.
    --output PATH    Root output directory. Per-study sub-folders (raw/,
                     harmonized/, rds/) are created inside it automatically.
                     Overrides [paths] output_root and OUTPUT_ROOT env var.
    --transgenes LIST
                     Comma-separated gene names to preserve in the reference
                     gene set (e.g. --transgenes EGFP,mCherry,Cre).
                     Use the exact gene name as it appears in your count matrix.
                     Can also be set via [transgenes] names = ... in the
                     config file.
    --to-rds         After harmonizing, convert each study's h5ad to a
                     native Seurat .rds object via rpy2 + anndata2ri.
                     Requires R, Seurat, SingleCellExperiment, rpy2,
                     and anndata2ri.
    --summary        Run the full pipeline, then generate summary plots
                     and dataset_summary.csv in <output>/<report_subdir>/.
    --summary-only   Skip the pipeline entirely and only (re-)generate
                     plots from *_harmonized.h5ad files already under
                     <output>.

Configuration precedence (highest -> lowest)
---------------------------------------------
    1. CLI flags
    2. harmonize.config  ([paths], [transgenes], [summary] sections)
    3. Environment variables (GTF_FILE, METADATA_CSV, OUTPUT_ROOT, WORKING_DIR)

Config file sections
--------------------
    [paths]
        working_dir   = ./
        gtf_file      = path/to/Homo_sapiens.GRCh38.104.gtf
        metadata_csv  = path/to/datasets_metadata.csv
        output_root   = path/to/results

    [transgenes]
        names = EGFP,MCHERRY,TTC21B   # preserved even if absent from GTF

    [summary]
        figure_dpi          = 300     # must be > 0
        figure_extensions   = png, pdf
        report_subdir       = reports/atlas_summary
        age_plot_top_n      = 15      # must be > 0

Outputs per study
-----------------
    <output>/<study>/raw/<study>.h5ad
    <output>/<study>/harmonized/<study>_harmonized.h5ad
    <output>/<study>/rds/<study>_harmonized.rds   (only with --to-rds)
    <working_dir>/harmonize_pipeline.log

Examples
--------
    # All paths from harmonize.config in CWD
    hkoca qc-filter

    # Override output; add custom transgenes
    hkoca qc-filter --output /data/results --transgenes EGFP,MCHERRY

    # Full run: harmonize (--to-rds) then QC
    hkoca qc-filter --csv config/meta.csv --gtf genes.gtf --output results --summary
"""

import os
import sys
import gzip
import collections
import configparser
import warnings
import logging

import pandas as pd
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import anndata
import h5py
from scipy.io import mmread
from scipy.sparse import csr_matrix

from hkoca.qc_filter.harmonize.config import get_summary_options

warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)
anndata.settings.allow_write_nullable_strings = True

# Initialize global logger
logger = logging.getLogger("harmonize_atlas")

def setup_logging(working_dir: str):
    """Set up dual logging: INFO to console, DEBUG to file."""
    if logger.hasHandlers():
        logger.handlers.clear()

    logger.setLevel(logging.DEBUG)
    
    # 1. Console Handler (Clean, INFO level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch_fmt = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%H:%M:%S')
    ch.setFormatter(ch_fmt)
    logger.addHandler(ch)

    # 2. File Handler (Detailed, DEBUG level)
    log_file = os.path.join(working_dir, "harmonize_pipeline.log")
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.DEBUG)
    fh_fmt = logging.Formatter('%(asctime)s [%(levelname)s] %(funcName)s:%(lineno)d - %(message)s')
    fh.setFormatter(fh_fmt)
    logger.addHandler(fh)
    
    logger.info(f"Detailed debug log initialized at: {log_file}")


# ─────────────────────────────────────────────────────────────────────────────
# GTF PARSING
# ─────────────────────────────────────────────────────────────────────────────

def load_allowed_genes(gtf_path: str) -> set:
    if not os.path.exists(gtf_path):
        logger.error(f"GTF file not found: {gtf_path}")
        raise FileNotFoundError(f"GTF file not found: {gtf_path!r}")

    open_fn = gzip.open if gtf_path.endswith(".gz") else open
    mode    = "rt"     if gtf_path.endswith(".gz") else "r"
    genes   = []

    logger.info(f"Parsing GTF: {os.path.basename(gtf_path)} ...")
    with open_fn(gtf_path, mode) as fh:
        for line in fh:
            if line.startswith("#"): continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene": continue
            
            info = {}
            # FIXED: Split by ';' and strip whitespace to handle Ensembl formatting inconsistencies
            for kv in cols[8].split(";"):
                kv = kv.strip()
                if not kv: continue
                parts = kv.split(" ", 1)
                if len(parts) == 2:
                    info[parts[0]] = parts[1].strip('";')
                    
            if info.get("gene_biotype") in ("protein_coding", "lncRNA"):
                g = info.get("gene_name")
                if g and pd.notna(g):
                    genes.append(g)

    allowed = set(genes)
    if not allowed:
        logger.error("allowed_genes is empty after parsing the GTF.")
        raise RuntimeError("allowed_genes is empty. Check GTF gene_biotype annotations.")
    
    logger.info(f"GTF parsed: {len(allowed):,} protein_coding + lncRNA genes")
    return allowed


# ─────────────────────────────────────────────────────────────────────────────
# FILE READERS
# ─────────────────────────────────────────────────────────────────────────────

def _decode(arr) -> list:
    return [x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x) for x in arr]

def read_h5_safe(path: str) -> sc.AnnData:
    try:
        return sc.read_10x_h5(path)
    except Exception as e1:
        logger.warning(f"Standard H5 read failed ({e1}). Attempting manual fallback...")
        try:
            with h5py.File(path, "r") as f:
                g = f["matrix"] if "matrix" in f else f[list(f.keys())[0]]
                barcodes = _decode(g["barcodes"][()])
                gene_names = None
                for loc, key in [("features", "name"), ("features", "gene_names"),
                                  ("features", "id"), (None, "gene_names"), (None, "genes")]:
                    try:
                        src = g[loc] if loc else g
                        gene_names = _decode(src[key][()]); break
                    except Exception: continue
                if gene_names is None:
                    raise RuntimeError("Gene names not found in H5")
                X = csr_matrix((g["data"][()], g["indices"][()], g["indptr"][()]),
                               shape=tuple(g["shape"][()])).T
            return sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=gene_names))
        except Exception as e2:
            logger.exception(f"Cannot read H5 fallback for {path}")
            raise RuntimeError(f"Cannot read H5: {path}\n  err1={e1}\n  err2={e2}")

def read_mtx_safe(folder: str, prefix: str = "") -> sc.AnnData:
    try:
        return sc.read_10x_mtx(folder, var_names="gene_symbols", prefix=prefix)
    except Exception as e:
        logger.warning(f"sc.read_10x_mtx failed ({e}). Attempting manual fallback load...")
        mtx_path      = os.path.join(folder, f"{prefix}matrix.mtx.gz")
        barcodes_path = os.path.join(folder, f"{prefix}barcodes.tsv.gz")
        features_path = None
        for fname in [f"{prefix}features.tsv.gz", f"{prefix}genes.tsv.gz"]:
            if os.path.exists(os.path.join(folder, fname)):
                features_path = os.path.join(folder, fname)
                break
        if not features_path:
            raise FileNotFoundError(f"No features/genes file for prefix={prefix!r}")
        with gzip.open(mtx_path, "rb") as fh:
            X = mmread(fh).tocsr().T
        feats    = pd.read_csv(features_path, sep="\t", header=None, compression="gzip")
        barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, compression="gzip")
        gene_names = feats[1].values if feats.shape[1] > 1 else feats[0].values
        return sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes[0].values), var=pd.DataFrame(index=gene_names))

def geo_redownload(folder: str, prefix: str, sample_id: str, feature_file: str = "features.tsv.gz") -> None:
    import urllib.request
    if len(sample_id) < 9:
        raise ValueError(f"GEO re-download requires sample_id of length >= 9, got {sample_id!r}")
    gsm_prefix = f"{sample_id[:-3]}nnn"
    base_url   = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm_prefix}/{sample_id}/suppl"
    for fname in [f"{prefix}matrix.mtx.gz", f"{prefix}barcodes.tsv.gz", f"{prefix}{feature_file}"]:
        url   = f"{base_url}/{fname}"
        fpath = os.path.join(folder, fname)
        if os.path.exists(fpath):
            os.replace(fpath, fpath + ".corrupt")
        logger.info(f"Redownloading {url}...")
        try:
            urllib.request.urlretrieve(url, fpath)
        except Exception as e:
            logger.error(f"GEO re-download failed for {fname}: {e}")
            raise RuntimeError(f"GEO re-download failed for {fname}: {e}") from e


# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE LOADER
# ─────────────────────────────────────────────────────────────────────────────

# Columns that every metadata CSV must contain.
MANDATORY_META_COLS = [
    "data_path", "sample_id", "study", "source",
    "diff_protocol", "sc_protocol", "sequencing", "genome_build", "Age", "type",
]

# Optional pipeline-control columns: recognised when present but never written to obs.
#   file_prefix : MTX folder prefix (e.g. "sample1_"); defaults to "" when absent.
#   output_dir  : per-study output path, only used when --output is not supplied on the CLI.
#   skip        : set to "True" to exclude a row from the run.
_OPTIONAL_PIPELINE_COLS = {"file_prefix", "output_dir", "skip"}

# All columns that must never be written to AnnData obs.
_NON_OBS_COLS = {"data_path"} | _OPTIONAL_PIPELINE_COLS

def _get_obs_cols(df: pd.DataFrame) -> list:
    """Return all columns from *df* that should be written to AnnData .obs.

    Mandatory columns that carry biological/metadata meaning are included first
    (in their declared order), followed by any additional optional columns the
    user has added, excluding pipeline-control columns (data_path, file_prefix,
    output_dir, skip).
    """
    mandatory_obs = [c for c in MANDATORY_META_COLS if c not in _NON_OBS_COLS]
    optional_obs  = [c for c in df.columns
                     if c not in MANDATORY_META_COLS and c not in _NON_OBS_COLS]
    # Preserve order: mandatory first, then optional, deduplicated
    seen = set()
    ordered = []
    for c in mandatory_obs + optional_obs:
        if c not in seen and c in df.columns:
            seen.add(c)
            ordered.append(c)
    return ordered

def load_sample(row: pd.Series, working_dir: str | None = None,
                obs_cols: list | None = None) -> sc.AnnData:
    path = str(row["data_path"]).strip()
    if not path:
        raise ValueError("data_path is empty")
    if working_dir and not os.path.isabs(path):
        path = os.path.normpath(os.path.join(working_dir, path))
    prefix    = str(row.get("file_prefix", "") or "").strip()
    if prefix.lower() in ("nan", "none"): prefix = ""
    sample_id = str(row.get("sample_id", "") or "").strip()

    if not os.path.exists(path):
        raise FileNotFoundError(f"data_path not found: {path}")

    if os.path.isfile(path):
        ext = path.lower()
        if ext.endswith(".h5ad"):
            logger.debug(f"[{sample_id}] Reading h5ad: {os.path.basename(path)}")
            adata = sc.read_h5ad(path)
        elif ext.endswith(".h5") or ext.endswith(".hdf5"):
            logger.debug(f"[{sample_id}] Reading H5: {os.path.basename(path)}")
            adata = read_h5_safe(path)
        else:
            raise ValueError(f"Unknown file extension: {path}")
    elif os.path.isdir(path):
        logger.debug(f"[{sample_id}] Reading MTX folder: {os.path.basename(path)} (prefix={prefix!r})")
        try:
            adata = read_mtx_safe(path, prefix=prefix)
        except Exception as e1:
            if sample_id.startswith("GSM"):
                logger.warning(f"MTX read failed for {sample_id}. Attempting GEO re-download...")
                try:
                    geo_redownload(path, prefix, sample_id)
                except Exception as e2:
                    raise RuntimeError(f"MTX read failed: {e1}; GEO re-download failed: {e2}") from e2
                adata = read_mtx_safe(path, prefix=prefix)
            else:
                raise
    else:
        raise FileNotFoundError(f"data_path not found: {path}")

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    cols_to_write = obs_cols if obs_cols is not None else [
        c for c in MANDATORY_META_COLS if c not in _NON_OBS_COLS
    ]
    for col in cols_to_write:
        val = row.get(col, "")
        if pd.notna(val) and str(val).strip():
            adata.obs[col] = str(val).strip()
        else:
            adata.obs[col] = "Unknown"
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# HARMONIZATION
# ─────────────────────────────────────────────────────────────────────────────

def harmonize_matrix_sparse(adata: sc.AnnData, allowed_genes: set) -> sc.AnnData:
    """
    Harmonize an AnnData to the allowed_genes reference set using native sparse ops.
    - Genes in both dataset and allowed_genes → counts preserved.
    - Genes in allowed_genes but not in dataset → zero-filled.
    - Genes in dataset but not in allowed_genes → dropped.
    Output matrix is float32 CSR.

    To preserve transgene counts, pass their names as part of allowed_genes
    (via --transgenes on the CLI or the [transgenes] section in harmonize.config).
    """
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)

    target_genes = sorted(allowed_genes)
    common_genes = sorted(set(adata.var_names.tolist()) & allowed_genes)

    gene_to_idx = {g: i for i, g in enumerate(target_genes)}
    sub         = adata[:, common_genes].copy()
    coo         = sp.coo_matrix(sub.X)
    new_cols    = np.array([gene_to_idx[g] for g in common_genes])[coo.col]
    new_X       = sp.coo_matrix(
        (coo.data, (coo.row, new_cols)),
        shape=(adata.shape[0], len(target_genes))
    ).tocsr().astype(np.float32)

    return sc.AnnData(
        X=new_X,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=target_genes),
    )


# ─────────────────────────────────────────────────────────────────────────────
# CSV LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_metadata_csv(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        logger.error(f"CSV not found: {csv_path}")
        raise FileNotFoundError(f"CSV not found: {csv_path!r}")

    with open(csv_path, encoding="utf-8") as fh:
        header = fh.readline()
    sep = ";" if header.count(";") > header.count(",") else ","
    if sep == ";":
        logger.warning("CSV uses semicolons — reading with sep=';'. Re-save as comma-separated to avoid issues.")

    df = pd.read_csv(csv_path, sep=sep, keep_default_na=False, encoding="utf-8")
    df.columns = df.columns.str.strip()

    required = list(MANDATORY_META_COLS)
    missing  = [c for c in required if c not in df.columns]
    if missing:
        logger.error(f"CSV missing required columns: {missing}")
        raise ValueError(f"CSV missing required columns: {missing}")

    # Warn (not error) when optional pipeline-control columns are absent so the
    # user knows which defaults will be applied.
    if "file_prefix" not in df.columns:
        logger.info("Column 'file_prefix' not found — defaulting to no prefix for all samples.")
    if "output_dir" not in df.columns:
        logger.info(
            "Column 'output_dir' not found — per-study output paths will be derived from "
            "--output / output_root config. Pass --output or add the column if needed."
        )

    optional_cols = [c for c in df.columns
                     if c not in MANDATORY_META_COLS and c not in _NON_OBS_COLS]
    if optional_cols:
        logger.info(f"Optional metadata columns detected and will be propagated to obs: {optional_cols}")

    if "skip" in df.columns:
        n_skip = (df["skip"].astype(str).str.strip().str.lower() == "true").sum()
        if n_skip:
            logger.info(f"Skipping {n_skip} row(s) with skip=True")
        df = df[df["skip"].astype(str).str.strip().str.lower() != "true"]

    if df.empty:
        logger.error("CSV has no rows left after skip filter.")
        raise ValueError("CSV has no rows left after skip filter.")

    empty_path = df["data_path"].astype(str).str.strip() == ""
    if empty_path.any():
        rows = (df.index[empty_path] + 2).tolist()
        logger.error(f"CSV has empty data_path in row(s) {rows}")
        raise ValueError(f"CSV has empty data_path in row(s) {rows}")

    empty_study = df["study"].astype(str).str.strip() == ""
    if empty_study.any():
        rows = (df.index[empty_study] + 2).tolist()
        logger.error(f"CSV has empty study in row(s) {rows}")
        raise ValueError(f"CSV has empty study in row(s) {rows}")

    id_series = df["sample_id"].astype(str).str.strip()
    id_dups = id_series[id_series.duplicated(keep=False)]
    if not id_dups.empty:
        dup_ids = sorted(id_dups.unique().tolist())
        logger.error(f"CSV has duplicate sample_id: {dup_ids}")
        raise ValueError(f"CSV has duplicate sample_id: {dup_ids}")

    return df

# ─────────────────────────────────────────────────────────────────────────────
# RDS CONVERSION (R/Seurat)
# ─────────────────────────────────────────────────────────────────────────────

def convert_h5ad_to_rds(h5ad_path: str, rds_path: str) -> None:
    """Converts a saved .h5ad file directly to a Seurat .rds object using rpy2."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import default_converter, pandas2ri
        from rpy2.robjects.conversion import localconverter
        import anndata2ri
        import gc

        # rpy2 >= 3.5.x removed pandas2ri.activate(); use a localconverter
        # context that combines the default, pandas, and anndata2ri converters.
        py2r_converter = default_converter + pandas2ri.converter + anndata2ri.converter
        logger.info(f"Starting Seurat conversion for {os.path.basename(h5ad_path)}...")
        
        adata = sc.read_h5ad(h5ad_path)
        if adata.n_obs == 0:
            logger.warning("0 cells found. Skipping RDS conversion.")
            return

        adata.obs.index = adata.obs.index.astype(str)
        adata.var.index = adata.var.index.astype(str)
        if not adata.obs.index.is_unique: adata.obs_names_make_unique()
        if not adata.var.index.is_unique: adata.var_names_make_unique()

        if sp.issparse(adata.X):
            X_clean = adata.X.tocsr().astype('float64')
        else:
            X_clean = sp.csr_matrix(adata.X).astype('float64')

        if np.isnan(X_clean.data).any():
            logger.warning("NAs found in count matrix — replacing with 0")
            X_clean.data = np.nan_to_num(X_clean.data, nan=0.0)

        obs_clean = pd.DataFrame(index=adata.obs.index)
        for col in adata.obs.columns:
            obs_clean[col] = adata.obs[col].astype(str)

        var_clean = pd.DataFrame(index=adata.var.index)
        for col in adata.var.columns:
            var_clean[col] = adata.var[col].astype(str)

        clean_adata = sc.AnnData(X=X_clean, obs=obs_clean, var=var_clean)

        # Guard: rds_path is interpolated directly into an R f-string.
        # Reject paths that contain characters which would break R string
        # literals or allow code injection (quotes, newlines).
        _unsafe = set('"\'\'\n\r')
        if any(c in rds_path for c in _unsafe):
            raise ValueError(
                f"rds_path contains characters unsafe for R interpolation: {rds_path!r}"
            )

        with localconverter(py2r_converter):
            ro.globalenv["adata_sce"] = clean_adata

        ro.r(f'''
            suppressPackageStartupMessages(library(Seurat))
            suppressPackageStartupMessages(library(SingleCellExperiment))
            suppressPackageStartupMessages(library(Matrix))

            counts_mat <- assay(adata_sce, assayNames(adata_sce)[1])

            # CreateSeuratObject already computes nCount_RNA and nFeature_RNA;
            # we pass meta.data directly so orig.ident and all other columns
            # are present from the start — no slot overwrite needed.
            cell_meta  <- as.data.frame(colData(adata_sce))
            seurat_obj <- CreateSeuratObject(
                counts    = counts_mat,
                assay     = "RNA",
                meta.data = cell_meta
            )

            saveRDS(seurat_obj, file = "{rds_path}")
        ''')
        
        logger.info(f"Successfully saved RDS: {rds_path}")

        del adata, clean_adata, X_clean, obs_clean, var_clean
        ro.r('rm(list = c("adata_sce", "seurat_obj", "counts_mat")); gc()')
        gc.collect()

    except ImportError:
        logger.error("rpy2 or anndata2ri not installed. Cannot convert to RDS.")
        raise
    except Exception as e:
        logger.exception("RDS conversion failed.")
        raise RuntimeError(f"RDS conversion failed: {e}")

# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(metadata_csv: str, gtf_file: str, output_root: str,
                 working_dir: str, to_rds: bool,
                 transgene_names: set | None = None) -> list:
    allowed_genes = load_allowed_genes(gtf_file)

    # Merge user-supplied transgene names into the reference gene set
    if transgene_names:
        allowed_genes = allowed_genes | set(transgene_names)
        logger.info(
            f"Reference gene set extended with {len(transgene_names)} transgene(s): "
            f"{sorted(transgene_names)}"
        )

    df_meta = load_metadata_csv(metadata_csv)
    logger.info(f"CSV loaded: {len(df_meta)} sample rows, {df_meta['study'].nunique()} unique studies")

    obs_cols = _get_obs_cols(df_meta)
    logger.debug(f"Columns propagated to AnnData obs: {obs_cols}")

    study_groups = collections.OrderedDict()
    for _, row in df_meta.iterrows():
        study_groups.setdefault(row["study"], []).append(row)

    failed = []

    for study, rows in study_groups.items():
        logger.info("-" * 64)
        logger.info(f"Processing study: {study} ({len(rows)} sample(s))")

        if output_root:
            out_dir = os.path.join(output_root, study)
        else:
            out_dir = str(rows[0].get("output_dir", study)).strip() or study
            if out_dir and not os.path.isabs(out_dir):
                out_dir = os.path.normpath(os.path.join(working_dir, out_dir))
        
        out_dir_raw  = os.path.join(out_dir, "raw")
        out_dir_harm = os.path.join(out_dir, "harmonized")
        os.makedirs(out_dir_raw,  exist_ok=True)
        os.makedirs(out_dir_harm, exist_ok=True)
        
        raw_h5ad        = os.path.join(out_dir_raw,  f"{study}.h5ad")
        harmonized_h5ad = os.path.join(out_dir_harm, f"{study}_harmonized.h5ad")
        rds_file: str | None = None  # assigned below only when to_rds=True

        if to_rds:
            out_dir_rds = os.path.join(out_dir, "rds")
            os.makedirs(out_dir_rds, exist_ok=True)
            rds_file = os.path.join(out_dir_rds, f"{study}_harmonized.rds")

        # ── Load all samples
        adatas = []
        for row in rows:
            sample_id = row.get('sample_id', '?')
            try:
                adata = load_sample(row, working_dir=working_dir, obs_cols=obs_cols)
                adatas.append(adata)
                logger.info(f"Loaded sample {sample_id} → {adata.shape}")
            except Exception as e:
                logger.error(f"FAILED to load {sample_id}: {e}")
                logger.debug("Exception traceback:", exc_info=True)
                failed.append({"study": study, "sample_id": sample_id, "error": str(e)})

        if not adatas:
            logger.warning(f"No valid samples loaded for {study}. Skipping.")
            continue

        # ── Concatenate
        if len(adatas) == 1:
            adata = adatas[0]
        else:
            adata = sc.concat(adatas, join="outer")
            adata.obs_names_make_unique()
            logger.info(f"Concatenated {len(adatas)} samples → {adata.shape}")
        del adatas  # free per-sample objects; adata now holds the study-level data

        for col in adata.obs.select_dtypes(include=["string", "object"]).columns:
            adata.obs[col] = adata.obs[col].astype("object")

        # ── Save raw
        try:
            adata.write(raw_h5ad)
            logger.info(f"Saved raw h5ad: {raw_h5ad}")
        except Exception as e:
            logger.error(f"Failed to write raw h5ad for {study}: {e}")
            logger.debug("Exception traceback:", exc_info=True)

        # ── Harmonize (passed directly from memory)
        try:
            logger.info(f"Harmonizing {study} to standard gene space...")

            adata_new = harmonize_matrix_sparse(adata, allowed_genes)
            logger.info(f"Harmonized shape: {adata_new.shape}")

            adata_new.write(harmonized_h5ad)
            logger.info(f"Saved harmonized h5ad: {harmonized_h5ad}")

            if to_rds:
                convert_h5ad_to_rds(harmonized_h5ad, rds_file)
                
        except Exception as e:
            logger.error(f"Harmonization/Conversion failed for {study}: {e}")
            logger.debug("Exception traceback:", exc_info=True)
            failed.append({"study": study, "sample_id": "harmonize_or_rds", "error": str(e)})

    logger.info("=" * 64)
    if failed:
        logger.warning(f"PIPELINE DONE with {len(failed)} failure(s).")
        for f in failed:
            logger.warning(f"  [{f['study']} / {f['sample_id']}] {f['error']}")
    else:
        logger.info("PIPELINE FULLY COMPLETE — No failures.")

    return failed


# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY & PLOTS
# ─────────────────────────────────────────────────────────────────────────────

def _get_obs_meta(obs: pd.DataFrame, col: str,
                  warned_missing: set | None = None) -> str:
    """Return a representative string value from one obs column of a backed h5ad.

    Returns ``"Unknown"`` when the column is absent or contains only NaN/nan.
    Returns the single unique value when all cells agree, otherwise a
    comma-separated sorted list of unique values.

    When *warned_missing* is supplied (a set shared across calls for the same
    file), a warning is emitted the first time a column is found to be absent,
    so the user knows which expected fields are missing from their data.
    """
    if col not in obs.columns:
        if warned_missing is not None and col not in warned_missing:
            logger.warning(
                f"Summary: expected obs column '{col}' not found in this dataset — "
                "it will appear as 'Unknown' in the report. "
                "Add it to your metadata CSV as an optional column if needed."
            )
            warned_missing.add(col)
        return "Unknown"
    vals = [
        str(v) for v in obs[col].unique()
        if pd.notna(v) and str(v).lower() != "nan"
    ]
    if not vals:
        return "Unknown"
    return vals[0] if len(vals) == 1 else ", ".join(sorted(vals))


def run_summary(scan_path: str, cfg: configparser.ConfigParser | None = None) -> None:
    import glob
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    opts = get_summary_options(cfg or configparser.ConfigParser())
    plt.style.use("default")
    sns.set_theme(
        style="whitegrid",
        context="paper",
        font_scale=1.35,
        rc={
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white",
        },
    )
    plt.rcParams["figure.dpi"]    = opts["figure_dpi"]
    plt.rcParams["savefig.dpi"]   = opts["figure_dpi"]
    plt.rcParams["savefig.bbox"]  = "tight"
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["savefig.edgecolor"] = "white"

    report_dir = os.path.join(scan_path, opts["report_subdir"])
    os.makedirs(report_dir, exist_ok=True)

    h5ad_files = sorted(set(glob.glob(os.path.join(scan_path, "**", "*_harmonized.h5ad"), recursive=True)))
    logger.info(f"Scanning {len(h5ad_files)} harmonized h5ad files for summary report...")

    rows = []
    sample_rows = []
    warned_missing: set = set()
    for fpath in h5ad_files:
        try:
            ad  = sc.read_h5ad(fpath, backed="r")
            obs = ad.obs
            study_id = _get_obs_meta(obs, "study", warned_missing) or os.path.basename(fpath).replace(".h5ad", "")

            # Per-sample cell counts within this study (defaults to a single
            # "Unknown" sample when the column is absent).
            if "sample_id" in obs.columns:
                sid_counts = (
                    obs["sample_id"].astype(str).fillna("Unknown")
                       .replace({"": "Unknown", "nan": "Unknown"})
                       .value_counts()
                )
            else:
                if "sample_id" not in warned_missing:
                    logger.warning(
                        "Summary: 'sample_id' not found — per-sample plot will "
                        "treat each dataset as a single sample."
                    )
                    warned_missing.add("sample_id")
                sid_counts = pd.Series({study_id: int(ad.n_obs)})

            for sid, n in sid_counts.items():
                sample_rows.append({
                    "Study":     study_id,
                    "SampleID":  str(sid),
                    "Cells":     int(n),
                })

            rows.append({
                "File":       os.path.basename(fpath),
                "Path":       fpath,
                "Study":      study_id,
                "Samples":    int(len(sid_counts)),
                "Cells":      int(ad.n_obs),
                "Genes":      int(ad.n_vars),
                "Protocol":   _get_obs_meta(obs, "diff_protocol", warned_missing),
                "Age":        _get_obs_meta(obs, "Age", warned_missing),
                "Technology": _get_obs_meta(obs, "sc_protocol", warned_missing),
                "Source":     _get_obs_meta(obs, "source", warned_missing),
                "Disease":    _get_obs_meta(obs, "disease", warned_missing),
            })
        except Exception as e:
            logger.warning(f"Skipping summary for {fpath}: {e}")

    df = pd.DataFrame(rows)
    if df.empty:
        logger.warning("No datasets found for summary.")
        return

    df = df.sort_values(["Cells", "Study"], ascending=[False, True]).reset_index(drop=True)
    summary_csv = os.path.join(report_dir, "dataset_summary.csv")
    df.to_csv(summary_csv, index=False)
    logger.info(f"Summary table saved: {summary_csv}")

    # Per-sample table: ordered by study (largest study first), then by cells.
    df_samples = pd.DataFrame(sample_rows)
    if not df_samples.empty:
        study_order = df["Study"].tolist()
        df_samples["_study_rank"] = df_samples["Study"].map(
            {s: i for i, s in enumerate(study_order)}
        ).fillna(len(study_order)).astype(int)
        df_samples = (df_samples
                      .sort_values(["_study_rank", "Cells", "SampleID"],
                                   ascending=[True, False, True])
                      .drop(columns="_study_rank")
                      .reset_index(drop=True))
        samples_csv = os.path.join(report_dir, "samples_summary.csv")
        df_samples.to_csv(samples_csv, index=False)
        logger.info(f"Per-sample table saved: {samples_csv}")

    total_cells = int(df["Cells"].sum())
    
    print(f"\n{'='*88}")
    print(f"HUMAN KIDNEY ORGANOID ATLAS SUMMARY  (N={len(df)} datasets)")
    print(f"{'='*88}")
    print(f"Total cells          : {total_cells:,}")
    print(f"Mean genes / dataset : {df['Genes'].mean():,.0f}")
    print("\nTop datasets by cell count:")
    print(df[["Study","Samples","Cells","Genes","Protocol","Age","Technology"]].head(10).to_string(index=False))
    print("\n")

    def save_fig(fig, stem):
        for ext in opts["figure_extensions"]:
            path = os.path.join(report_dir, f"{stem}.{ext}")
            fig.savefig(path, facecolor="white", edgecolor="white", transparent=False)
        logger.info("Saved figure: %s.[%s]", stem, ", ".join(opts["figure_extensions"]))

    # ── Plot 1: cells per dataset
    plot_data = df.sort_values("Cells", ascending=True)
    colors    = sns.color_palette("tab20", n_colors=len(plot_data))
    fig, ax   = plt.subplots(figsize=(14, max(6, 0.35 * len(plot_data))))
    ax.barh(plot_data["Study"], plot_data["Cells"], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_title(f"Total Cells per Dataset  (Total: {total_cells:,})", fontweight="bold")
    ax.set_xlabel("Number of cells")
    ax.set_ylabel("Dataset / study")
    fig.tight_layout()
    save_fig(fig, "cells_per_dataset")
    plt.close(fig)

    # ── Plot 2: cells by organoid age (top N categories by cell count)
    age_top = (
        df.groupby("Age")["Cells"].sum()
          .sort_values(ascending=False)
          .reset_index()
          .head(opts["age_plot_top_n"])
    )
    if not age_top.empty:
        age_colors = sns.color_palette("crest", n_colors=len(age_top))
        fig, ax = plt.subplots(figsize=(12, max(5, 0.4 * len(age_top))))
        ax.barh(age_top["Age"][::-1], age_top["Cells"][::-1],
                color=age_colors, edgecolor="black", linewidth=0.4)
        ax.set_title(
            f"Cell Contribution by Organoid Age  "
            f"(Top {opts['age_plot_top_n']})",
            fontweight="bold",
        )
        ax.set_xlabel("Number of cells")
        ax.set_ylabel("Age / Timepoint")
        ax.grid(axis="x", linestyle="--", alpha=0.4)
        fig.tight_layout()
        save_fig(fig, "cells_by_age")
        plt.close(fig)

    # ── Plot 3: cells per sample, grouped & colored by study
    if not df_samples.empty:
        studies      = list(dict.fromkeys(df_samples["Study"].tolist()))
        study_palette = dict(zip(studies, sns.color_palette("tab20", n_colors=len(studies))))
        bar_colors   = [study_palette[s] for s in df_samples["Study"]]

        # Reverse so the first study (largest) appears at the top of the chart.
        labels = [f"{sid}  ({st})" for st, sid in
                  zip(df_samples["Study"][::-1], df_samples["SampleID"][::-1])]
        values = df_samples["Cells"][::-1]
        colors_rev = bar_colors[::-1]

        fig, ax = plt.subplots(figsize=(14, max(6, 0.28 * len(df_samples))))
        ax.barh(labels, values, color=colors_rev, edgecolor="black", linewidth=0.4)
        ax.set_title(
            f"Cells per Sample  "
            f"({len(df_samples)} samples across {len(studies)} datasets)",
            fontweight="bold",
        )
        ax.set_xlabel("Number of cells")
        ax.set_ylabel("Sample  (Study)")
        ax.grid(axis="x", linestyle="--", alpha=0.4)

        # One legend entry per study so the color mapping is readable.
        legend_handles = [
            plt.Rectangle((0, 0), 1, 1, color=study_palette[s], ec="black", lw=0.4)
            for s in studies
        ]
        ax.legend(legend_handles, studies, title="Study",
                  loc="lower right", frameon=True, fontsize="small")

        fig.tight_layout()
        save_fig(fig, "cells_per_sample")
        plt.close(fig)

    logger.info(f"All report artifacts written to: {report_dir}/")
