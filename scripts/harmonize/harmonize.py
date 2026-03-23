#!/usr/bin/env python3
"""
scRNA-seq Atlas Harmonization Script
=======================================
Reads a metadata CSV, loads each sample (H5 / H5AD / MTX),
concatenates samples per study, then harmonizes gene space to
the 24,100 protein-coding + lncRNA genes from GRCh38.104.

Usage
-----
    python harmonize.py [-w DIR] [--csv PATH] [--gtf PATH] [--output PATH] [--summary]

Outputs:
- Harmonized .h5ad files per study
- harmonize_pipeline.log (Detailed debug log)
"""

import os
import sys
import gzip
import argparse
import configparser
import collections
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


# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════

CONFIG_FILENAME = "harmonize.config"

def _config_path(cli_config: str | None) -> str:
    if cli_config and cli_config.strip():
        return cli_config.strip()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(script_dir, CONFIG_FILENAME)
    if os.path.isfile(p):
        return p
    return os.path.join(os.getcwd(), CONFIG_FILENAME)

def load_config(config_path: str) -> configparser.ConfigParser:
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    if os.path.isfile(config_path):
        cfg.read(config_path, encoding="utf-8")
    return cfg

def resolve_paths(args, cfg: configparser.ConfigParser) -> dict:
    def get_path(cli_val, cfg_key: str, env_key: str, default: str = "") -> str:
        if cli_val is not None and str(cli_val).strip():
            return str(cli_val).strip()
        if cfg.has_section("paths") and cfg.has_option("paths", cfg_key):
            v = cfg.get("paths", cfg_key).strip()
            if v:
                return v
        v = os.environ.get(env_key, "").strip()
        if v:
            return v
        return default

    working_dir = get_path(getattr(args, "working_dir", None), "working_dir", "WORKING_DIR", os.getcwd())
    if not os.path.isabs(working_dir):
        working_dir = os.path.abspath(working_dir)

    return {
        "gtf_file":     get_path(getattr(args, "gtf", None), "gtf_file", "GTF_FILE"),
        "metadata_csv": get_path(getattr(args, "csv", None), "metadata_csv", "METADATA_CSV"),
        "output_root":  get_path(getattr(args, "output", None), "output_root", "OUTPUT_ROOT"),
        "working_dir":  working_dir,
    }

def get_summary_options(cfg: configparser.ConfigParser) -> dict:
    defaults = {
        "figure_dpi": 300,
        "figure_extensions": ["png", "pdf"],
        "report_subdir": "reports/atlas_summary",
        "age_plot_top_n": 15,
    }
    if not cfg.has_section("summary"): return defaults
    out = dict(defaults)
    if cfg.has_option("summary", "figure_dpi"):
        try: out["figure_dpi"] = cfg.getint("summary", "figure_dpi")
        except ValueError: pass
    if cfg.has_option("summary", "figure_extensions"):
        raw = cfg.get("summary", "figure_extensions").strip()
        if raw: out["figure_extensions"] = [x.strip() for x in raw.split(",") if x.strip()] or defaults["figure_extensions"]
    if cfg.has_option("summary", "report_subdir"):
        v = cfg.get("summary", "report_subdir").strip()
        if v: out["report_subdir"] = v
    if cfg.has_option("summary", "age_plot_top_n"):
        try: out["age_plot_top_n"] = cfg.getint("summary", "age_plot_top_n")
        except ValueError: pass
    return out


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

META_COLS = [
    "sample_id", "study", "source","species", "tissue", "diff_protocol", "sc_protocol",
    "sequencing", "genome_build", "Age", "type", "disease", "condition",
]

def load_sample(row: pd.Series, working_dir: str | None = None) -> sc.AnnData:
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

    for col in META_COLS:
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
    
    # === NEW TRANSGENE DISCOVERY STEP ===
    sample_genes = set(adata.var_names.tolist())
    common_genes = sorted(sample_genes & allowed_genes)
    
    extra_genes = sample_genes.difference(allowed_genes)
    if extra_genes:
        transgene_patterns = ('AAV', 'GFP', 'MCHERRY', 'TDTOMATO', 'ERCC-', 'HTO', 'CRE')
        potential_transgenes = [
            g for g in extra_genes 
            if any(pattern in g.upper() for pattern in transgene_patterns)
        ]
        if potential_transgenes:
            logger.warning(
                f"⚠️ ALERT: Dropping {len(potential_transgenes)} potential transgenes/reporters "
                f"because they are absent from the GTF: {potential_transgenes[:10]}..."
            )
            logger.warning(
                "If you need to preserve these, restart the pipeline and pass them "
                f"via the CLI: --transgenes {','.join(potential_transgenes[:3])}"
            )
    # ====================================

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

    required = ["data_path", "sample_id", "study"]
    missing  = [c for c in required if c not in df.columns]
    if missing:
        logger.error(f"CSV missing required columns: {missing}")
        raise ValueError(f"CSV missing required columns: {missing}")

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
        from rpy2.robjects import pandas2ri
        import anndata2ri
        import gc
        
        pandas2ri.activate()
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
        
        with ro.conversion.localconverter(anndata2ri.converter):
            ro.globalenv["adata_sce"] = clean_adata

        ro.r(f'''
            suppressPackageStartupMessages(library(Seurat))
            suppressPackageStartupMessages(library(SingleCellExperiment))
            suppressPackageStartupMessages(library(Matrix))

            counts_mat <- assay(adata_sce, assayNames(adata_sce)[1])
            
            seurat_obj <- CreateSeuratObject(counts = counts_mat, assay = "RNA")
            seurat_obj[[]] <- as.data.frame(colData(adata_sce))
            seurat_obj[["nCount_RNA"]] <- Matrix::colSums(counts_mat)
            seurat_obj[["nFeature_RNA"]] <- Matrix::colSums(counts_mat > 0)

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
        
        if to_rds:
            out_dir_rds = os.path.join(out_dir, "rds")
            os.makedirs(out_dir_rds, exist_ok=True)
            rds_file = os.path.join(out_dir_rds, f"{study}_harmonized.rds")

        # ── Load all samples
        adatas = []
        for row in rows:
            sample_id = row.get('sample_id', '?')
            try:
                adata = load_sample(row, working_dir=working_dir)
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

def run_summary(scan_path: str, cfg: configparser.ConfigParser | None = None) -> None:
    import glob
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    opts = get_summary_options(cfg or configparser.ConfigParser())
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.35)
    plt.rcParams["figure.dpi"]    = opts["figure_dpi"]
    plt.rcParams["savefig.dpi"]   = opts["figure_dpi"]
    plt.rcParams["savefig.bbox"]  = "tight"

    report_dir = os.path.join(scan_path, opts["report_subdir"])
    os.makedirs(report_dir, exist_ok=True)

    h5ad_files = sorted(set(glob.glob(os.path.join(scan_path, "**", "*_harmonized.h5ad"), recursive=True)))
    logger.info(f"Scanning {len(h5ad_files)} harmonized h5ad files for summary report...")

    rows = []
    for fpath in h5ad_files:
        try:
            ad  = sc.read_h5ad(fpath, backed="r")
            obs = ad.obs

            def get_meta(col):
                if col not in obs.columns: return "Unknown"
                vals = [str(v) for v in obs[col].unique() if pd.notna(v) and str(v).lower() != "nan"]
                return vals[0] if len(vals) == 1 else (", ".join(sorted(vals)) if vals else "Unknown")

            study_id = get_meta("study") or os.path.basename(fpath).replace(".h5ad", "")
            rows.append({
                "File":       os.path.basename(fpath),
                "Path":       fpath,
                "Study":      study_id,
                "Cells":      int(ad.n_obs),
                "Genes":      int(ad.n_vars),
                "Protocol":   get_meta("diff_protocol"),
                "Age":        get_meta("Age"),
                "Technology": get_meta("sc_protocol"),
                "Source":     get_meta("source"),
                "Disease":    get_meta("disease"),
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

    total_cells = int(df["Cells"].sum())
    
    print(f"\n{'='*88}")
    print(f"HUMAN KIDNEY ORGANOID ATLAS SUMMARY  (N={len(df)} datasets)")
    print(f"{'='*88}")
    print(f"Total cells          : {total_cells:,}")
    print(f"Mean genes / dataset : {df['Genes'].mean():,.0f}")
    print("\nTop datasets by cell count:")
    print(df[["Study","Cells","Genes","Protocol","Age","Technology"]].head(10).to_string(index=False))
    print("\n")

    def save_fig(fig, stem):
        for ext in opts["figure_extensions"]:
            path = os.path.join(report_dir, f"{stem}.{ext}")
            fig.savefig(path)
        logger.info(f"Saved figure: {stem}.{opts['figure_extensions']}")

    plot_data = df.sort_values("Cells", ascending=True)
    colors    = plt.cm.viridis(np.linspace(0.15, 0.9, len(plot_data)))
    fig, ax   = plt.subplots(figsize=(14, max(6, 0.35 * len(plot_data))))
    ax.barh(plot_data["Study"], plot_data["Cells"], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_title(f"Total Cells per Dataset  (Total: {total_cells:,})", fontweight="bold")
    ax.set_xlabel("Number of cells")
    ax.set_ylabel("Dataset / study")
    fig.tight_layout()
    save_fig(fig, "cells_per_dataset")
    plt.close(fig)

    logger.info(f"All report artifacts written to: {report_dir}/")


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="scRNA-seq atlas harmonization pipeline")
    p.add_argument("--config",  default=None, help="Path to harmonize.config")
    p.add_argument("--working-dir", "-w", dest="working_dir", default=None, help="Working directory")
    p.add_argument("--csv",     default=None, help="Path to metadata CSV")
    p.add_argument("--gtf",     default=None, help="Path to GRCh38 GTF")
    p.add_argument("--output",  default=None, help="Output root directory")
    p.add_argument("--summary", action="store_true", help="Generate summary plots after pipeline")
    p.add_argument("--summary-only", action="store_true", help="Skip pipeline, only generate plots")
    p.add_argument("--to-rds", action="store_true", help="Convert harmonized .h5ad to Seurat .rds")
    p.add_argument("--transgenes", default=None,
                   help="Comma-separated transgene names to add to the reference gene set "
                        "(e.g. --transgenes GFP,Cre,tdTomato). "
                        "Can also be set in harmonize.config under [transgenes] names = ...")
    return p.parse_args()

def main():
    args = parse_args()
    
    cfg = load_config(_config_path(args.config))
    paths = resolve_paths(args, cfg)
    working_dir  = paths["working_dir"]
    gtf_file     = paths["gtf_file"]
    metadata_csv = paths["metadata_csv"]
    output_root  = paths["output_root"]

    setup_logging(working_dir)

    missing = []
    if not gtf_file:
        missing.append("GTF path")
    if not args.summary_only:
        if not metadata_csv: missing.append("Metadata CSV")
        if not output_root:  missing.append("Output directory")
    elif not output_root:
        missing.append("Output directory for --summary-only")

    if missing:
        logger.error(f"Missing required inputs: {', '.join(missing)}")
        sys.exit(1)

    if metadata_csv and not os.path.isabs(metadata_csv):
        metadata_csv = os.path.normpath(os.path.join(working_dir, metadata_csv))
    if gtf_file and not os.path.isabs(gtf_file):
        gtf_file = os.path.normpath(os.path.join(working_dir, gtf_file))
    if output_root and not os.path.isabs(output_root):
        output_root = os.path.normpath(os.path.join(working_dir, output_root))

    logger.info("=== Atlas Harmonization Initialized ===")
    logger.info(f"Working dir  : {working_dir}")
    logger.info(f"GTF file     : {gtf_file}")
    logger.info(f"Metadata CSV : {metadata_csv}")
    logger.info(f"Output root  : {output_root}")

    # ── Resolve transgene names: CLI > config file
    transgene_names = set()
    if getattr(args, "transgenes", None):
        transgene_names = {t.strip() for t in args.transgenes.split(",") if t.strip()}
    elif cfg.has_section("transgenes") and cfg.has_option("transgenes", "names"):
        raw = cfg.get("transgenes", "names").strip()
        transgene_names = {t.strip() for t in raw.split(",") if t.strip()}
    if transgene_names:
        logger.info(f"Transgenes added to reference: {sorted(transgene_names)}")

    if not args.summary_only:
        failed = run_pipeline(
            metadata_csv=metadata_csv,
            gtf_file=gtf_file,
            output_root=output_root,
            working_dir=working_dir,
            to_rds=args.to_rds,
            transgene_names=transgene_names,
        )
        if failed:
            sys.exit(1)

    if args.summary or args.summary_only:
        logger.info("Generating summary plots...")
        run_summary(scan_path=output_root or ".", cfg=cfg)

if __name__ == "__main__":
    main()