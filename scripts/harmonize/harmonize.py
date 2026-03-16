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

Only GTF path is in harmonize.config (stable). CSV and output are user-provided (--csv, --output or env):
Relative paths in the CSV’s first column are relative to working dir (see --working-dir).

"""

import os
import sys
import gzip
import argparse
import configparser
import collections
import warnings

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


# ══════════════════════════════════════════════════════════════════════════════
# CONFIG — harmonize.config (same dir as script) or --config PATH. CLI and env override.
# ══════════════════════════════════════════════════════════════════════════════

CONFIG_FILENAME = "harmonize.config"


def _config_path(cli_config: str | None) -> str:
    """Path to config file: CLI --config, or same dir as script, or cwd."""
    if cli_config and cli_config.strip():
        return cli_config.strip()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(script_dir, CONFIG_FILENAME)
    if os.path.isfile(p):
        return p
    return os.path.join(os.getcwd(), CONFIG_FILENAME)


def load_config(config_path: str) -> configparser.ConfigParser:
    """Load INI config. Returns ConfigParser with optionxform=str (preserve key case)."""
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    if os.path.isfile(config_path):
        cfg.read(config_path, encoding="utf-8")
    return cfg


def resolve_paths(args, cfg: configparser.ConfigParser) -> dict:
    """
    gtf_file: CLI > env > config (required in config if not passed). metadata_csv, output_root: CLI or env only (required, no fallback).
    working_dir: CLI --working_dir > env WORKING_DIR > current directory.
    """
    def get_gtf() -> str:
        if getattr(args, "gtf", None) and str(args.gtf).strip():
            return str(args.gtf).strip()
        v = os.environ.get("GTF_FILE", "").strip()
        if v:
            return v
        if cfg.has_section("paths") and cfg.has_option("paths", "gtf_file"):
            v = cfg.get("paths", "gtf_file").strip()
            if v:
                return v
        return ""

    def get_required(cli_val, env_key: str) -> str:
        """CSV and output: only from CLI or env; no fallback."""
        if cli_val is not None and str(cli_val).strip():
            return str(cli_val).strip()
        return os.environ.get(env_key, "").strip()

    working_dir = get_required(getattr(args, "working_dir", None), "WORKING_DIR") or os.getcwd()
    if not os.path.isabs(working_dir):
        working_dir = os.path.abspath(working_dir)

    return {
        "gtf_file":     get_gtf(),
        "metadata_csv": get_required(getattr(args, "csv", None), "METADATA_CSV"),
        "output_root":  get_required(getattr(args, "output", None), "OUTPUT_ROOT"),
        "working_dir":  working_dir,
    }


def get_summary_options(cfg: configparser.ConfigParser) -> dict:
    """Read [summary] options from config; return dict with defaults for missing keys."""
    defaults = {
        "figure_dpi": 300,
        "figure_extensions": ["png", "pdf"],
        "report_subdir": "reports/atlas_summary",
        "age_plot_top_n": 15,
    }
    if not cfg.has_section("summary"):
        return defaults
    out = dict(defaults)
    if cfg.has_option("summary", "figure_dpi"):
        try:
            out["figure_dpi"] = cfg.getint("summary", "figure_dpi")
        except ValueError:
            pass
    if cfg.has_option("summary", "figure_extensions"):
        raw = cfg.get("summary", "figure_extensions").strip()
        if raw:
            out["figure_extensions"] = [x.strip() for x in raw.split(",") if x.strip()] or defaults["figure_extensions"]
    if cfg.has_option("summary", "report_subdir"):
        v = cfg.get("summary", "report_subdir").strip()
        if v:
            out["report_subdir"] = v
    if cfg.has_option("summary", "age_plot_top_n"):
        try:
            out["age_plot_top_n"] = cfg.getint("summary", "age_plot_top_n")
        except ValueError:
            pass
    return out


# ══════════════════════════════════════════════════════════════════════════════


# ─────────────────────────────────────────────────────────────────────────────
# GTF PARSING
# ─────────────────────────────────────────────────────────────────────────────

def load_allowed_genes(gtf_path: str) -> set:
    """
    Parse Ensembl GTF and return the set of protein_coding + lncRNA gene names.
    Handles both plain .gtf and .gtf.gz files.
    """
    if not os.path.exists(gtf_path):
        raise FileNotFoundError(
            f"GTF file not found: {gtf_path!r}\n"
            "Set gtf_file in harmonize.config or GTF_FILE env var, or pass --gtf.\n"
            "Download from:\n"
            "  https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/"
            "Homo_sapiens.GRCh38.104.gtf.gz"
        )

    open_fn = gzip.open if gtf_path.endswith(".gz") else open
    mode    = "rt"     if gtf_path.endswith(".gz") else "r"
    genes   = []

    print(f"Parsing GTF: {os.path.basename(gtf_path)} ...")
    with open_fn(gtf_path, mode) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            info = {}
            for kv in cols[8].split("; "):
                parts = kv.strip().split(" ", 1)
                if len(parts) == 2:
                    info[parts[0]] = parts[1].strip('";')
            if info.get("gene_biotype") in ("protein_coding", "lncRNA"):
                g = info.get("gene_name")
                if g and pd.notna(g):
                    genes.append(g)

    allowed = set(genes)
    if not allowed:
        raise RuntimeError(
            "allowed_genes is empty after parsing the GTF. "
            "Check that the GTF has gene_biotype annotations."
        )
    print(f"GTF parsed: {len(allowed):,} protein_coding + lncRNA genes")
    return allowed


# ─────────────────────────────────────────────────────────────────────────────
# FILE READERS
# ─────────────────────────────────────────────────────────────────────────────

def _decode(arr) -> list:
    return [x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x) for x in arr]


def read_h5_safe(path: str) -> sc.AnnData:
    """Read a 10x H5 file with fallback for non-standard layouts."""
    try:
        return sc.read_10x_h5(path)
    except Exception as e1:
        print(f"  Standard H5 read failed ({e1}), trying fallback...")
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
                    except Exception:
                        continue
                if gene_names is None:
                    raise RuntimeError("Gene names not found in H5")
                X = csr_matrix(
                    (g["data"][()], g["indices"][()], g["indptr"][()]),
                    shape=tuple(g["shape"][()])
                ).T
            return sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes),
                              var=pd.DataFrame(index=gene_names))
        except Exception as e2:
            raise RuntimeError(f"Cannot read H5: {path}\n  err1={e1}\n  err2={e2}")


def read_mtx_safe(folder: str, prefix: str = "") -> sc.AnnData:
    """Read a 10x MTX triplet (with optional filename prefix) from a folder."""
    try:
        return sc.read_10x_mtx(folder, var_names="gene_symbols", prefix=prefix)
    except Exception as e:
        print(f"  sc.read_10x_mtx failed ({e}), trying manual load...")
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
        return sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes[0].values),
                          var=pd.DataFrame(index=gene_names))


def geo_redownload(folder: str, prefix: str, sample_id: str,
                   feature_file: str = "features.tsv.gz") -> None:
    """Re-download MTX triplet files from GEO FTP when local files are corrupted.
    Raises on first download failure so the caller does not retry read with bad files.
    GEO sample IDs are typically GSM + 6+ digits (e.g. GSM1234567).
    """
    import urllib.request
    if len(sample_id) < 9:
        raise ValueError(f"GEO re-download requires sample_id of length >= 9 (e.g. GSM1234567), got {sample_id!r}")
    gsm_prefix = f"{sample_id[:-3]}nnn"
    base_url   = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm_prefix}/{sample_id}/suppl"
    for fname in [f"{prefix}matrix.mtx.gz", f"{prefix}barcodes.tsv.gz",
                  f"{prefix}{feature_file}"]:
        url   = f"{base_url}/{fname}"
        fpath = os.path.join(folder, fname)
        if os.path.exists(fpath):
            os.replace(fpath, fpath + ".corrupt")
        print(f"  Redownloading {url}...")
        try:
            urllib.request.urlretrieve(url, fpath)
        except Exception as e:
            raise RuntimeError(f"GEO re-download failed for {fname}: {e}") from e


# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE LOADER
# ─────────────────────────────────────────────────────────────────────────────

META_COLS = [
    "sample_id", "study", "source", "diff_protocol", "sc_protocol",
    "sequencing", "genome_build", "Age", "type", "disease", "condition",
]


def load_sample(row: pd.Series, working_dir: str | None = None) -> sc.AnnData:
    """
    Load one sample from a CSV metadata row.
    Auto-detects format from data_path:
        .h5ad  → read_h5ad
        .h5    → read_h5_safe (10x H5)
        folder → read_mtx_safe (MTX triplet, with optional prefix)
    If data_path is relative, it is resolved against working_dir (where data is stored).
    Returns AnnData with all metadata columns attached to obs.
    """
    path = str(row["data_path"]).strip()
    if not path:
        raise ValueError("data_path is empty")
    if working_dir and not os.path.isabs(path):
        path = os.path.normpath(os.path.join(working_dir, path))
    prefix    = str(row.get("file_prefix", "") or "").strip()
    if prefix.lower() in ("nan", "none"):
        prefix = ""
    sample_id = str(row.get("sample_id", "") or "").strip()

    if not os.path.exists(path):
        raise FileNotFoundError(f"data_path not found: {path}")

    if os.path.isfile(path):
        ext = path.lower()
        if ext.endswith(".h5ad"):
            print(f"  [{sample_id}] Reading h5ad: {os.path.basename(path)}")
            adata = sc.read_h5ad(path)
        elif ext.endswith(".h5") or ext.endswith(".hdf5"):
            print(f"  [{sample_id}] Reading H5: {os.path.basename(path)}")
            adata = read_h5_safe(path)
        else:
            raise ValueError(f"Unknown file extension: {path}")
    elif os.path.isdir(path):
        print(f"  [{sample_id}] Reading MTX: {os.path.basename(path)} (prefix={prefix!r})")
        try:
            adata = read_mtx_safe(path, prefix=prefix)
        except Exception as e1:
            if sample_id.startswith("GSM"):
                print("  Read failed, attempting GEO re-download...")
                try:
                    geo_redownload(path, prefix, sample_id)
                except Exception as e2:
                    raise RuntimeError(
                        f"MTX read failed: {e1}; GEO re-download failed: {e2}"
                    ) from e2
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

    return adata


# ─────────────────────────────────────────────────────────────────────────────
# HARMONIZATION
# ─────────────────────────────────────────────────────────────────────────────

def harmonize_matrix_sparse(adata: sc.AnnData, allowed_genes: set) -> sc.AnnData:
    """
    Harmonize an AnnData to the allowed_genes reference set using native sparse ops.
    - Genes in the dataset but not in allowed_genes are dropped.
    - Genes in allowed_genes but not in the dataset are zero-filled.
    - Output matrix is float32 CSR.
    """
    target_genes  = sorted(allowed_genes)
    dataset_genes = adata.var_names.tolist()
    common_genes  = sorted(set(dataset_genes) & set(target_genes))

    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)

    adata_common = adata[:, common_genes].copy()
    gene_to_idx  = {g: i for i, g in enumerate(target_genes)}
    coo          = sp.coo_matrix(adata_common.X)
    new_cols     = np.array([gene_to_idx[g] for g in common_genes])[coo.col]
    new_X        = sp.coo_matrix(
        (coo.data, (coo.row, new_cols)),
        shape=(adata.shape[0], len(target_genes))
    ).tocsr().astype(np.float32)

    return sc.AnnData(
        X=new_X,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=target_genes),
    )


def memory_efficient_harmonize_sparse(input_file: str, output_file: str,
                                       allowed_genes: set) -> None:
    """Read a raw h5ad, harmonize it, and write the harmonized h5ad."""
    print(f"  Harmonizing {os.path.basename(input_file)} ...")
    adata = sc.read_h5ad(input_file)
    print(f"    Initial shape : {adata.shape}")
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata_new = harmonize_matrix_sparse(adata, allowed_genes)
    print(f"    Harmonized shape: {adata_new.shape}")
    adata_new.write(output_file)
    print(f"    Saved to {output_file}")


# ─────────────────────────────────────────────────────────────────────────────
# CSV LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_metadata_csv(csv_path: str) -> pd.DataFrame:
    """
    Load and validate the metadata CSV.
    - Auto-detects comma vs semicolon separator.
    - Strips column name whitespace.
    - Keeps empty strings as '' (not NaN).
    - Validates required columns.
    - Applies skip column if present.
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path!r}  (cwd: {os.getcwd()})")

    with open(csv_path, encoding="utf-8") as fh:
        header = fh.readline()
    sep = ";" if header.count(";") > header.count(",") else ","
    if sep == ";":
        print("WARNING: CSV uses semicolons — reading with sep=';'. "
              "Re-save as comma-separated to avoid issues.")

    df = pd.read_csv(csv_path, sep=sep, keep_default_na=False, encoding="utf-8")
    df.columns = df.columns.str.strip()

    required = ["data_path", "file_prefix", "sample_id", "output_dir", "study"]
    missing  = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"CSV missing required columns: {missing}\n"
            f"Found: {list(df.columns)}\nFile: {csv_path!r}"
        )

    if "skip" in df.columns:
        n_skip = (df["skip"].astype(str).str.strip().str.lower() == "true").sum()
        if n_skip:
            print(f"Skipping {n_skip} row(s) with skip=True")
        df = df[df["skip"].astype(str).str.strip().str.lower() != "true"]

    if df.empty:
        raise ValueError(f"CSV has no rows left after skip filter. File: {csv_path!r}")

    # Require non-empty data_path (first column)
    empty_path = df["data_path"].astype(str).str.strip() == ""
    if empty_path.any():
        rows = (df.index[empty_path] + 2).tolist()  # 2 = header + 1-based
        raise ValueError(
            f"CSV has empty data_path in row(s) {rows}\n"
            f"File: {csv_path!r}"
        )

    # Require non-empty study
    empty_study = df["study"].astype(str).str.strip() == ""
    if empty_study.any():
        rows = (df.index[empty_study] + 2).tolist()
        raise ValueError(
            f"CSV has empty study in row(s) {rows}\n"
            f"File: {csv_path!r}"
        )

    # Require unique data_path
    path_series = df["data_path"].astype(str).str.strip()
    path_dups = path_series[path_series.duplicated(keep=False)]
    if not path_dups.empty:
        dup_paths = sorted(path_dups.unique().tolist())
        raise ValueError(
            f"CSV has duplicate data_path: {dup_paths[:5]}{'...' if len(dup_paths) > 5 else ''}\n"
            f"File: {csv_path!r}"
        )

    # Require unique sample_id
    id_series = df["sample_id"].astype(str).str.strip()
    id_dups = id_series[id_series.duplicated(keep=False)]
    if not id_dups.empty:
        dup_ids = sorted(id_dups.unique().tolist())
        raise ValueError(
            f"CSV has duplicate sample_id: {dup_ids}\n"
            f"File: {csv_path!r}"
        )

    return df


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(metadata_csv: str, gtf_file: str,
                 output_root: str, working_dir: str) -> list:
    """
    Main pipeline loop.
    working_dir: base for relative data_path entries in the CSV.
    Returns list of dicts describing any failures.
    """
    allowed_genes = load_allowed_genes(gtf_file)

    df_meta = load_metadata_csv(metadata_csv)
    print(f"CSV loaded: {len(df_meta)} sample rows, "
          f"{df_meta['study'].nunique()} unique studies")
    print(f"Columns: {list(df_meta.columns)}\n")

    study_groups = collections.OrderedDict()
    for _, row in df_meta.iterrows():
        study_groups.setdefault(row["study"], []).append(row)

    failed = []

    for study, rows in study_groups.items():
        print("=" * 64)
        print(f"Processing study: {study} ({len(rows)} sample(s))")

        if output_root:
            out_dir = os.path.join(output_root, study)
        else:
            out_dir = str(rows[0].get("output_dir", study)).strip() or study
            if out_dir and not os.path.isabs(out_dir):
                out_dir = os.path.normpath(os.path.join(working_dir, out_dir))
        os.makedirs(out_dir, exist_ok=True)
        raw_h5ad        = os.path.join(out_dir, f"{study}.h5ad")
        harmonized_h5ad = os.path.join(out_dir, f"{study}_harmonized.h5ad")

        # ── Load all samples (relative data_path resolved against working_dir)
        adatas = []
        for row in rows:
            try:
                adata = load_sample(row, working_dir=working_dir)
                adatas.append(adata)
                print(f"  Loaded {row.get('sample_id', '?')} → {adata.shape}")
            except Exception as e:
                print(f"  FAILED {row.get('sample_id', '?')}: {e}")
                failed.append({"study": study,
                               "sample_id": row.get("sample_id", "?"),
                               "error": str(e)})

        if not adatas:
            print(f"  No samples loaded for {study}, skipping.")
            continue

        # ── Concatenate
        if len(adatas) == 1:
            adata = adatas[0]
        else:
            adata = sc.concat(adatas, join="outer")
            adata.obs_names_make_unique()
            print(f"  Concatenated {len(adatas)} samples → {adata.shape}")

        # Fix string dtypes
        for col in adata.obs.select_dtypes(include=["string", "object"]).columns:
            adata.obs[col] = adata.obs[col].astype("object")

        # ── Save raw
        adata.write(raw_h5ad)
        print(f"  Saved raw: {raw_h5ad}  shape={adata.shape}")

        # ── Harmonize
        try:
            memory_efficient_harmonize_sparse(raw_h5ad, harmonized_h5ad, allowed_genes)
        except Exception as e:
            print(f"  Harmonization failed: {e}")
            failed.append({"study": study, "sample_id": "harmonize", "error": str(e)})

    print()
    print("=" * 64)
    if failed:
        print(f"DONE with {len(failed)} failure(s):")
        for f in failed:
            print(f"  [{f['study']} / {f['sample_id']}] {f['error']}")
    else:
        print("ALL DONE — no failures.")

    return failed


# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY & PLOTS
# ─────────────────────────────────────────────────────────────────────────────

def run_summary(scan_path: str, cfg: configparser.ConfigParser | None = None) -> None:
    """Scan harmonized h5ad files and produce summary table + publication plots."""
    import glob
    import matplotlib
    matplotlib.use("Agg")          # non-interactive backend — safe on servers
    import matplotlib.pyplot as plt
    import seaborn as sns

    opts = get_summary_options(cfg or configparser.ConfigParser())
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.35)
    plt.rcParams["figure.dpi"]    = opts["figure_dpi"]
    plt.rcParams["savefig.dpi"]   = opts["figure_dpi"]
    plt.rcParams["savefig.bbox"]  = "tight"

    report_dir = os.path.join(scan_path, opts["report_subdir"])
    os.makedirs(report_dir, exist_ok=True)

    # ── Scan files (harmonized only; raw .h5ad have variable gene counts)
    h5ad_files = sorted(set(glob.glob(
        os.path.join(scan_path, "**", "*_harmonized.h5ad"), recursive=True
    )))
    print(f"Scanning {len(h5ad_files)} harmonized h5ad files...")

    rows = []
    for fpath in h5ad_files:
        try:
            ad  = sc.read_h5ad(fpath, backed="r")
            obs = ad.obs

            def get_meta(col):
                if col not in obs.columns:
                    return "Unknown"
                vals = [str(v) for v in obs[col].unique()
                        if pd.notna(v) and str(v).lower() != "nan"]
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
            print(f"  Skipping {fpath}: {e}")

    df = pd.DataFrame(rows)
    if df.empty:
        print("No datasets found.")
        return

    df = df.sort_values(["Cells", "Study"], ascending=[False, True]).reset_index(drop=True)

    # ── Save summary CSV
    summary_csv = os.path.join(report_dir, "dataset_summary.csv")
    df.to_csv(summary_csv, index=False)
    print(f"Saved: {summary_csv}")

    total_cells = int(df["Cells"].sum())
    print(f"\n{'='*88}")
    print(f"HUMAN KIDNEY ORGANOID ATLAS SUMMARY  (N={len(df)} datasets)")
    print(f"{'='*88}")
    print(f"Total cells          : {total_cells:,}")
    print(f"Mean genes / dataset : {df['Genes'].mean():,.0f}")
    print("\nTop datasets by cell count:")
    print(df[["Study","Cells","Genes","Protocol","Age","Technology"]].head(10).to_string(index=False))

    def save_fig(fig, stem):
        for ext in opts["figure_extensions"]:
            path = os.path.join(report_dir, f"{stem}.{ext}")
            fig.savefig(path)
        print(f"Saved figure: {os.path.join(report_dir, stem)}.{', '.join(opts['figure_extensions'])}")

    # ── Plot 1: cells per dataset
    plot_data = df.sort_values("Cells", ascending=True)
    colors    = plt.cm.viridis(np.linspace(0.15, 0.9, len(plot_data)))
    fig, ax   = plt.subplots(figsize=(14, max(6, 0.35 * len(plot_data))))
    ax.barh(plot_data["Study"], plot_data["Cells"], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_title(f"Total Cells per Dataset  (Total: {total_cells:,})", fontweight="bold")
    ax.set_xlabel("Number of cells")
    ax.set_ylabel("Dataset / study")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    for y, v in enumerate(plot_data["Cells"].values):
        ax.text(v, y, f" {v:,}", va="center", fontsize=8)
    fig.tight_layout()
    save_fig(fig, "cells_per_dataset")
    plt.close(fig)

    # ── Plot 2: cells by differentiation protocol
    proto = df.groupby("Protocol")["Cells"].sum().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.pie(proto.values, labels=proto.index, autopct="%1.1f%%", startangle=140,
           colors=sns.color_palette("pastel", n_colors=len(proto)),
           wedgeprops={"edgecolor": "black", "linewidth": 0.6},
           textprops={"fontsize": 10})
    ax.add_artist(plt.Circle((0, 0), 0.70, fc="white"))
    ax.set_title("Cell Contribution by Differentiation Protocol", fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "cells_by_protocol")
    plt.close(fig)

    # ── Plot 3: cells by organoid age
    age_top = (df.groupby("Age")["Cells"].sum()
                 .sort_values(ascending=False)
                 .reset_index()
                 .head(opts["age_plot_top_n"]))
    colors  = plt.cm.magma(np.linspace(0.15, 0.9, len(age_top)))
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.barh(age_top["Age"][::-1], age_top["Cells"][::-1],
            color=colors, edgecolor="black", linewidth=0.4)
    ax.set_title("Cell Contribution by Organoid Age (Timepoint)", fontweight="bold")
    ax.set_xlabel("Number of cells")
    ax.set_ylabel("Age")
    ax.grid(axis="x", linestyle="--", alpha=0.4)
    fig.tight_layout()
    save_fig(fig, "cells_by_age")
    plt.close(fig)

    print(f"\nAll report artifacts written to: {report_dir}/")


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="scRNA-seq atlas harmonization pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--config",  default=None, help="Path to harmonize.config (default: same dir as script)")
    p.add_argument("--working-dir", "-w", dest="working_dir", default=None, help="Working directory; CSV path and paths in CSV are relative to this (default: cwd)")
    p.add_argument("--csv",     default=None, help="Path to metadata CSV (relative to working dir if not absolute)")
    p.add_argument("--gtf",     default=None, help="Path to GRCh38 GTF (overrides config/env)")
    p.add_argument("--output",  default=None, help="Output root directory (relative to working dir if not absolute)")
    p.add_argument("--summary", action="store_true", help="Generate summary plots after pipeline")
    p.add_argument("--summary-only", action="store_true",
                   help="Skip pipeline, only generate summary plots from existing h5ad files")
    return p.parse_args()


def main():
    args = parse_args()

    cfg = load_config(_config_path(args.config))
    paths = resolve_paths(args, cfg)
    working_dir = paths["working_dir"]
    gtf_file = paths["gtf_file"]
    metadata_csv = paths["metadata_csv"]
    output_root = paths["output_root"]

    # Require inputs: GTF from config (or --gtf / GTF_FILE), CSV and output from user
    missing = []
    if not gtf_file:
        missing.append("GTF path (set gtf_file in harmonize.config or pass --gtf or set GTF_FILE)")
    if not args.summary_only:
        if not metadata_csv:
            missing.append("Metadata CSV (pass --csv or set METADATA_CSV)")
        if not output_root:
            missing.append("Output directory (pass --output or set OUTPUT_ROOT)")
    elif not output_root:
        missing.append("Output directory for --summary-only (pass --output or set OUTPUT_ROOT)")

    if missing:
        print("Error: required input is missing:", file=sys.stderr)
        for m in missing:
            print(f"  - {m}", file=sys.stderr)
        sys.exit(1)

    # Resolve CSV and output relative to working directory
    if metadata_csv and not os.path.isabs(metadata_csv):
        metadata_csv = os.path.normpath(os.path.join(working_dir, metadata_csv))
    if output_root and not os.path.isabs(output_root):
        output_root = os.path.normpath(os.path.join(working_dir, output_root))

    print(f"Working dir   : {working_dir}")
    print(f"GTF file     : {gtf_file}")
    print(f"Metadata CSV : {metadata_csv}")
    print(f"Output root  : {output_root}\n")

    if not args.summary_only:
        failed = run_pipeline(
            metadata_csv=metadata_csv,
            gtf_file=gtf_file,
            output_root=output_root,
            working_dir=working_dir,
        )
        if failed:
            sys.exit(1)

    if args.summary or args.summary_only:
        print("\nGenerating summary plots...")
        run_summary(scan_path=output_root or ".", cfg=cfg)


if __name__ == "__main__":
    main()
