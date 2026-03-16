#!/usr/bin/env python3
"""
scRNA-seq Atlas Harmonization Pipeline
=======================================
Reads a metadata CSV, loads each sample (H5 / H5AD / MTX),
concatenates samples per study, then harmonizes gene space to
the 24,100 protein-coding + lncRNA genes from GRCh38.104.

Usage
-----
    python harmonize_pipeline.py [--csv PATH] [--gtf PATH] [--output PATH] [--summary]

All arguments are optional — defaults are read from the USER CONFIG section below,
or from environment variables:
    METADATA_CSV   path to the metadata CSV
    GTF_FILE       path to the Ensembl GRCh38.104 GTF
    DATA_ROOT      working directory (relative CSV paths resolve from here)
    OUTPUT_ROOT    where to write output h5ad files

Examples
--------
    # Use defaults from USER CONFIG
    python harmonize_pipeline.py

    # Override CSV and output dir on the command line
    python harmonize_pipeline.py --csv my_datasets.csv --output /data/results

    # Run pipeline then generate summary plots
    python harmonize_pipeline.py --summary
"""

import os
import sys
import gzip
import argparse
import collections
import warnings
import traceback

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
# USER CONFIG — edit these four paths, or set them as environment variables
# ══════════════════════════════════════════════════════════════════════════════

DATA_ROOT = os.environ.get(
    "DATA_ROOT",
    "/data-master/pure-workspace/labss/hmami/new_data/data",
)

GTF_FILE = os.environ.get(
    "GTF_FILE",
    "/data-master/pure-workspace/labss/hmami/new_data/meta/gtf/Homo_sapiens.GRCh38.104.gtf",
)

METADATA_CSV = os.environ.get(
    "METADATA_CSV",
    "/data-master/pure-workspace/labss/hmami/new_data/notebooks/data_loaders/datasets_metadata.csv",
)

OUTPUT_ROOT = os.environ.get(
    "OUTPUT_ROOT",
    "/data-master/pure-workspace/labss/hmami/new_data/notebooks/data_loaders/CSV_driven_results",
)

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
            "Set GTF_FILE in USER CONFIG or via environment variable.\n"
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
    """Re-download MTX triplet files from GEO FTP when local files are corrupted."""
    import urllib.request
    gsm_prefix = f"{sample_id[:-3]}nnn"
    base_url   = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm_prefix}/{sample_id}/suppl"
    for fname in [f"{prefix}matrix.mtx.gz", f"{prefix}barcodes.tsv.gz",
                  f"{prefix}{feature_file}"]:
        url   = f"{base_url}/{fname}"
        fpath = os.path.join(folder, fname)
        if os.path.exists(fpath):
            os.replace(fpath, fpath + ".corrupt")
        try:
            print(f"  Redownloading {url}...")
            urllib.request.urlretrieve(url, fpath)
        except Exception as e:
            print(f"  Download failed: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE LOADER
# ─────────────────────────────────────────────────────────────────────────────

META_COLS = [
    "sample_id", "study", "source", "diff_protocol", "sc_protocol",
    "sequencing", "genome_build", "Age", "type", "disease", "condition",
]


def load_sample(row: pd.Series) -> sc.AnnData:
    """
    Load one sample from a CSV metadata row.
    Auto-detects format from data_path:
        .h5ad  → read_h5ad
        .h5    → read_h5_safe (10x H5)
        folder → read_mtx_safe (MTX triplet, with optional prefix)
    Returns AnnData with all metadata columns attached to obs.
    """
    path      = str(row["data_path"]).strip()
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
        except Exception:
            if sample_id.startswith("GSM"):
                print("  Read failed, attempting GEO re-download...")
                geo_redownload(path, prefix, sample_id)
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

    with open(csv_path) as fh:
        header = fh.readline()
    sep = ";" if header.count(";") > header.count(",") else ","
    if sep == ";":
        print("WARNING: CSV uses semicolons — reading with sep=';'. "
              "Re-save as comma-separated to avoid issues.")

    df = pd.read_csv(csv_path, sep=sep, keep_default_na=False)
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

    return df


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def run_pipeline(metadata_csv: str, gtf_file: str,
                 output_root: str) -> list:
    """
    Main pipeline loop.
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

        out_dir = os.path.join(output_root, study) if output_root \
                  else str(rows[0].get("output_dir", study)).strip() or study
        os.makedirs(out_dir, exist_ok=True)
        raw_h5ad        = os.path.join(out_dir, f"{study}.h5ad")
        harmonized_h5ad = os.path.join(out_dir, f"{study}_harmonized.h5ad")

        # ── Load all samples
        adatas = []
        for row in rows:
            try:
                adata = load_sample(row)
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

def run_summary(scan_path: str) -> None:
    """Scan harmonized h5ad files and produce summary table + publication plots."""
    import glob
    import matplotlib
    matplotlib.use("Agg")          # non-interactive backend — safe on servers
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="paper", font_scale=1.35)
    plt.rcParams["figure.dpi"]    = 300
    plt.rcParams["savefig.dpi"]   = 300
    plt.rcParams["savefig.bbox"]  = "tight"

    report_dir = os.path.join(scan_path, "reports", "atlas_summary")
    os.makedirs(report_dir, exist_ok=True)

    # ── Scan files
    h5ad_files = sorted(set(glob.glob(
        os.path.join(scan_path, "**", "*.h5ad"), recursive=True
    )))
    print(f"Scanning {len(h5ad_files)} h5ad files...")

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
        for ext in ("png", "pdf"):
            path = os.path.join(report_dir, f"{stem}.{ext}")
            fig.savefig(path)
        print(f"Saved figure: {os.path.join(report_dir, stem)}.png/pdf")

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
                 .head(15))
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
    p.add_argument("--csv",     default=None, help="Path to metadata CSV (overrides USER CONFIG)")
    p.add_argument("--gtf",     default=None, help="Path to GRCh38 GTF file (overrides USER CONFIG)")
    p.add_argument("--data",    default=None, help="Data root directory (overrides USER CONFIG)")
    p.add_argument("--output",  default=None, help="Output root directory (overrides USER CONFIG)")
    p.add_argument("--summary", action="store_true", help="Generate summary plots after pipeline")
    p.add_argument("--summary-only", action="store_true",
                   help="Skip pipeline, only generate summary plots from existing h5ad files")
    return p.parse_args()


def main():
    args = parse_args()

    # Resolve config — CLI args > env vars > USER CONFIG defaults
    data_root    = args.data   or DATA_ROOT
    gtf_file     = args.gtf    or GTF_FILE
    metadata_csv = args.csv    or METADATA_CSV
    output_root  = args.output or OUTPUT_ROOT

    # Change to data root so relative paths in CSV resolve correctly
    os.chdir(data_root)
    print(f"Working directory : {os.getcwd()}")
    print(f"GTF file          : {gtf_file}")
    print(f"Metadata CSV      : {metadata_csv}")
    print(f"Output root       : {output_root}\n")

    if not args.summary_only:
        failed = run_pipeline(
            metadata_csv=metadata_csv,
            gtf_file=gtf_file,
            output_root=output_root,
        )
        if failed:
            sys.exit(1)

    if args.summary or args.summary_only:
        print("\nGenerating summary plots...")
        run_summary(scan_path=output_root or ".")


if __name__ == "__main__":
    main()
