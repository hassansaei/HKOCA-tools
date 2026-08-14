"""Load projection query objects from h5ad or Seurat RDS (RNA counts)."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from hkoca.conda_env import ENV_INTEGRATION, r_env_for_rscript, resolve_env_prefix

logger = logging.getLogger("hkoca.projection")

DEFAULT_INTEGRATION_ENV = ENV_INTEGRATION
BATCH_ALIASES = ("sample_id", "sample", "donor_id", "batch", "orig.ident", "library_id")
RDS_SUFFIXES = {".rds", ".Rds", ".RDS"}


def rds_converter_script() -> Path:
    return Path(resources.files("hkoca.projection.r").joinpath("seurat_rds_to_h5ad.R")).resolve()


def find_rscript() -> str:
    override = os.environ.get("HKOCA_RSCRIPT", "").strip()
    if override:
        path = Path(override)
        if not path.is_file():
            raise FileNotFoundError(f"HKOCA_RSCRIPT is set but not a file: {override}")
        return str(path.resolve())

    env_name = os.environ.get("HKOCA_INTEGRATION_ENV", DEFAULT_INTEGRATION_ENV).strip() or DEFAULT_INTEGRATION_ENV
    prefix = resolve_env_prefix(env_name, "Rscript")
    if prefix is not None:
        rscript = prefix / "bin" / "Rscript"
        logger.info("Using integration conda env Rscript for RDS conversion: %s", rscript)
        return str(rscript)

    exe = shutil.which("Rscript")
    if exe is None:
        raise FileNotFoundError(
            "Rscript not found (needed to convert query .rds). Create conda env "
            f"'{env_name}' from conda/environment_integration.yaml, or set "
            "HKOCA_INTEGRATION_ENV / HKOCA_RSCRIPT."
        )
    logger.warning("Integration env '%s' not found; using Rscript on PATH (%s).", env_name, exe)
    return exe


def _subprocess_env_for_rscript(rscript: str) -> dict[str, str]:
    return r_env_for_rscript(rscript, harmonize_python=False)


def convert_seurat_rds_to_h5ad(
    rds_path: Path,
    output_h5ad: Path,
    *,
    force: bool = False,
) -> Path:
    """Export RNA counts from a Seurat RDS and write a counts h5ad."""
    import anndata as ad

    rds_path = Path(rds_path).expanduser().resolve()
    output_h5ad = Path(output_h5ad).expanduser().resolve()
    if output_h5ad.is_file() and output_h5ad.stat().st_size > 0 and not force:
        logger.info("Using existing converted query h5ad: %s", output_h5ad)
        return output_h5ad

    work = output_h5ad.parent / f".{output_h5ad.stem}_mtx"
    work.mkdir(parents=True, exist_ok=True)
    rscript = find_rscript()
    script = rds_converter_script()
    cmd = [
        rscript,
        str(script),
        "--input_rds",
        str(rds_path),
        "--output_dir",
        str(work),
    ]
    logger.info("Converting Seurat RDS to RNA-counts h5ad: %s", rds_path)
    proc = subprocess.run(
        cmd,
        check=False,
        env=_subprocess_env_for_rscript(rscript),
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        msg = (proc.stderr or proc.stdout or "").strip() or f"exit {proc.returncode}"
        raise RuntimeError(f"RDS conversion failed: {msg}")
    if proc.stdout:
        logger.info(proc.stdout.strip())

    mtx_path = work / "counts.mtx"
    var_path = work / "var.tsv"
    obs_path = work / "obs.csv"
    for p in (mtx_path, var_path, obs_path):
        if not p.is_file():
            raise FileNotFoundError(f"RDS converter did not write {p}")

    from scipy.io import mmread

    genes_x_cells = sp.csr_matrix(mmread(mtx_path))
    x = genes_x_cells.T.tocsr().astype(np.float32)
    var = pd.read_csv(var_path, sep="\t")
    gene_col = "gene" if "gene" in var.columns else var.columns[0]
    var = var.set_index(gene_col)
    var.index = var.index.astype(str)
    obs = pd.read_csv(obs_path)
    if "barcode" not in obs.columns:
        raise ValueError("RDS converter obs.csv is missing a barcode column.")
    obs["barcode"] = obs["barcode"].astype(str)
    obs = obs.set_index("barcode")
    if x.shape[0] != obs.shape[0] or x.shape[1] != var.shape[0]:
        raise ValueError(
            f"Converted matrix shape {x.shape} does not match obs={obs.shape[0]} var={var.shape[0]}."
        )
    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad, compression="gzip")
    logger.info("Wrote converted query h5ad: %s (%s cells x %s genes)", output_h5ad, f"{adata.n_obs:,}", f"{adata.n_vars:,}")
    return output_h5ad


def ensure_gene_names(adata) -> None:
    if "features" in adata.var.columns:
        adata.var_names = adata.var["features"].astype(str)
    elif "gene_symbols" in adata.var.columns:
        adata.var_names = adata.var["gene_symbols"].astype(str)
    adata.var_names_make_unique()


def ensure_sample_id(adata, *, fallback: str, batch_key: str = "sample_id"):
    if batch_key in adata.obs.columns:
        adata.obs[batch_key] = adata.obs[batch_key].astype(str)
        return adata
    for alt in BATCH_ALIASES:
        if alt == batch_key:
            continue
        if alt in adata.obs.columns:
            adata.obs[batch_key] = adata.obs[alt].astype(str)
            logger.warning("Created obs[%r] from obs[%r]", batch_key, alt)
            return adata
    adata.obs[batch_key] = str(fallback)
    logger.warning("No batch column found; set all cells to %s=%r", batch_key, fallback)
    return adata


def load_query_adata(
    query_path: Path,
    *,
    converted_h5ad: Path | None = None,
    force: bool = False,
):
    """Load query as AnnData from .h5ad or Seurat .rds (RNA counts)."""
    import scanpy as sc

    query_path = Path(query_path).expanduser().resolve()
    if not query_path.is_file():
        raise FileNotFoundError(f"Query file not found: {query_path}")

    suffix = query_path.suffix
    if suffix in RDS_SUFFIXES:
        out = converted_h5ad or (query_path.parent / f"{query_path.stem}_rna_counts.h5ad")
        h5ad_path = convert_seurat_rds_to_h5ad(query_path, out, force=force)
        adata = sc.read_h5ad(h5ad_path)
    elif suffix.lower() == ".h5ad":
        adata = sc.read_h5ad(query_path)
    else:
        raise ValueError(
            f"Unsupported query file type: {query_path.suffix}. Use .h5ad or Seurat .rds."
        )
    ensure_gene_names(adata)
    return adata
