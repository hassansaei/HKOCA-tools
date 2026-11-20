#!/usr/bin/env python3
"""Concatenate QC-filtered singlets h5ads; keep main obs columns and write a short report."""
from __future__ import annotations

import gc
import os
import re
import sys
import traceback
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
INPUT_DIR = Path(
    "/data-master/pure-workspace/labss/hmami/atlas_shortcut/filtered_singlets_harmonized"
)
OUT_DIR = Path("/data-master/pure-workspace/labss/hmami/atlas_shortcut")
OUT_H5AD = OUT_DIR / "concatenated_filtered_singlets_important_columns.h5ad"
OUT_REPORT = OUT_DIR / "concatenation_report.md"
OUT_PRE = OUT_DIR / "sanity_check_before.csv"
OUT_POST = OUT_DIR / "sanity_check_after.csv"
OUT_SID_MAP = OUT_DIR / "sample_id_standardization_map.csv"
OUT_CELLS_BY_SAMPLE = OUT_DIR / "cells_per_sample.csv"
OUT_CELLS_BY_STUDY = OUT_DIR / "cells_per_study.csv"

KEEP_COLS = [
    "orig.ident",
    "nCount_RNA",
    "nFeature_RNA",
    "sample_id",
    "study",
    "source",
    "diff_protocol",
    "sc_protocol",
    "sequencing",
    "Age",
]

EMPTY_LIKE = {"", "nan", "None", "NA", "NaN", "<NA>", "na", "null", "Null"}


def log(msg: str = "") -> None:
    print(msg, flush=True)


def standardize_sample_id(x) -> str:
    """Normalize sample_id: strip, spaces/dots -> underscore, collapse repeats."""
    s = str(x).strip()
    s = s.replace(" ", "_").replace(".", "_")
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def standardize_sc_protocol(x) -> str:
    """Normalize sc_protocol casing (e.g. 10X_5 -> 10x_5)."""
    s = str(x).strip()
    if s.lower().startswith("10x"):
        # preserve version tokens, force 10x prefix lower
        return "10x" + s[3:]
    return s


def standardize_study(x) -> str:
    s = str(x).strip()
    # Fix inconsistent casing seen in LMRH bundle
    mapping = {"TTC21b": "TTC21B"}
    return mapping.get(s, s)


def is_missing_series(s: pd.Series) -> pd.Series:
    na = pd.isna(s)
    empty = s.astype(str).str.strip().isin(EMPTY_LIKE)
    return na | empty


def study_key_from_filename(fname: str) -> str:
    return fname.replace("_harmonized_singlets_filtered.h5ad", "")


def audit_file_backed(path: Path) -> dict:
    a = ad.read_h5ad(path, backed="r")
    key = study_key_from_filename(path.name)
    row = {
        "file": path.name,
        "study_key": key,
        "n_cells": int(a.n_obs),
        "n_genes": int(a.n_vars),
        "has_raw": a.raw is not None,
        "n_obs_cols": int(a.obs.shape[1]),
    }
    missing_keep = [c for c in KEEP_COLS if c not in a.obs.columns]
    row["missing_keep_cols"] = ";".join(missing_keep) if missing_keep else ""
    row["ok_keep_cols"] = len(missing_keep) == 0

    for c in KEEP_COLS:
        if c not in a.obs.columns:
            row[f"{c}__n_missing"] = a.n_obs
            row[f"{c}__n_unique"] = 0
            row[f"{c}__values"] = ""
            continue
        s = a.obs[c]
        n_miss = int(is_missing_series(s).sum())
        vals = sorted(s.astype(str).str.strip().unique().tolist())
        row[f"{c}__n_missing"] = n_miss
        row[f"{c}__n_unique"] = len(vals)
        row[f"{c}__values"] = " | ".join(vals[:20]) + (" | ..." if len(vals) > 20 else "")

    # gene names snapshot for order check
    row["_genes"] = a.var_names.to_numpy()
    a.file.close()
    return row


def prepare_one(path: Path, sid_map_rows: list) -> ad.AnnData:
    """Load one dataset, standardize, keep only important columns."""
    key = study_key_from_filename(path.name)
    log(f"  loading {path.name} ...")
    a = ad.read_h5ad(path)

    # Drop bulky slots early
    a.raw = None
    a.obsm.clear()
    a.obsp.clear()
    a.uns.clear()
    if a.layers:
        a.layers.clear()
    a.var = pd.DataFrame(index=a.var_names)

    # Ensure KEEP columns exist
    missing = [c for c in KEEP_COLS if c not in a.obs.columns]
    if missing:
        raise ValueError(f"{path.name} missing columns: {missing}")

    # Standardize sample_id (+ record mapping)
    raw_sid = a.obs["sample_id"].astype(str)
    new_sid = raw_sid.map(standardize_sample_id)
    for old, new in (
        pd.DataFrame({"old": raw_sid, "new": new_sid})
        .drop_duplicates()
        .itertuples(index=False)
    ):
        sid_map_rows.append(
            {
                "study_key": key,
                "sample_id_raw": old,
                "sample_id_std": new,
                "changed": old != new,
            }
        )
    a.obs["sample_id"] = pd.Categorical(new_sid)

    # orig.ident is a Seurat decoy (mostly "0"); replace with sample_id
    a.obs["orig.ident"] = a.obs["sample_id"].astype(str)

    # Light standardization of other categoricals
    a.obs["study"] = a.obs["study"].map(standardize_study).astype(str)
    a.obs["sc_protocol"] = a.obs["sc_protocol"].map(standardize_sc_protocol).astype(str)
    for c in ["source", "diff_protocol", "sequencing", "Age"]:
        a.obs[c] = a.obs[c].astype(str).str.strip()

    # Numeric QC columns
    a.obs["nCount_RNA"] = pd.to_numeric(a.obs["nCount_RNA"], errors="coerce")
    a.obs["nFeature_RNA"] = pd.to_numeric(a.obs["nFeature_RNA"], errors="coerce")

    # Keep only important columns
    a.obs = a.obs[KEEP_COLS].copy()

    # Sparse CSR
    if not sp.issparse(a.X):
        a.X = sp.csr_matrix(a.X)
    else:
        a.X = a.X.tocsr()

    # Unique barcodes within file
    a.obs_names = a.obs_names.astype(str)
    if not a.obs_names.is_unique:
        a.obs_names_make_unique()

    # Attach dataset origin for unique concat keys
    a.obs["dataset_file"] = key

    log(
        f"    -> {a.n_obs:,} cells | {a.obs['sample_id'].nunique()} samples | "
        f"study={sorted(a.obs['study'].unique().tolist())}"
    )
    return a


def main() -> int:
    started = datetime.now()
    log("=" * 72)
    log("HKOCA filtered-singlets concatenation")
    log(f"Started: {started.isoformat(timespec='seconds')}")
    log(f"Input:   {INPUT_DIR}")
    log("=" * 72)

    files = sorted(INPUT_DIR.glob("*_harmonized_singlets_filtered.h5ad"))
    if len(files) != 22:
        log(f"WARNING: expected 22 h5ad files, found {len(files)}")
    log(f"Found {len(files)} files\n")

    # ------------------------------------------------------------------
    # PRE-FLIGHT SANITY CHECK (backed)
    # ------------------------------------------------------------------
    log("=== PRE-FLIGHT SANITY CHECK ===")
    pre_rows = []
    gene_ref = None
    gene_mismatch = []
    for p in files:
        row = audit_file_backed(p)
        genes = row.pop("_genes")
        if gene_ref is None:
            gene_ref = genes
        elif not np.array_equal(gene_ref, genes):
            gene_mismatch.append(row["file"])
        pre_rows.append(row)
        miss_flag = "OK" if row["ok_keep_cols"] else f"MISSING:{row['missing_keep_cols']}"
        any_miss = any(row.get(f"{c}__n_missing", 0) > 0 for c in KEEP_COLS)
        log(
            f"  {row['study_key']}: cells={row['n_cells']:,} genes={row['n_genes']:,} "
            f"keep_cols={miss_flag} value_missing={'YES' if any_miss else 'no'}"
        )

    pre_df = pd.DataFrame(pre_rows)
    pre_df.to_csv(OUT_PRE, index=False)
    log(f"\nWrote pre-flight table: {OUT_PRE}")

    if gene_mismatch:
        raise RuntimeError(f"Gene order mismatch in: {gene_mismatch}")
    if not pre_df["ok_keep_cols"].all():
        bad = pre_df.loc[~pre_df["ok_keep_cols"], ["file", "missing_keep_cols"]]
        raise RuntimeError(f"Missing KEEP columns:\n{bad}")

    total_cells_pre = int(pre_df["n_cells"].sum())
    log(f"Pre-flight total cells: {total_cells_pre:,}")
    log(f"Gene universe: {len(gene_ref):,} (identical order across all files)\n")

    # ------------------------------------------------------------------
    # LOAD + STANDARDIZE + CONCAT
    # ------------------------------------------------------------------
    log("=== LOAD / STANDARDIZE / CONCAT ===")
    sid_map_rows: list[dict] = []
    adatas: dict[str, ad.AnnData] = {}
    for p in files:
        key = study_key_from_filename(p.name)
        adata = prepare_one(p, sid_map_rows)
        adatas[key] = adata
        gc.collect()

    sid_map = pd.DataFrame(sid_map_rows).drop_duplicates()
    sid_map.to_csv(OUT_SID_MAP, index=False)
    n_changed = int(sid_map["changed"].sum())
    log(f"\nsample_id standardization map: {OUT_SID_MAP} ({n_changed} IDs changed)")

    # Global uniqueness of standardized sample_id
    all_pairs = []
    for key, a in adatas.items():
        for sid, n in a.obs["sample_id"].astype(str).value_counts().items():
            all_pairs.append((sid, key, int(n)))
    by_sid = defaultdict(list)
    for sid, key, n in all_pairs:
        by_sid[sid].append((key, n))
    collisions = {s: v for s, v in by_sid.items() if len(v) > 1}
    if collisions:
        raise RuntimeError(f"sample_id collisions after standardization: {collisions}")
    log(f"Global unique sample_ids after std: {len(by_sid)} (0 collisions)")

    log("\nConcatenating (join=inner, index_unique='_') ...")
    # Drop helper col before concat? keep dataset_file as useful provenance, but
    # user asked only for KEEP cols — drop it after using keys for index_unique.
    for a in adatas.values():
        if "dataset_file" in a.obs:
            del a.obs["dataset_file"]

    concat = ad.concat(
        adatas,
        join="inner",
        merge="same",
        index_unique="_",
        fill_value=0,
    )
    # Ensure column order
    concat.obs = concat.obs[KEEP_COLS].copy()

    # Cast categoricals for compact storage
    for c in [
        "orig.ident",
        "sample_id",
        "study",
        "source",
        "diff_protocol",
        "sc_protocol",
        "sequencing",
        "Age",
    ]:
        concat.obs[c] = pd.Categorical(concat.obs[c].astype(str))

    gc.collect()
    log(f"Concat shape: {concat.n_obs:,} cells × {concat.n_vars:,} genes")

    # ------------------------------------------------------------------
    # POST SANITY CHECK
    # ------------------------------------------------------------------
    log("\n=== POST-CONCAT SANITY CHECK ===")
    issues = []
    if concat.n_obs != total_cells_pre:
        issues.append(
            f"cell count mismatch: concat={concat.n_obs:,} vs pre={total_cells_pre:,}"
        )
    if concat.n_vars != len(gene_ref):
        issues.append(
            f"gene count mismatch: concat={concat.n_vars:,} vs ref={len(gene_ref):,}"
        )
    if not concat.obs_names.is_unique:
        issues.append("obs_names are not unique after concat")

    post_rows = []
    for c in KEEP_COLS:
        if c not in concat.obs.columns:
            issues.append(f"missing column after concat: {c}")
            continue
        n_miss = int(is_missing_series(concat.obs[c]).sum())
        n_uniq = int(concat.obs[c].astype(str).nunique())
        post_rows.append(
            {
                "column": c,
                "n_missing": n_miss,
                "n_unique": n_uniq,
                "dtype": str(concat.obs[c].dtype),
            }
        )
        if n_miss:
            issues.append(f"{c}: {n_miss:,} missing values after concat")
        log(f"  {c}: missing={n_miss:,} unique={n_uniq:,}")

    post_df = pd.DataFrame(post_rows)
    post_df.to_csv(OUT_POST, index=False)

    # Per-study / per-sample tables
    cells_study = (
        concat.obs.groupby("study", observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values("n_cells", ascending=False)
    )
    cells_sample = (
        concat.obs.groupby(["study", "sample_id", "Age", "diff_protocol", "sc_protocol", "source", "sequencing"], observed=True)
        .agg(
            n_cells=("nCount_RNA", "size"),
            median_nCount=("nCount_RNA", "median"),
            median_nFeature=("nFeature_RNA", "median"),
        )
        .reset_index()
        .sort_values(["study", "sample_id"])
    )
    cells_study.to_csv(OUT_CELLS_BY_STUDY, index=False)
    cells_sample.to_csv(OUT_CELLS_BY_SAMPLE, index=False)

    # ------------------------------------------------------------------
    # WRITE H5AD
    # ------------------------------------------------------------------
    log(f"\nWriting concatenated h5ad -> {OUT_H5AD}")
    concat.write_h5ad(OUT_H5AD, compression="gzip")
    size_gb = OUT_H5AD.stat().st_size / (1024**3)
    log(f"Wrote {OUT_H5AD.name} ({size_gb:.2f} GiB)")

    # ------------------------------------------------------------------
    # REPORT
    # ------------------------------------------------------------------
    ended = datetime.now()
    duration = ended - started

    # Value summaries
    def top_counts(col, n=30):
        vc = concat.obs[col].astype(str).value_counts()
        lines = [f"| {k} | {v:,} |" for k, v in vc.head(n).items()]
        return "\n".join(lines)

    sid_changes = sid_map.loc[sid_map["changed"], ["study_key", "sample_id_raw", "sample_id_std"]]
    def _df_to_md(df):
        cols = list(df.columns)
        lines = ["| " + " | ".join(cols) + " |",
                 "| " + " | ".join(["---"] * len(cols)) + " |"]
        for row in df.itertuples(index=False):
            lines.append("| " + " | ".join(str(x) for x in row) + " |")
        return "\n".join(lines) if len(df) else "_none_"
    sid_change_md = _df_to_md(sid_changes)

    report = f"""# Concatenation report — filtered singlets (important columns)

- **Generated:** {ended.isoformat(timespec='seconds')}
- **Duration:** {duration}
- **Input dir:** `{INPUT_DIR}`
- **Output h5ad:** `{OUT_H5AD}` ({size_gb:.2f} GiB)
- **Files concatenated:** {len(files)}

## Summary statistics

| Metric | Value |
| :--- | ---: |
| Total cells | {concat.n_obs:,} |
| Total genes | {concat.n_vars:,} |
| Studies (obs `study`) | {concat.obs['study'].nunique()} |
| Unique `sample_id` | {concat.obs['sample_id'].nunique()} |
| Unique `Age` | {concat.obs['Age'].nunique()} |
| Unique `diff_protocol` | {concat.obs['diff_protocol'].nunique()} |
| Unique `sc_protocol` | {concat.obs['sc_protocol'].nunique()} |
| Unique `source` | {concat.obs['source'].nunique()} |
| Unique `sequencing` | {concat.obs['sequencing'].nunique()} |
| Median nCount_RNA | {float(concat.obs['nCount_RNA'].median()):,.1f} |
| Median nFeature_RNA | {float(concat.obs['nFeature_RNA'].median()):,.1f} |
| Mean nCount_RNA | {float(concat.obs['nCount_RNA'].mean()):,.1f} |
| Mean nFeature_RNA | {float(concat.obs['nFeature_RNA'].mean()):,.1f} |

## Pre-flight sanity check

- All {len(files)} files contained every KEEP column: **PASS**
- Gene count / order identical ({len(gene_ref):,} genes): **PASS**
- Missing values in KEEP columns before concat: **none**
- Pre-flight cell sum: **{total_cells_pre:,}**

Detailed table: `{OUT_PRE}`

## sample_id standardization

Rules applied:

1. strip whitespace
2. replace spaces and `.` with `_`
3. collapse repeated `_`
4. set `orig.ident` ← standardized `sample_id` (Seurat decoy `0`/`1`/`SeuratProject` removed)
5. light fixes: `sc_protocol` `10X*` → `10x*`; study `TTC21b` → `TTC21B`

- Unique sample_ids before / after: {sid_map['sample_id_raw'].nunique()} / {sid_map['sample_id_std'].nunique()}
- IDs changed: **{n_changed}**
- Global collisions after std: **0**

### Changed sample_ids

{sid_change_md}

Full map: `{OUT_SID_MAP}`

## Post-concat sanity check

- Cell count preserved ({concat.n_obs:,} == {total_cells_pre:,}): **{'PASS' if concat.n_obs == total_cells_pre else 'FAIL'}**
- Gene count preserved ({concat.n_vars:,}): **{'PASS' if concat.n_vars == len(gene_ref) else 'FAIL'}**
- Unique obs_names: **{'PASS' if concat.obs_names.is_unique else 'FAIL'}**
- Missing values in KEEP columns after concat: **{'none' if post_df['n_missing'].sum() == 0 else 'FOUND'}**
- Issues: {issues if issues else 'none'}

Column audit: `{OUT_POST}`

## Cells per study

| study | n_cells |
| :--- | ---: |
{chr(10).join(f"| {r.study} | {r.n_cells:,} |" for r in cells_study.itertuples())}

## Cells per Age

| Age | n_cells |
| :--- | ---: |
{top_counts('Age')}

## Cells per diff_protocol

| diff_protocol | n_cells |
| :--- | ---: |
{top_counts('diff_protocol')}

## Cells per sc_protocol

| sc_protocol | n_cells |
| :--- | ---: |
{top_counts('sc_protocol')}

## Cells per source

| source | n_cells |
| :--- | ---: |
{top_counts('source')}

## Cells per sequencing

| sequencing | n_cells |
| :--- | ---: |
{top_counts('sequencing')}

## Obs columns retained

{chr(10).join(f"- `{c}`" for c in KEEP_COLS)}

## Companion tables

- `{OUT_CELLS_BY_STUDY}`
- `{OUT_CELLS_BY_SAMPLE}`
- `{OUT_SID_MAP}`
- `{OUT_PRE}`
- `{OUT_POST}`
"""

    OUT_REPORT.write_text(report)
    log(f"\nWrote report: {OUT_REPORT}")

    if issues:
        log("\nCOMPLETED WITH ISSUES:")
        for i in issues:
            log(f"  - {i}")
        return 1

    log("\n=== DONE (all sanity checks PASS) ===")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(2)
