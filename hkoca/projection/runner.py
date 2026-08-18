"""Project query datasets onto the HKOCA atlas with scPoli surgery."""

from __future__ import annotations

import json
import logging
import random
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from hkoca.projection.config import load_projection_config, packaged_prototype_map_path
from hkoca.projection.plots import (
    composition_table,
    knn_project_onto_atlas_umap,
    plot_composition,
    plot_composition_by_sample,
    plot_confusion,
    plot_joint_umap_panels,
    plot_overlay,
    plot_query_categorical,
    plot_similarity,
    plot_uncertainty,
    subsample_rows,
)
from hkoca.projection.query_io import ensure_sample_id, load_query_adata
from hkoca.projection.scpoli import (
    align_to_reference_genes,
    classify_safe,
    ensure_query_condition_columns,
    ensure_raw_counts,
    get_latent_safe,
    load_atlas_obs,
    load_atlas_subsample_h5py,
    load_prototype_bridge,
    load_query_model,
    load_reference_genes,
    majority_map,
    patch_anndata_for_scarches,
    read_atlas_obsm,
    scpoli_placeholder_col,
    train_query_model,
    validate_model_dir,
)

from hkoca.projection.stack import require_projection_stack

logger = logging.getLogger("hkoca.projection")

REQUIRED_ATLAS_OBS = ("Level_3_Integrated", "Level_2_Integrated", "Level_1_Integrated")


def _setup_logging_file(log_path: Path) -> None:
    if any(isinstance(h, logging.FileHandler) and h.baseFilename == str(log_path) for h in logger.handlers):
        return
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"))
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)


def _assign_predicted_labels(
    query,
    *,
    preds_proto: np.ndarray,
    uncert: np.ndarray,
    prototype_to_l3: dict[str, str],
    l3_to_l2: dict[str, str],
    l3_to_l1: dict[str, str],
) -> None:
    proto = pd.Series(np.asarray(preds_proto).astype(str), index=query.obs_names)
    mapped = proto.map(prototype_to_l3)
    n_unmapped = int(mapped.isna().sum())
    if n_unmapped:
        missing = sorted(proto[mapped.isna()].unique().tolist())
        logger.warning(
            "%d cells have prototype labels missing from the bridge map: %s",
            n_unmapped,
            ", ".join(missing[:12]),
        )
        mapped = mapped.fillna(proto)
    query.obs["Level_3_prototype"] = pd.Categorical(proto)
    query.obs["Level_3_pred"] = pd.Categorical(mapped.astype(str))
    query.obs["Level_3_uncert"] = np.asarray(uncert, dtype=np.float32)
    query.obs["Level_2_pred"] = query.obs["Level_3_pred"].astype(str).map(l3_to_l2)
    query.obs["Level_1_pred"] = query.obs["Level_3_pred"].astype(str).map(l3_to_l1)
    n_l2 = int(query.obs["Level_2_pred"].isna().sum())
    n_l1 = int(query.obs["Level_1_pred"].isna().sum())
    if n_l2 or n_l1:
        logger.warning(
            "Rollup produced missing labels (Level_2_pred NA=%d, Level_1_pred NA=%d).",
            n_l2,
            n_l1,
        )


def _write_prediction_tables(query, tables_dir: Path, stem: str, batch_key: str) -> None:
    cols = [c for c in (batch_key, "study", "Level_3_prototype", "Level_3_pred", "Level_3_uncert", "Level_2_pred", "Level_1_pred") if c in query.obs.columns]
    pred_path = tables_dir / f"query_predictions_{stem}.tsv"
    query.obs[cols].to_csv(pred_path, sep="\t")
    logger.info("Wrote %s", pred_path)
    if "Level_2_pred" in query.obs.columns:
        comp = (
            query.obs.groupby(["Level_2_pred"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
            .sort_values("n_cells", ascending=False)
        )
        comp_path = tables_dir / f"query_Level2_composition_{stem}.tsv"
        comp.to_csv(comp_path, sep="\t", index=False)


def _encode_atlas_anchor_with_surgery(
    *,
    model,
    atlas_path: Path,
    atlas_obs: pd.DataFrame,
    ref_genes: list[str],
    placeholder_col: str,
    unknown_label: str,
    cfg: dict[str, Any],
):
    """Encode an atlas subsample with the query surgery model.

    Stored atlas X_emb is from the pre-surgery reference encoder. After
    surgery those coordinates are a different space than query X_scpoli.
    Re-encoding atlas cells with the surgery model (original atlas
    conditions) puts both on the same latent so kNN can transfer frozen
    X_umap_scpoli coordinates. Atlas UMAP is never recomputed.
    """
    rng = np.random.default_rng(cfg["seed"])
    n_atlas = atlas_obs.shape[0]
    n_take = min(int(cfg["atlas_umap_subsample"]), n_atlas)
    atlas_sub_idx = np.sort(rng.choice(n_atlas, size=n_take, replace=False))
    atlas_sub_barcodes = atlas_obs.index.to_numpy()[atlas_sub_idx]
    atlas_obs_sub = atlas_obs.iloc[atlas_sub_idx].copy()
    atlas_obs_sub["dataset_role"] = "atlas"
    logger.info(
        "Encoding atlas subsample with query surgery model (%s cells) for frozen UMAP overlay",
        f"{n_take:,}",
    )
    atlas_sub = load_atlas_subsample_h5py(
        atlas_path, atlas_sub_idx, atlas_obs_sub, atlas_sub_barcodes, ref_genes
    )
    atlas_sub.obs[placeholder_col] = unknown_label
    atlas_sub.obs[placeholder_col] = atlas_sub.obs[placeholder_col].astype("category")
    for cond in model.condition_keys_:
        if cond not in atlas_sub.obs.columns:
            raise KeyError(
                f"Atlas missing scPoli condition column {cond!r} needed to encode "
                "atlas anchors with the surgery model."
            )
        atlas_sub.obs[cond] = atlas_sub.obs[cond].astype(str).astype("category")
    atlas_sub.obsm["X_scpoli"] = np.asarray(
        get_latent_safe(model, atlas_sub, mean=True), dtype=np.float32
    )
    return atlas_sub, atlas_sub_idx


def _maybe_joint_umap(
    *,
    query,
    model,
    atlas_path: Path,
    atlas_obs: pd.DataFrame,
    ref_genes: list[str],
    placeholder_col: str,
    condition_key: str,
    celltype_key: str,
    unknown_label: str,
    cfg: dict[str, Any],
    projected_dir: Path,
    plots_dir: Path,
    stem: str,
    atlas_sub=None,
):
    import anndata as ad
    import scanpy as sc

    if atlas_sub is None:
        atlas_sub, _ = _encode_atlas_anchor_with_surgery(
            model=model,
            atlas_path=atlas_path,
            atlas_obs=atlas_obs,
            ref_genes=ref_genes,
            placeholder_col=placeholder_col,
            unknown_label=unknown_label,
            cfg=cfg,
        )
    lat = ad.AnnData(X=np.vstack([atlas_sub.obsm["X_scpoli"], query.obsm["X_scpoli"]]))
    lat.obs = pd.concat(
        [
            atlas_sub.obs[["dataset_role", condition_key, "Level_2_Integrated", celltype_key]].rename(
                columns={"Level_2_Integrated": "Level_2_plot", celltype_key: "Level_3_plot"}
            ),
            query.obs[["dataset_role", condition_key, "Level_2_pred", "Level_3_pred"]].rename(
                columns={"Level_2_pred": "Level_2_plot", "Level_3_pred": "Level_3_plot"}
            ),
        ],
        axis=0,
    )
    lat.obs["Level_2_plot"] = lat.obs["Level_2_plot"].astype(str)
    lat.obs["Level_3_plot"] = lat.obs["Level_3_plot"].astype(str)
    lat.obsm["X_scpoli"] = lat.X.copy()
    sc.pp.neighbors(lat, use_rep="X_scpoli", n_neighbors=15, random_state=cfg["seed"])
    sc.tl.umap(lat, random_state=cfg["seed"])
    q_mask = lat.obs["dataset_role"].to_numpy() == "query"
    query.obsm["X_umap_joint"] = np.asarray(lat.obsm["X_umap"][q_mask], dtype=np.float32)
    umap_path = projected_dir / f"joint_latent_umap_{stem}.h5ad"
    lat.write_h5ad(umap_path, compression="gzip")
    logger.info("Wrote exploratory joint UMAP: %s", umap_path)
    if cfg.get("save_plots", True):
        plot_joint_umap_panels(
            lat,
            query.obs["Level_3_uncert"].to_numpy(),
            plots_dir / f"joint_umap_{stem}.png",
            dpi=cfg["dpi"],
        )


def _atlas_umap_figures(
    *,
    query,
    atlas_path: Path,
    atlas_obs: pd.DataFrame,
    cfg: dict[str, Any],
    plots_dir: Path,
    tables_dir: Path,
    stem: str,
    batch_key: str,
    query_label_column: str | None,
    atlas_sub,
    atlas_sub_idx: np.ndarray,
) -> dict[str, Any]:
    umap_key = cfg["atlas_umap_key"]
    atlas_umap = read_atlas_obsm(atlas_path, umap_key)
    if atlas_umap is None:
        logger.warning(
            "Atlas missing %s; skipping atlas-UMAP overlay plots.",
            umap_key,
        )
        return {}

    atlas_lat = np.asarray(atlas_sub.obsm["X_scpoli"], dtype=np.float32)
    atlas_umap_anchor = np.asarray(atlas_umap[atlas_sub_idx], dtype=np.float32)
    q_lat = np.asarray(query.obsm["X_scpoli"], dtype=np.float32)
    if q_lat.shape[1] != atlas_lat.shape[1]:
        raise ValueError(
            f"Query scPoli latent dim {q_lat.shape[1]} does not match surgery-"
            f"encoded atlas latent dim {atlas_lat.shape[1]}."
        )
    logger.info(
        "Projecting query onto frozen atlas UMAP %s by kNN in surgery latent "
        "(query %s, atlas anchors %s)",
        umap_key,
        q_lat.shape,
        atlas_lat.shape,
    )
    q_umap, sim, d_ref = knn_project_onto_atlas_umap(
        atlas_lat, atlas_umap_anchor, q_lat, k=cfg["knn_neighbors"]
    )
    query.obsm["X_umap_scpoli_projected"] = q_umap
    query.obs["atlas_similarity"] = sim

    bg_idx = subsample_rows(atlas_umap.shape[0], cfg["atlas_bg_max"], cfg["seed"])
    bg_umap = atlas_umap[bg_idx]
    dpi = cfg["dpi"]

    plot_overlay(
        bg_umap,
        q_umap,
        query.obs["Level_3_pred"],
        plots_dir / f"umap_overlay_Level_3_pred_{stem}.png",
        title="Atlas + projected query (Level_3_pred)",
        dpi=dpi,
        palette_key="Level_3_Integrated",
    )
    plot_similarity(
        bg_umap,
        q_umap,
        sim,
        d_ref,
        plots_dir / f"query_on_atlas_umap_similarity_{stem}.png",
        dpi=dpi,
    )
    plot_query_categorical(
        bg_umap,
        q_umap,
        query.obs["Level_3_pred"].astype(str).to_numpy(),
        plots_dir / f"query_on_atlas_umap_Level3_pred_{stem}.png",
        title="Query projected on atlas UMAP · Level_3_pred",
        dpi=dpi,
        palette_key="Level_3_Integrated",
    )
    plot_query_categorical(
        bg_umap,
        q_umap,
        query.obs["Level_2_pred"].astype(str).to_numpy(),
        plots_dir / f"query_on_atlas_umap_Level2_pred_{stem}.png",
        title="Query projected on atlas UMAP · Level_2_pred",
        dpi=dpi,
        palette_key="Level_2_Integrated",
    )
    plot_query_categorical(
        bg_umap,
        q_umap,
        query.obs["Level_1_pred"].astype(str).to_numpy(),
        plots_dir / f"query_on_atlas_umap_Level1_pred_{stem}.png",
        title="Query projected on atlas UMAP · Level_1_pred",
        dpi=dpi,
        palette_key="Level_1_Integrated",
    )
    plot_uncertainty(
        bg_umap,
        q_umap,
        query.obs["Level_3_uncert"].to_numpy(dtype=float),
        plots_dir / f"query_on_atlas_umap_uncertainty_{stem}.png",
        dpi=dpi,
    )
    plot_composition_by_sample(
        query.obs,
        plots_dir / f"query_Level2_composition_by_sample_{stem}.png",
        dpi=dpi,
        batch_key=batch_key,
    )

    if "Level_3_Integrated" in atlas_obs.columns:
        ref_comp = composition_table(atlas_obs["Level_3_Integrated"], "atlas")
        query_comp = composition_table(query.obs["Level_3_pred"], "projected_query")
        pd.concat([ref_comp, query_comp], ignore_index=True).to_csv(
            tables_dir / f"composition_Level_3_Integrated_{stem}.csv", index=False
        )
        plot_composition(ref_comp, query_comp, plots_dir / f"composition_Level_3_Integrated_{stem}.png", dpi=dpi)

    if query_label_column and query_label_column in query.obs.columns:
        conf_df = plot_confusion(
            query.obs[query_label_column],
            query.obs["Level_3_pred"],
            plots_dir / f"confusion_{query_label_column}_vs_Level_3_pred_{stem}.png",
            dpi=dpi,
        )
        conf_df.to_csv(tables_dir / f"confusion_{query_label_column}_vs_Level_3_pred_{stem}.csv")

    sim_stats = {
        "n_query_cells": int(len(sim)),
        "n_atlas_anchor_cells": int(atlas_lat.shape[0]),
        "k_neighbors": int(cfg["knn_neighbors"]),
        "atlas_latent_key": "surgery_X_scpoli",
        "atlas_umap_key": umap_key,
        "atlas_latent_calib_dist_p95": float(d_ref),
        "mean_similarity": float(sim.mean()),
        "median_similarity": float(np.median(sim)),
        "pct_ge_50": float((sim >= 0.50).mean() * 100),
        "pct_ge_70": float((sim >= 0.70).mean() * 100),
        "pct_ge_85": float((sim >= 0.85).mean() * 100),
        "note": (
            "Query and an atlas subsample are encoded with the query surgery "
            "model. Query cells are placed on the fixed atlas UMAP "
            "(X_umap_scpoli) by kNN in that shared latent. Atlas UMAP is "
            "never recomputed. Level_3_uncert is scaled per query batch."
        ),
    }
    (tables_dir / f"atlas_query_similarity_summary_{stem}.json").write_text(
        json.dumps(sim_stats, indent=2), encoding="utf-8"
    )
    return sim_stats


def project_query(
    *,
    atlas_h5ad: str | Path,
    query_path: str | Path,
    model_dir: str | Path,
    output_dir: str | Path,
    query_label_column: str | None = None,
    query_batch_key: str | None = None,
    cfg: dict[str, Any] | None = None,
    force: bool = False,
    joint_umap: bool | None = None,
) -> Path:
    """Run scPoli surgery, prototype classify, and write labels + figures."""
    require_projection_stack()
    import torch

    atlas_path = Path(atlas_h5ad).expanduser().resolve()
    query_src = Path(query_path).expanduser().resolve()
    model_path = Path(model_dir).expanduser().resolve()
    out_root = Path(output_dir).expanduser().resolve()
    cfg = cfg or load_projection_config()
    if joint_umap is not None:
        cfg = dict(cfg)
        cfg["joint_umap"] = bool(joint_umap)

    batch_key = query_batch_key or cfg.get("query_batch_key") or "sample_id"
    condition_key = cfg.get("condition_key") or "sample_id"
    celltype_key = cfg.get("celltype_key") or "Level_3_Integrated"
    unknown_label = cfg.get("unknown_label") or "Unknown"
    counts_layer = cfg.get("counts_layer") or "counts"
    query_label_column = query_label_column or cfg.get("query_label_column")

    projected_dir = out_root / cfg["projected_subdir"]
    plots_dir = out_root / cfg["plots_subdir"]
    tables_dir = out_root / cfg["tables_subdir"]
    models_dir = out_root / cfg["models_subdir"] / "scpoli_query_surgery"
    converted_dir = out_root / "query_converted"
    out_root.mkdir(parents=True, exist_ok=True)
    for d in (projected_dir, plots_dir, tables_dir, models_dir, converted_dir):
        d.mkdir(parents=True, exist_ok=True)

    stem = query_src.stem
    out_h5ad = projected_dir / f"{stem}_projected.h5ad"
    meta_json = projected_dir / f".{stem}_projected.meta.json"
    _setup_logging_file(out_root / "projection.log")

    if not atlas_path.is_file():
        raise FileNotFoundError(f"Atlas h5ad not found: {atlas_path}")
    validate_model_dir(model_path)
    if not query_src.is_file():
        raise FileNotFoundError(f"Query file not found: {query_src}")

    if not force and cfg.get("skip_existing") and out_h5ad.is_file() and out_h5ad.stat().st_size > 0:
        logger.info("Skipping existing projected output: %s", out_h5ad)
        return out_h5ad

    if not torch.cuda.is_available():
        logger.warning("CUDA not available; scPoli surgery will run on CPU and may be very slow.")

    patch_anndata_for_scarches()
    random.seed(cfg["seed"])
    np.random.seed(cfg["seed"])

    logger.info("Loading query: %s", query_src)
    converted_h5ad = converted_dir / f"{stem}_rna_counts.h5ad"
    query_raw = load_query_adata(query_src, converted_h5ad=converted_h5ad, force=force)
    logger.info("Query: %s cells x %s genes", f"{query_raw.n_obs:,}", f"{query_raw.n_vars:,}")
    query_raw = ensure_sample_id(query_raw, fallback=stem, batch_key=batch_key)
    ensure_query_condition_columns(
        query_raw, model_path, batch_key=batch_key, study_fallback=stem
    )
    if condition_key != batch_key and condition_key not in query_raw.obs.columns:
        query_raw.obs[condition_key] = query_raw.obs[batch_key].astype(str)
    query_raw.obs[condition_key] = query_raw.obs[condition_key].astype(str).astype("category")

    ref_genes = load_reference_genes(model_path)
    logger.info("Reference model genes (var_names.csv): %s", f"{len(ref_genes):,}")
    placeholder_col = scpoli_placeholder_col(model_path)
    query = ensure_raw_counts(query_raw, prefer_layer=counts_layer)
    query = align_to_reference_genes(query, ref_genes)
    query.obs[placeholder_col] = unknown_label
    query.obs[placeholder_col] = query.obs[placeholder_col].astype("category")
    query.obs["dataset_role"] = "query"
    query.obs["conditions_combined"] = query.obs[condition_key].astype(str)

    logger.info("Loading atlas metadata: %s", atlas_path)
    atlas_obs = load_atlas_obs(atlas_path)
    for col in REQUIRED_ATLAS_OBS:
        if col not in atlas_obs.columns:
            raise KeyError(f"Atlas missing required obs column: {col!r}")
    l3_to_l2 = majority_map(atlas_obs, celltype_key, "Level_2_Integrated")
    l3_to_l1 = majority_map(atlas_obs, celltype_key, "Level_1_Integrated")
    map_df = pd.DataFrame(
        {
            celltype_key: list(l3_to_l2.keys()),
            "Level_2_Integrated": [l3_to_l2[k] for k in l3_to_l2],
            "Level_1_Integrated": [l3_to_l1[k] for k in l3_to_l2],
        }
    ).sort_values(celltype_key)
    map_df.to_csv(tables_dir / "atlas_Level3Integrated_to_Level2_Level1_majority_map.tsv", sep="\t", index=False)

    prototype_map = Path(cfg.get("prototype_map") or packaged_prototype_map_path())
    prototype_to_l3 = load_prototype_bridge(prototype_map, celltype_key)
    pd.read_csv(prototype_map, sep="\t").to_csv(
        tables_dir / "atlas_prototype_to_Level3Integrated_map.tsv", sep="\t", index=False
    )

    scpoli_query = load_query_model(query, model_path, unknown_label=unknown_label)
    logger.info(
        "Query model ready. condition_keys=%s input_dim=%s latent_dim=%s",
        scpoli_query.condition_keys_,
        scpoli_query.input_dim_,
        scpoli_query.latent_dim_,
    )
    train_query_model(
        scpoli_query,
        n_epochs=int(cfg["n_epochs"]),
        pretrain_epochs=int(cfg["pretrain_epochs"]),
        eta=float(cfg["eta"]),
    )
    scpoli_query.save(str(models_dir), overwrite=True)
    logger.info("Saved query surgery model: %s", models_dir)

    scpoli_query.model.eval()
    q_latent = np.asarray(get_latent_safe(scpoli_query, query, mean=True), dtype=np.float32)
    query.obsm["X_scpoli"] = q_latent
    results = classify_safe(scpoli_query, query, scale_uncertainties=True)
    proto_key = scpoli_query.cell_type_keys_[0]
    _assign_predicted_labels(
        query,
        preds_proto=results[proto_key]["preds"],
        uncert=results[proto_key]["uncert"],
        prototype_to_l3=prototype_to_l3,
        l3_to_l2=l3_to_l2,
        l3_to_l1=l3_to_l1,
    )
    logger.info("Level_3_pred counts:\n%s", query.obs["Level_3_pred"].value_counts().head(20).to_string())
    logger.info("Mean Level_3_uncert: %.4f", float(np.mean(query.obs["Level_3_uncert"])))

    sim_stats: dict[str, Any] = {}
    atlas_sub = None
    atlas_sub_idx = None
    if cfg.get("save_plots", True) or cfg.get("joint_umap"):
        atlas_sub, atlas_sub_idx = _encode_atlas_anchor_with_surgery(
            model=scpoli_query,
            atlas_path=atlas_path,
            atlas_obs=atlas_obs,
            ref_genes=ref_genes,
            placeholder_col=placeholder_col,
            unknown_label=unknown_label,
            cfg=cfg,
        )
    if cfg.get("save_plots", True):
        sim_stats = _atlas_umap_figures(
            query=query,
            atlas_path=atlas_path,
            atlas_obs=atlas_obs,
            cfg=cfg,
            plots_dir=plots_dir,
            tables_dir=tables_dir,
            stem=stem,
            batch_key=batch_key,
            query_label_column=query_label_column,
            atlas_sub=atlas_sub,
            atlas_sub_idx=atlas_sub_idx,
        )

    if cfg.get("joint_umap"):
        _maybe_joint_umap(
            query=query,
            model=scpoli_query,
            atlas_path=atlas_path,
            atlas_obs=atlas_obs,
            ref_genes=ref_genes,
            placeholder_col=placeholder_col,
            condition_key=condition_key,
            celltype_key=celltype_key,
            unknown_label=unknown_label,
            cfg=cfg,
            projected_dir=projected_dir,
            plots_dir=plots_dir,
            stem=stem,
            atlas_sub=atlas_sub,
        )

    query.uns["projection_atlas"] = str(atlas_path)
    query.uns["projection_query"] = str(query_src)
    query.uns["projection_model"] = str(model_path)
    query.uns["projection_method"] = "scpoli_surgery"

    summary = {
        "atlas_h5ad": str(atlas_path),
        "query": str(query_src),
        "model_dir": str(model_path),
        "projected_h5ad": str(out_h5ad),
        "n_query_cells": int(query.n_obs),
        "n_query_genes": int(query.n_vars),
        "n_sample_id": int(query.obs[condition_key].nunique()),
        "celltype_key": celltype_key,
        "label_transfer": "scpoli_prototype_classify",
        "cuda": bool(torch.cuda.is_available()),
        "mean_Level_3_uncert": float(np.mean(query.obs["Level_3_uncert"])),
        "similarity": sim_stats,
        "atlas_umap_frozen": True,
    }

    min_sim = float(cfg.get("min_mean_similarity") or 0.0)
    mean_sim = sim_stats.get("mean_similarity") if sim_stats else None
    if min_sim > 0 and mean_sim is not None and float(mean_sim) < min_sim:
        (tables_dir / "projection_summary.json").write_text(
            json.dumps(summary, indent=2), encoding="utf-8"
        )
        raise RuntimeError(
            f"Query did not land on the frozen atlas UMAP (mean similarity "
            f"{float(mean_sim):.4f} < {min_sim}). Overlay kNN uses surgery-"
            f"encoded atlas anchors mapped to {sim_stats.get('atlas_umap_key')}. "
            "Atlas X_umap_scpoli is never recomputed."
        )

    logger.info("Saving projected query: %s", out_h5ad)
    query.write_h5ad(out_h5ad, compression="gzip")
    _write_prediction_tables(query, tables_dir, stem, batch_key)
    (tables_dir / "projection_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    meta_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info("Projection completed. Outputs under: %s", out_root)
    return out_h5ad
