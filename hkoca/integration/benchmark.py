"""scIB metrics for Harmony / RPCA / CCA after hkoca integration run."""

from __future__ import annotations

import argparse
import logging
import sys
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger("hkoca.integration.benchmark")

DEFAULT_METHODS = ("unintegrated", "harmony", "rpca", "cca")


def packaged_benchmark_config() -> Path:
    return Path(resources.files("hkoca.integration").joinpath("benchmark_config.yaml")).resolve()


def load_benchmark_config(path: str | Path | None = None) -> dict[str, Any]:
    cfg_path = Path(path) if path else packaged_benchmark_config()
    with cfg_path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def benchmark_dir(output_dir: str | Path) -> Path:
    return Path(output_dir) / "benchmark"


def metrics_csv_path(output_dir: str | Path, cfg: dict[str, Any] | None = None) -> Path:
    metrics = (cfg or {}).get("metrics") or {}
    name = metrics.get("output_file") or "scib_raw_metrics.csv"
    return benchmark_dir(output_dir) / name


def benchmark_complete(output_dir: str | Path, cfg: dict[str, Any] | None = None) -> bool:
    path = metrics_csv_path(output_dir, cfg)
    return path.is_file() and path.stat().st_size > 0


def _resolve_label_key(obs: pd.DataFrame, cfg: dict[str, Any]) -> str:
    meta = cfg.get("metadata") or {}
    preferred = str(meta.get("label_key") or "Level_3")
    fallbacks = [preferred, *list(meta.get("label_key_fallbacks") or [])]
    for key in fallbacks:
        if key in obs.columns and obs[key].notna().any():
            n = obs[key].astype(str).replace({"": np.nan, "nan": np.nan, "None": np.nan}).notna().sum()
            if n > 0:
                return key
    raise ValueError(
        "No cell-type label column found for scIB. Expected Level_3 (from --annotated-h5ad) "
        f"in {list(obs.columns)}."
    )


def _overall_score(row: pd.Series, bio_keys: list[str], batch_keys: list[str],
                   bio_weight: float, batch_weight: float) -> float:
    bio = np.nanmean([row.get(k, np.nan) for k in bio_keys])
    batch = np.nanmean([row.get(k, np.nan) for k in batch_keys])
    return float(bio_weight * bio + batch_weight * batch)


def _load_method_adata(bench: Path, method: str, meta: pd.DataFrame, embed_key: str):
    import anndata as ad

    csv_path = bench / f"{method}_embeddings.csv"
    if not csv_path.is_file() or csv_path.stat().st_size <= 0:
        return None
    emb = pd.read_csv(csv_path, index_col=0)
    shared = meta.index.intersection(emb.index)
    if len(shared) < 10:
        logger.warning("[%s] Too few overlapping cells (%d); skipping.", method, len(shared))
        return None
    emb = emb.reindex(shared)
    obs = meta.loc[shared].copy()
    adata = ad.AnnData(obs=obs)
    adata.obsm[embed_key] = np.asarray(emb.values, dtype=np.float32)
    return adata


def _plot_results(df: pd.DataFrame, cfg: dict[str, Any], out_dir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    metrics_cfg = cfg["metrics"]
    score_cols = [c for c in metrics_cfg["bio_keys"] + metrics_cfg["batch_keys"] + ["overall"] if c in df.columns]
    plot_df = df.loc[:, score_cols].astype(float)

    heat_path = out_dir / metrics_cfg.get("heatmap_plot", "scib_metrics_heatmap.png")
    fig, ax = plt.subplots(figsize=(max(8, 0.9 * len(score_cols) + 3), max(3.5, 0.7 * len(plot_df) + 2)))
    sns.heatmap(plot_df, annot=True, fmt=".3f", cmap="YlGnBu", vmin=0, vmax=1, ax=ax)
    ax.set_title("scIB metrics (Harmony / RPCA / CCA)")
    fig.tight_layout()
    fig.savefig(heat_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved heatmap: %s", heat_path)

    rank_path = out_dir / metrics_cfg.get("ranking_plot", "scib_method_ranking.png")
    order = plot_df["overall"].sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = ["#0072B2" if m != "unintegrated" else "#BBBBBB" for m in order.index]
    ax.bar(order.index.astype(str), order.values, color=colors)
    ax.set_ylabel("Overall score (0.6 bio + 0.4 batch)")
    ax.set_ylim(0, 1)
    ax.set_title("Integration ranking")
    fig.tight_layout()
    fig.savefig(rank_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved ranking plot: %s", rank_path)


def run_scib_benchmark(
    output_dir: str | Path,
    *,
    config_path: str | Path | None = None,
    methods: list[str] | None = None,
) -> Path:
    try:
        import scanpy as sc
        import scib
    except ImportError as exc:
        raise ImportError(
            "scIB benchmark requires scanpy and scib. Install in the annotation env: "
            "pip install scib   (or update conda/environment_harmonize.yaml)"
        ) from exc

    cfg = load_benchmark_config(config_path)
    metrics_cfg = cfg["metrics"]
    int_cfg = cfg["integration"]
    meta_cfg = cfg["metadata"]
    bench = benchmark_dir(output_dir)
    meta_csv = bench / "cell_metadata.csv"
    if not meta_csv.is_file():
        raise FileNotFoundError(f"Missing benchmark metadata: {meta_csv}")

    meta = pd.read_csv(meta_csv, index_col=0)
    batch_key = str(meta_cfg.get("batch_key") or "sample_id")
    if batch_key not in meta.columns:
        if "orig.ident" in meta.columns:
            batch_key = "orig.ident"
        else:
            raise ValueError(f"Batch key '{batch_key}' not in metadata.")
    label_key = _resolve_label_key(meta, cfg)
    logger.info("scIB batch_key=%s label_key=%s cells=%d", batch_key, label_key, len(meta))

    method_list = list(methods) if methods else list(metrics_cfg.get("methods") or DEFAULT_METHODS)
    # Always score unintegrated PCA when the file exists, plus the three integration methods.
    wanted = []
    for name in ("unintegrated", *method_list):
        if name not in wanted:
            wanted.append(name)
    wanted = [m for m in wanted if m in ("unintegrated", "harmony", "rpca", "cca")]

    embed_key = str(metrics_cfg.get("embed_key") or "X_emb")
    cluster_key = str(metrics_cfg.get("cluster_key") or "leiden_scib")
    bio_keys = list(metrics_cfg["bio_keys"])
    batch_keys = list(metrics_cfg["batch_keys"])
    n_neighbors = int(int_cfg.get("n_neighbors") or 15)
    random_state = int(int_cfg.get("random_state") or 42)

    adata_pre = _load_method_adata(bench, "unintegrated", meta, embed_key)
    rows: list[dict[str, Any]] = []

    for method in wanted:
        adata = _load_method_adata(bench, method, meta, embed_key)
        if adata is None:
            logger.warning("No embeddings for method '%s'; skipping.", method)
            continue
        logger.info("Computing scIB metrics: %s", method)
        sc.pp.neighbors(adata, use_rep=embed_key, n_neighbors=n_neighbors, random_state=random_state)
        scib.metrics.cluster_optimal_resolution(
            adata,
            label_key=label_key,
            cluster_key=cluster_key,
            cluster_function=sc.tl.leiden,
            metric=scib.metrics.nmi,
        )
        row: dict[str, Any] = {"method": method}
        row["ARI"] = scib.me.ari(adata, cluster_key=cluster_key, label_key=label_key)
        row["NMI"] = scib.me.nmi(adata, cluster_key=cluster_key, label_key=label_key)
        row["ASW_celltype"] = scib.me.silhouette(adata, label_key=label_key, embed=embed_key)
        try:
            row["iso_label_F1"] = scib.me.isolated_labels_f1(
                adata, label_key=label_key, batch_key=batch_key, embed=embed_key, verbose=False
            )
        except Exception as exc:
            logger.warning("[%s] iso_label_F1 failed: %s", method, exc)
            row["iso_label_F1"] = float("nan")
        try:
            row["iso_label_ASW"] = scib.me.isolated_labels_asw(
                adata, label_key=label_key, batch_key=batch_key, embed=embed_key, verbose=False
            )
        except Exception as exc:
            logger.warning("[%s] iso_label_ASW failed: %s", method, exc)
            row["iso_label_ASW"] = float("nan")
        try:
            row["cLISI"] = scib.me.clisi_graph(
                adata, label_key=label_key, type_="embed", use_rep=embed_key,
                k0=int(metrics_cfg.get("lisi_k0") or 90),
                subsample=int(metrics_cfg.get("lisi_subsample") or 100),
                scale=bool(metrics_cfg.get("lisi_scale", True)),
                verbose=False,
            )
        except Exception as exc:
            logger.warning("[%s] cLISI failed: %s", method, exc)
            row["cLISI"] = float("nan")
        try:
            row["ASW_batch"] = scib.me.silhouette_batch(
                adata, batch_key=batch_key, label_key=label_key, embed=embed_key, verbose=False
            )
        except Exception as exc:
            logger.warning("[%s] ASW_batch failed: %s", method, exc)
            row["ASW_batch"] = float("nan")
        try:
            row["iLISI"] = scib.me.ilisi_graph(
                adata, batch_key=batch_key, type_="embed", use_rep=embed_key,
                k0=int(metrics_cfg.get("lisi_k0") or 90),
                subsample=int(metrics_cfg.get("lisi_subsample") or 100),
                scale=bool(metrics_cfg.get("lisi_scale", True)),
                verbose=False,
            )
        except Exception as exc:
            logger.warning("[%s] iLISI failed: %s", method, exc)
            row["iLISI"] = float("nan")
        try:
            row["graph_connectivity"] = scib.me.graph_connectivity(adata, label_key=label_key)
        except Exception as exc:
            logger.warning("[%s] graph_connectivity failed: %s", method, exc)
            row["graph_connectivity"] = float("nan")
        try:
            row["kBET"] = scib.me.kBET(
                adata, batch_key=batch_key, label_key=label_key, type_="embed", embed=embed_key
            )
        except Exception as exc:
            logger.warning("[%s] kBET failed: %s", method, exc)
            row["kBET"] = float("nan")
        try:
            if adata_pre is not None:
                row["PCR"] = scib.me.pcr_comparison(
                    adata_pre, adata, covariate=batch_key, embed=embed_key,
                    n_comps=int(metrics_cfg.get("pcr_n_comps") or 8),
                    scale=bool(metrics_cfg.get("pcr_scale", True)),
                    verbose=False,
                )
            else:
                row["PCR"] = float("nan")
        except Exception as exc:
            logger.warning("[%s] PCR failed: %s", method, exc)
            row["PCR"] = float("nan")

        row["overall"] = _overall_score(
            pd.Series(row), bio_keys, batch_keys,
            float(metrics_cfg.get("bio_weight") or 0.6),
            float(metrics_cfg.get("batch_weight") or 0.4),
        )
        logger.info("[%s] overall=%.3f", method, row["overall"])
        rows.append(row)
        del adata

    if not rows:
        raise RuntimeError("scIB benchmark produced no method scores.")

    df = pd.DataFrame(rows).set_index("method")
    out_csv = metrics_csv_path(output_dir, cfg)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv)
    logger.info("Wrote scIB metrics: %s", out_csv)

    integrated = df.drop(index=["unintegrated"], errors="ignore")
    if not integrated.empty and "overall" in integrated.columns:
        winner = str(integrated["overall"].astype(float).idxmax())
        score = float(integrated.loc[winner, "overall"])
        logger.info("Best integration method: %s (overall=%.3f)", winner, score)
        (bench / "best_method.txt").write_text(f"{winner}\t{score:.4f}\n", encoding="utf-8")

    try:
        _plot_results(df, cfg, bench)
    except Exception as exc:
        logger.warning("Could not write benchmark plots: %s", exc)
    return out_csv


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="scIB metrics for Harmony / RPCA / CCA")
    parser.add_argument("--output-dir", required=True, help="Integration output directory")
    parser.add_argument("--config", default=None, help="benchmark_config.yaml")
    parser.add_argument(
        "--methods",
        default="harmony,rpca,cca",
        help="Comma-separated methods (harmony,rpca,cca)",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
    methods = [m.strip().lower() for m in str(args.methods).split(",") if m.strip()]
    run_scib_benchmark(args.output_dir, config_path=args.config, methods=methods)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
