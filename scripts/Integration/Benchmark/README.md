# HKOCA integration benchmark (Level_2)

Benchmark of single-cell integration methods on the HKOCA organoid atlas, scored with **scIB** against **Level_2** cell-type labels.

---

## Framework and citation

We follow the **scIB** (single-cell Integration Benchmarking) framework of Luecken et al.:

> Luecken, M.D., Büttner, M., Chaichoompu, K. *et al.* Benchmarking atlas-level data integration in single-cell genomics. *Nat Methods* **19**, 41–50 (2022).  
> https://doi.org/10.1038/s41592-021-01336-8  
> Software: https://github.com/theislab/scib

Overall score (as used here):

```text
overall = 0.6 × mean(biological conservation metrics)
        + 0.4 × mean(batch-correction metrics)
```

---

## What this benchmark does

1. Load the clean annotated atlas (`paths.master_h5ad` in `config.yaml`).
2. Split `sample_id` into **N non-overlapping windows** (default 10; shuffle seed **42**).
3. Per iteration: recompute **3000 HVGs** (`seurat_v3` flavor, batched by `sample_id`).
4. Integrate with every method set to `true` in `config.yaml`.
5. Score those methods with **scIB** (`label_key=Level_2`, `batch_key=sample_id`).
6. Aggregate mean ± SD across iterations and plot the master scorecard.

---

## Config (`config.yaml`)

All runtime settings live in `config.yaml` (read directly by the orchestrator and each script):

- **paths** — atlas h5ad, Python binary, `Rscript`
- **benchmark** — `n_iters`, seeds, HVG count, batch/label keys
- **methods** — per-method `true` / `false` toggles

Example: CPU-only classical methods

```yaml
methods:
  pca: true
  harmony: true
  scanorama: true
  combat: true
  bbknn: true
  scvi: false
  scanvi: false
  scpoli: false
  seurat_v3: true
  seurat_v4: true
  css_pearson: true
  css_spearman: true
```

Notes:

- Disabled methods are skipped (not trained, not scored).
- If `scanvi: true` but `scvi: false`, scVI is still trained once as a warm-start dependency (not scored).
- With all of `scvi` / `scanvi` / `scpoli` set to `false`, no GPU is required.

---

## Methods

| Method | Type | Notes |
|--------|------|--------|
| **pca** | Unintegrated baseline | Scaled PCA + neighbors |
| **harmony** | Batch correction | Harmony on PCA |
| **scanorama** | Batch correction | Scanorama integrate |
| **combat** | Batch correction | ComBat → PCA |
| **bbknn** | Graph integration | BBKNN on PCA |
| **scvi** | Deep generative | Label-free scVI (GPU preferred) |
| **scanvi** | Deep generative | Warm-start from scVI; labels = `Level_2` |
| **scpoli** | Deep generative | scArches scPoli; labels = `Level_2` |
| **seurat_v3** | Anchor (R) | CCA integration |
| **seurat_v4** | Anchor (R) | RPCA integration |
| **css_pearson** | Spectrum (R) | Cluster similarity spectrum (Pearson) |
| **css_spearman** | Spectrum (R) | Cluster similarity spectrum (Spearman) |

---

## scIB metrics

### Biological conservation (weight 0.6)

| Metric | Meaning |
|--------|---------|
| **ARI** / **NMI** | Agreement of Leiden clusters with Level_2 labels |
| **ASW_celltype** | Silhouette of Level_2 in the embedding |
| **iso_label_F1** / **iso_label_ASW** | Recovery of batch-isolated labels |
| **cLISI** | Cell-type local inverse Simpson index |

### Batch correction (weight 0.4)

| Metric | Meaning |
|--------|---------|
| **ASW_batch** | Batch mixing by silhouette (higher = better mixing) |
| **iLISI** | Batch local inverse Simpson index |
| **graph_connectivity** | Connectivity of same-label cells across batches |
| **kBET** | k-nearest-neighbour batch effect test |
| **PCR** | Principal-component regression (batch variance removed vs pre-integration) |

---

## Directory layout

```text
benchmark/
├── README.md
├── config.yaml                 # only config file (paths, seeds, method toggles)
├── scripts/
│   ├── run_benchmark.sh        # local orchestrator (entry point)
│   ├── run_integrations.py     # Python methods + MTX export
│   ├── run_r_integrations.R    # Seurat / CSS
│   ├── run_metrics.py          # scIB scoring
│   └── visualisation_results.py
├── results/                    # mean/SD tables + scorecard PNGs
├── checkpoints/                # per-iteration outputs
├── logs/
└── tmp/
```

All scripts and the orchestrator read `config.yaml` directly (no separate Python config module).

---

## How to run

Requires the `atlas` conda env (Python + R Seurat/simspec). Edit `config.yaml`, then:

```bash
bash Organization/integration/benchmark/scripts/run_benchmark.sh
```

The orchestrator loops iterations locally:

1. Python integrations (enabled methods only)
2. R integrations (if any R method is enabled)
3. scIB metrics
4. After all iterations: aggregate CSVs + figures → `results/`

Per-iteration logs: `logs/iter_N.log`.

### Outputs per iteration

```text
checkpoints/subsample_iter_N/reference.h5ad
checkpoints/subsample_iter_N/checkpoint_<method>.h5ad
checkpoints/subsample_iter_N/scib_raw_metrics_full.csv
```

### Aggregated results

```text
results/scib_mean_metrics.csv
results/scib_std_metrics.csv
results/MASTER_SCORECARD_size.png
results/scib_overall_errorbars_size.png
results/scib_mean_sd_heatmap_size.png
```

Archived tables/figures in `results/` come from the completed Level_2 run (10 iterations, all methods).

---

## Reproducibility note

Subsampling, HVGs, and classical methods use fixed seeds (`subsample_seed` / `random_state` in `config.yaml`). GPU deep-learning methods can still vary slightly across hardware (CUDA non-determinism).
