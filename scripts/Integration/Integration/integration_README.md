# scPoli integration

This folder holds the production integration of the HKOCA organoid atlas: one annotated object in, one scPoli-integrated atlas out (plus the trained reference model).

It sits next to the multi-method scIB benchmark. That benchmark compared classical and deep methods on Level_2 labels; scPoli came out on top on the overall scorecard, which is why it was chosen for the final atlas.

## Why scPoli

The atlas mixes many studies, protocols, and samples. Integration has to do two things at once: mix batches so technical structure fades, and keep biology so cell types stay separable.

scPoli is a deep generative model (scArches) built for atlas-scale, population-level integration. It learns a shared latent space conditioned on study and sample, while using Level_2 cell-type labels to shape prototypes. That combination is a good fit for a multi-cohort organoid atlas: conditions absorb batch, prototypes keep cell-type structure.

In our Level_2 scIB benchmark, scPoli edged the field on the weighted overall score. Relative to strong peers such as scANVI, it mixed samples more convincingly (higher batch metrics such as iLISI, kBET, and graph connectivity) while still preserving Level_2 biology. That balance is what we need for a reference atlas that later supports projection, reannotation, and cross-dataset comparison.

## What happens in the run

1. Load the annotated atlas counts (raw counts in `layers['counts']`; log-normalized `X` is not used for training).
2. Clean required metadata (`study`, `sample_id`, Level_1 / Level_2) and subset to 4000 highly variable genes (`seurat_v3`, batched by study).
3. Build a PCA baseline (neighbors + UMAP) only for before/after comparison plots.
4. Train scPoli in two phases on the count matrix:
   - pretraining (80 epochs) without prototypes
   - prototype fine-tuning (20 epochs) with Level_2 labels
   Conditions are `study` and `sample_id`. Reconstruction uses a negative-binomial loss.
5. Extract the latent embedding (`X_emb`), compute neighbors and UMAP on that space, and write the integrated object plus model checkpoints.

The integrated object carries the scPoli latent (`X_emb`), PCA and scPoli UMAPs, and the Snapseed Level_1 / Level_2 / Level_3 labels that entered training. The saved model is the frozen reference used for later query surgery and projection.

## Layout

```text
scpoli/
├── README.md
├── config.yaml
└── run_scpoli_integration.py
```

`config.yaml` sets input/output locations, HVG and training hyperparameters, and plotting style. `run_scpoli_integration.py` is the single entry point (color palettes live inside the script).

## How to run

Use the `atlas` conda environment on a GPU node:

```bash
python run_scpoli_integration.py
# or
python run_scpoli_integration.py --config config.yaml
```

## Outputs

Configured under `paths` in `config.yaml`:

- integrated `.h5ad` with `X_emb` and UMAPs
- scPoli reference model (and a post-pretrain checkpoint)
- comparison PNGs (PCA vs scPoli, colored by sample, study, Level_1/2/3, and covariates)
