# 00 — Reference atlases and UMAPs

## Methods

We start from two primary kidney references that later stages use as the ground truth for presence, DE, and similarity.

**Adult.** The adult map comes from the Human Kidney Cell Atlas (HCKA). We keep only cells labelled `disease == normal` with `sampled_site_condition` in `{healthy, adjacent}` and drop diseased tissue. That leaves roughly 599k cells. Broad annotations (`broad_annotation`) are taken as provided by the atlas; we do not re-annotate. The UMAP shown here is the frozen HCKA scANVI embedding, subset to that healthy/adjacent cohort.

**Fetal.** The fetal map merges three public datasets (Stewart, Hochane, Tran) restricted to post-fertilization weeks **9–17**. Author labels are mapped to a shared `harmonized_celltype` vocabulary. Tran was processed from GSE124472 10x matrices, QC’d, and annotated with DevKidCC; Stewart and Hochane enter as prebuilt objects. Integration is Harmony on 50 PCs from 3,000 HVGs (batch = dataset), followed by neighbors, UMAP, and Leiden clustering. The object used downstream has 47,999 cells.

**Organoid soft landing.** To place organoid cells in the same UMAP space as each reference, we compute a weighted k-nearest-neighbor graph in the shared latent space (k = 100, Jaccard-square weights, query→reference only), transfer reference labels, and land each organoid cell by a soft average of its neighbors’ UMAP coordinates. Figures use a slightly sharpened landing (top-10 neighbors, weight power 2) so structure is readable; the full-k soft coordinates are kept for analysis.

Presence scores (blue→red UMAPs) are **not** computed here — that lives in section 01.

## Scripts

```bash
# Upstream processing (adult filter audit; fetal figure sync)
python scripts/00_reference_umaps/upstream/process_adult.py
python scripts/00_reference_umaps/upstream/process_fetal.py

# Annotation UMAPs
python scripts/00_reference_umaps/plot_broad_umaps.py

# Organoid-on-reference overlays (Level_2, confidence, transferred labels)
python scripts/00_reference_umaps/plot_mapping_umaps.py

# Rebuild wkNN / soft landing only if needed (slow)
# python scripts/00_reference_umaps/compute_mapping.py
```

Outputs: `00_reference_umaps/figures/`, `00_reference_umaps/full_projection/{adult,fetal}/`.
