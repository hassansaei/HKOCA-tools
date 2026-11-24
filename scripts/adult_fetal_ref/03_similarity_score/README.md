# 03 — Transcriptional similarity

## Methods

Presence and DE say which niches are occupied and which genes differ. Similarity asks how close organoid and primary profiles are when they *do* match.

**Pseudobulk profiles.** For each organoid Level_2 type (optionally split by protocol) and each primary broad / harmonized type, we sum raw counts across cells, convert to log₂-CPM, and work in that space.

**Gene sets.** Two complementary gene panels are used:

- **HVG** — top 2,000 highly variable genes on the reference log₂-CPM matrix (with a light expression floor), so the correlation reflects the main axes of variation in primary tissue  
- **TF** — human transcription-factor genes from the Lambert list (local copy: `TF_names_Lambert_v1.01.csv` / `.txt` in this folder), to emphasize regulatory state rather than effector programs  

**Similarity.** For every organoid type × primary type pair we compute Spearman correlation of the two log₂-CPM vectors on the shared genes in each panel. Full matrices are shown as heatmaps (HVG and TF) against adult and fetal.

**Matched niches.** Many organoid labels map to more than one primary subtype (e.g. collecting-duct-like organoids vs PC / IC-A / IC-B / CNT). For those pairs we build a matched primary profile as the equal-weight mean of the relevant primary columns, then score organoid vs that matched vector. Protocol-wise matched scores (HVG) are plotted as bars; a fetal-vs-adult comparison shows, for shared organoid types, which reference they resemble more closely.

Off-target and undifferentiated pluripotent cells are kept on the heatmaps for context but are not treated as matched renal niches.

## Scripts

```bash
python scripts/03_similarity_score/plot_similarity.py
```

Typical outputs under `03_similarity_score/`:

- `*_similarity_heatmap_HVG.png`, `*_similarity_heatmap_TF.png`  
- `*_similarity_matched_by_protocol.png`  
- fetal-vs-adult matched fidelity panel (when both references are run)
