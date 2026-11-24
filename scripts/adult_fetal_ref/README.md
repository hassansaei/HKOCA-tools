# Organoid vs primary atlas — analysis overview

Scripts that compare the kidney organoid atlas to adult and fetal primary kidney references. The pipeline is ordered: first build the maps and land organoids on them, then ask which niches are occupied, which genes and pathways diverge, and how close matched niches are transcriptionally.

Figures and tables land under `organoid_vs_primary_atlas/{00,01,02,03}_*/`. Each numbered folder below has its own README with methods and the exact commands used to regenerate outputs.

---

## Sections

### 00 — Reference atlases and UMAPs

Adult reference is the Human Kidney Cell Atlas (HCKA) restricted to healthy/adjacent normal tissue (~599k cells), using the atlas’s broad labels and frozen scANVI UMAP. Fetal reference merges Stewart, Hochane, and Tran for post-fertilization weeks 9–17 (~48k cells), with shared `harmonized_celltype` labels and Harmony integration. Organoid cells are soft-landed onto each reference via a weighted k-nearest-neighbor graph in latent space (label transfer plus soft UMAP coordinates). Presence scoring is not done here — that belongs in 01.

Details: [`00_reference_umaps/README.md`](00_reference_umaps/README.md)

### 01 — Presence / coverage

Organoids are placed on the adult and fetal references, and each primary niche is scored for how strongly it is occupied: confident wkNN composition by protocol, then a size-normalized presence score, shown as blue→red UMAPs. Cell types are called present, low, or absent from those scores.

Details: [`01_presence_coverage/README.md`](01_presence_coverage/README.md)

### 02 — Differential expression and pathways

Pseudobulk edgeR compares organoids to matched adult/fetal tissue per cell type × protocol (Morizane, Takasato, Freedman). Genes that differ—and those Shared Up/Down in all three protocols—are tested with pathway ORA, a kidney-focused panel, and preranked GSEA. ORA asks whether DEG lists are enriched in a pathway; GSEA asks whether that pathway shifts coherently along the full organoid–primary ranking.

Details: [`02_differential_expression/README.md`](02_differential_expression/README.md)

### 03 — Transcriptional similarity

Organoid and primary average profiles are compared by Spearman correlation on log₂-CPM: all HVGs for global fidelity and transcription factors (Lambert list, local copy in this folder) for core identity, with matched-niche and protocol views against adult and fetal.

Details: [`03_similarity_score/README.md`](03_similarity_score/README.md)

---

## Shared helpers

| Path | Role |
|---|---|
| `lib_config.py` | Paths, cell-type maps, colors, TF loader, figure defaults |
| `mapping_data.py` | Load reference / organoid objects and saved wkNN artifacts |
| `wknn_mapping/` | Weighted kNN soft landing and label transfer |
| `r/edger_qlf_org_vs_primary.R` | edgeR QLF model called from the DE Python wrapper |

Run scripts from the project root (`organoid_vs_primary_atlas/`’s parent, or as documented in each section) so relative imports resolve. Prefer the same conda env used for the atlas analyses (`atlas`).
