# 01 — Presence and coverage

## Methods

This section asks a simple question: which primary kidney niches do organoids reach, and how strongly?

**Composition heatmap.** Organoid cells are mapped onto each reference with a prior wkNN label transfer. We keep confident assignments (confidence ≥ 0.5) and summarize, per primary cell type and per differentiation protocol (Morizane, Takasato, Freedman, plus all organoids pooled), what fraction of organoid cells land there. Immune lineages are collapsed into one “immune cells” row so rare subtypes do not dominate the display. Adjacent bars show how abundant each type is in the primary atlas itself.

**Presence score.** Separately, we score every *reference* cell for how strongly organoids “visit” its neighborhood. Using the saved query→reference wkNN graph, we sum organoid edge weight onto each reference cell and smooth with a random walk with restart on the reference graph (α = 0.1, 100 rounds), then log1p and scale to [0, 1]. That per-cell score is what colors the blue→red UMAP.

**Size-normalized bars.** Raw mean presence favors tiny niches that happen to sit under a few high-scoring organoids. For each primary cell type we therefore report  

`presence_norm = mean_presence × n / (n + 100)`  

so rare types need stronger local evidence to look “present.” Types are labelled Present / Low / Absent with fixed cutoffs (0.25 / 0.05).

Adult and fetal are run the same way; adult immune broad types are collapsed for the bar plot as for the heatmap.

## Scripts

```bash
python scripts/01_presence_coverage/plot_presence.py
```

One script writes all three figure families for adult and fetal:

| Figure | Role |
|---|---|
| `*_presence_by_protocol_heatmap.png` | confident composition × protocol |
| `*_presence_by_celltype.png` | size-normalized presence bars |
| `*_presence_umap.png` | blue→red presence on the reference UMAP |

Tables: `01_presence_coverage/tables/`. The wkNN graph is reused from `00_reference_umaps/full_projection/{adult,fetal}/tables/wknn.npz`.
