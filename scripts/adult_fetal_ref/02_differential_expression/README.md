# 02 — Differential expression and pathways

## Methods

Having asked *where* organoids sit on the primary maps, we next ask *how* their expression differs from matched primary cells of the same niche.

**Pseudobulk edgeR.** For each organoid Level_2 type paired to adult or fetal, and for each of the three major protocols (Morizane, Takasato, Freedman), we build pseudobulks (sample_id for organoid/adult; specimen for fetal), require at least three replicates and enough cells/counts, and fit an edgeR quasi-likelihood model with age (and QC covariates) centered within group. DEGs are called at FDR < 0.05 and |log₂FC| > 0.5. Genes significant in the same direction in all three protocols are tagged Shared Up or Shared Down (3/3); the rest fall into protocol-specific or inconsistent categories. Summary figures show DEG counts per protocol and the Shared Up / Shared Down / inconsistent breakdown.

**Broad pathway ORA.** Shared Up and Shared Down gene lists (per cell type, separately vs adult and vs fetal) are tested for over-representation against Enrichr libraries (MSigDB Hallmark, GO Biological Process, Reactome, KEGG) with a hypergeometric test and BH FDR. The background is the set of genes actually tested in that cell type’s DEG runs. Results are summarized as up/down pathway dotplots.

**Kidney-focused view.** The same consensus DEGs are then restricted to pathways that either (i) match kidney / transport / metabolism terms by name, or (ii) sit under kidney-related GO ontology parents (development, filtration, tubular transport, podocyte programs, and related branches). Sparse kidney GO terms missing from the genome-wide ORA are re-tested on that smaller gene-set universe and merged in. Per cell type we plot diverging bars: Shared Down to the left, Shared Up to the right (−log₁₀ p).

**Preranked GSEA.** ORA only sees hard DEG cutoffs. To use the full expression gradient, we also run preranked GSEA (gseapy). For each cell type we average logFC across protocols (genes present in ≥2 of 3), rank all tested genes, and enrich Hallmark, GO Biological Process, and Reactome. Positive NES means the pathway tracks organoid-high genes; negative NES tracks primary-high genes. Summary figures are an NES heatmap across cell types and a bar plot of the strongest pathways.

ORA answers “are these DEG lists over-represented in pathway X?”; GSEA asks whether pathway X is coherently shifted along the organoid–primary ranking even when individual genes sit near the cutoff.

## Scripts

```bash
python scripts/02_differential_expression/deg_edger_protocols.py
python scripts/02_differential_expression/pathway_enrichment.py   # ORA
python scripts/02_differential_expression/pathway_kidney.py       # kidney ORA view
python scripts/02_differential_expression/pathway_gsea.py         # preranked GSEA
```

| Script | Figures |
|---|---|
| `deg_edger_protocols.py` | `figures/*_deg_counts_by_protocol.png`, `*_deg_protocol_consensus.png` |
| `pathway_enrichment.py` | `figures/*_pathway_up_down_dotplot.png` |
| `pathway_kidney.py` | `pathways/figures/*_pathway_kidney_diverging_bars_*.png` |
| `pathway_gsea.py` | `gsea/figures/*_gsea_nes_heatmap.png`, `*_gsea_top_pathways.png` |
