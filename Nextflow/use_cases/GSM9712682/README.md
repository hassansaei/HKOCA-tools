# GSM9712682

Example query dataset used to check the Nextflow HKOCA pipeline end to end.

## What this sample is

Single-cell RNA-seq from **human iPSC-derived ureteric bud (UB) organoids**, wild-type, collected at **day 12** of differentiation. The deposit is public on GEO as sample **GSM9712682** within series **GSE329881** (contact: Veronika Sander, University of Auckland).

Differentiation protocol is tagged `vanichapol` in the atlas metadata. Libraries were prepared with **Parse Biosciences Evercode WT Mini v3** and sequenced on a **MGI DNBSEQ-G400** (platform GPL28038). Gene annotation uses **GRCh38**. Counts in the GEO deposit are already TrailMaker-processed, so CellBender was not run (`run_cellbender=False`).

Condition label in the sample sheet: **WT**.

## Metadata (from `sample_info`)

| Field | Value |
|---|---|
| `sample_id` / `geo_gsm` | GSM9712682 |
| `study` / `geo_gse` | GSE329881 |
| `source` | hiPSC |
| `cell_line` | hiPSC_derived |
| `cell_type` | UB_organoids |
| `diff_protocol` | vanichapol |
| `Age` | d12 |
| `genotype` / `condition` | wildtype / WT |
| `sc_protocol` | Parse_Evercode_WT_Mini_v3 |
| `library_kit` | Evercode_WT_Mini_v3 |
| `sequencing` | DNBSEQ-G400 |
| `platform_id` | GPL28038 |
| `genome_build` | GRCh38 |
| `type` | scRNAseq |
| `run_cellbender` | False |

Template for a local re-run: [`sample_info.template.csv`](sample_info.template.csv) (set `data_path` to your h5ad).

## What happened in this run

| Stage | Result |
|---|---|
| Harmonize → QC | 4,568 cells in; **4,246** kept (~7% removed) |
| Annotation | Snapseed Level 1–3 labels |
| Integration | SCT prep (single sample; Harmony/RPCA/CCA skipped) |
| Projection | 4,246 cells onto the HKOCA atlas; **24,100 / 24,100** genes shared with the reference model |

## Folder contents

Lightweight outputs only (plots, small tables, Nextflow reports). No expression objects or model weights.

| Path | Contents |
|---|---|
| `01_harmonize/` | Harmonize metadata + versions |
| `02_qc/` | QC / doublet summaries and violin plots |
| `03_annotation/` | Level UMAPs (res 0.6) + Snapseed assignments |
| `04_integration/` | Elbow / silhouette / Level UMAPs and proportions |
| `05_projection/` | Atlas overlay plots + projection summary JSON/CSV |
| `pipeline_info/` | Nextflow trace, report, timeline, DAG |
