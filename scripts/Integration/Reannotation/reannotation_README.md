# Post-integration reannotation

Snapseed assigns Level_1 / Level_2 / Level_3 from marker scores on pre-integration clusters. That works well at atlas scale, but it is not 100% accurate: some clusters mix lineages or states that markers alone mislabel. After scPoli, Leiden within each lineage resolves those pockets more clearly, so a manual reannotation step is needed to correct them on the integrated object (without retraining the model).

`run_reannotation.py` applies the agreed remaps and writes a new reannotated `.h5ad`. The source integrated file is left untouched.

## Cluster remaps

| Lineage → cluster | n | To (L1 / L2 / L3) | Why |
|---|---:|---|---|
| Stromal_cells · 13 | 14,413 | Nephron / NPC / NPC | Strong Nephron markers (PAX2/8, EPCAM, LYPD1); low stromal |
| Off_target_cells · 10 | 8,252 | Nephron / Proliferative_nephron | PAX2/8 + high cycle markers; not off-target |
| Nephron · 2_sub0 | 7,439 | Stromal_cells / S2 | Stromal signature (PDGFRA, COL14A1); negative Nephron score |
| Nephron · 14 | 6,143 | Nephron / Proliferative_nephron | Highest proliferation in Nephron + PAX2/8 |
| Nephron · 2_sub1 | 5,581 | Nephron / Proliferative_nephron | PAX2/8 + cycle markers; no stromal/segmental signal |
| Nephron · 2_sub2 | 5,237 | Nephron / NPC | Retained as NPC (PAX2/8, EPCAM, LYPD1) |
| Stromal_cells · 20 | 4,305 | Nephron / NPC | Strong Nephron score; kept NPC (not prolif.) due to CDH6/LYPD1 |
| Nephron · 19 | 3,969 | Nephron / Proliferative_nephron | Clean PAX2/8 + cycle markers; segmental ~0 |
| Nephron · 25 | 2,006 | Off_target_cells | Negative Nephron score; mixed neural/off-target |
| Off_target_cells · 29 | 1,448 | Stromal_cells / S2 | Strong PDGFRA / COL14A1 / PRRX1 |
| Stromal_cells · 26 | 110 | Ureteric_Epithelium / Collecting_duct | AQP2 / GATA3 / EPCAM |
| undifferentiated · 32 | 93 | Endothelial_cells | KDR / PECAM1; no pluripotency markers |
| Off_target_cells · 34 | 37 | Nephron / Proximal_tubule / Dev_PT | CUBN / LRP2 / PAX2/8 |
| Off_target_cells · 35 | 27 | Endothelial_cells | KDR / PECAM1 |

Nephron cluster 2 was split on `X_emb` into three subclusters because stromal and proliferative programs were segregated inside that cluster.

## Layout

```text
reannotation/
├── README.md
└── run_reannotation.py
```

## How to run

```bash
python run_reannotation.py
```
