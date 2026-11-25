#!/usr/bin/env Rscript
# Classify Tran GSE124472 fetal matrices with DevKidCC v0.4.0 and export
# predictions + counts for Python h5ad assembly.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DevKidCC)
})

# Local paths under organoid_vs_primary_atlas (no other project folders)
ROOT <- "/data-master/pure-workspace/labss/hmami/atlas_shortcut_new_start/organoid_vs_primary_atlas/00_reference_umaps/inputs/fetal"
RAW <- file.path(ROOT, "tran_raw")
OUT <- file.path(ROOT, "devkidcc_tran")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("DevKidCC version:", as.character(packageVersion("DevKidCC")), "\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")

samples <- list(
  list(id = "GSM3534656_H17w_Z1", week = "H17w", zone = "Z1"),
  list(id = "GSM3534657_H17w_Z2", week = "H17w", zone = "Z2"),
  list(id = "GSM3534658_H15w_Z1", week = "H15w", zone = "Z1"),
  list(id = "GSM3534659_H15w_Z2", week = "H15w", zone = "Z2")
)

# Light QC consistent with prior Tran processing
MIN_GENES <- 200
MAX_GENES <- 6000
MIN_COUNTS <- 500
MAX_MITO <- 20

read_one <- function(sample) {
  path <- file.path(RAW, sample$id)
  cat("Reading", sample$id, "...\n")
  counts <- Read10X(path)
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample$id,
    min.cells = 0,
    min.features = 0
  )
  # prefix barcodes with sample id
  obj <- RenameCells(obj, add.cell.id = sample$id)
  obj$sample_id <- sample$id
  obj$week <- sample$week
  obj$zone <- sample$zone
  obj$dataset <- "Tran"
  obj$author <- "Tran et al."
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  keep <- which(
    obj$nFeature_RNA >= MIN_GENES &
      obj$nFeature_RNA <= MAX_GENES &
      obj$nCount_RNA >= MIN_COUNTS &
      obj$percent.mt <= MAX_MITO
  )
  cat("  ", ncol(obj), "->", length(keep), "cells after QC\n")
  obj[, keep]
}

objs <- lapply(samples, read_one)
# Merge all (DevKidCC is per-cell; no integration required)
seu <- merge(objs[[1]], y = objs[-1], project = "Tran_GSE124472")
# Seurat v5 can keep split layers; join counts for a single matrix
if (inherits(seu[["RNA"]], "Assay5")) {
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
}
cat("Merged object:", ncol(seu), "cells x", nrow(seu), "genes\n")

# Required for DevKidCC / Seurat v5
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Patch GeneSummary for Seurat v5 (object[[col]] returns data.frame)
source(file.path(ROOT, "fix_GeneSummary.R"), local = TRUE)
assignInNamespace("GeneSummary", GeneSummary, ns = "DevKidCC")
cat("Injected Seurat-v5-safe GeneSummary into DevKidCC namespace\n")

cat("Running DKCC()...\n")
seu <- DKCC(seu)

md <- seu[[]]
stopifnot("DKCC" %in% colnames(md), "LineageID" %in% colnames(md))

cat("\nLineageID:\n")
print(sort(table(md$LineageID), decreasing = TRUE))
cat("\nDKCC:\n")
print(sort(table(md$DKCC), decreasing = TRUE))

# Export metadata
write.csv(md, file.path(OUT, "devkidcc_metadata.csv"), row.names = TRUE)

# Export raw counts (genes x cells) for h5ad
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
# Market Matrix triplet for Python
Matrix::writeMM(counts, file.path(OUT, "counts.mtx"))
write.table(
  rownames(counts),
  file.path(OUT, "genes.tsv"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
write.table(
  colnames(counts),
  file.path(OUT, "barcodes.tsv"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Also save Seurat RDS
saveRDS(seu, file.path(OUT, "tran_devkidcc.rds"))

# Compact prediction table
pred <- data.frame(
  barcode = rownames(md),
  sample_id = md$sample_id,
  week = md$week,
  zone = md$zone,
  LineageID = md$LineageID,
  DKCC = md$DKCC,
  NephronID = if ("NephronID" %in% colnames(md)) md$NephronID else NA,
  SegmentID = if ("SegmentID" %in% colnames(md)) md$SegmentID else NA,
  StromaID = if ("StromaID" %in% colnames(md)) md$StromaID else NA,
  stringsAsFactors = FALSE
)
write.csv(pred, file.path(OUT, "devkidcc_predictions.csv"), row.names = FALSE)

sink(file.path(OUT, "devkidcc_summary.txt"))
cat("DevKidCC version:", as.character(packageVersion("DevKidCC")), "\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("n_cells:", ncol(seu), "\n")
cat("n_genes:", nrow(seu), "\n\n")
cat("LineageID\n")
print(sort(table(md$LineageID), decreasing = TRUE))
cat("\nDKCC\n")
print(sort(table(md$DKCC), decreasing = TRUE))
sink()

cat("\nWrote outputs under", OUT, "\n")
