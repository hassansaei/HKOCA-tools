#!/usr/bin/env Rscript
# edgeR quasi-likelihood F-test: organoid (Query) vs primary (Control)
# Enforces >= 3 replicates per group and includes age (within-group centered/scaled)
# plus pseudobulk QC covariates (log1p(n_cells), median detected genes).
#
# Usage:
#   Rscript edger_qlf_org_vs_primary.R <counts.tsv> <meta.tsv> <outdir>
#
# Positive logFC = higher in organoid (Query) vs primary (Control).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: edger_qlf_org_vs_primary.R <counts.tsv> <meta.tsv> <outdir>")
}
counts_path <- args[[1]]
meta_path <- args[[2]]
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
})

counts <- as.matrix(read.delim(counts_path, row.names = 1, check.names = FALSE))
meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- as.character(meta$sample)
meta <- meta[colnames(counts), , drop = FALSE]
meta$group <- factor(meta$group, levels = c("Control", "Query"))

n_ctrl <- sum(meta$group == "Control")
n_query <- sum(meta$group == "Query")
if (n_ctrl < 3) stop(sprintf("Need >=3 Control replicates (got %d)", n_ctrl))
if (n_query < 3) stop(sprintf("Need >=3 Query replicates (got %d)", n_query))

# Helper: impute non-finite with within-group median (fallback 0) and scale within group
# Centering age within group adjusts for developmental/donor age variation without
# introducing collinearity between age and the group (organoid vs primary) contrast.
within_group_z <- function(x, grp) {
  x <- as.numeric(x)
  out <- rep(0, length(x))
  for (lv in levels(grp)) {
    idx <- which(grp == lv)
    vals <- x[idx]
    fin <- is.finite(vals)
    if (!any(fin)) {
      out[idx] <- 0
      next
    }
    vals[!fin] <- median(vals[fin])
    s <- sd(vals)
    if (!is.finite(s) || s < 1e-8) {
      out[idx] <- 0
    } else {
      out[idx] <- (vals - mean(vals)) / s
    }
  }
  out
}

scale_safe <- function(x) {
  x <- as.numeric(x)
  fin <- is.finite(x)
  if (!any(fin)) return(rep(0, length(x)))
  x[!fin] <- median(x[fin])
  s <- sd(x)
  if (!is.finite(s) || s < 1e-8) return(rep(0, length(x)))
  as.numeric(scale(x))
}

for (col in c("age_val", "n_cells", "med_genes")) {
  if (!(col %in% colnames(meta))) meta[[col]] <- NA_real_
}

meta$age_z <- within_group_z(meta$age_val, meta$group)
meta$ncells_z <- scale_safe(log1p(as.numeric(meta$n_cells)))
meta$medg_z <- scale_safe(as.numeric(meta$med_genes))

covars <- c("age_z", "ncells_z", "medg_z")
covars <- covars[sapply(covars, function(c) sd(meta[[c]]) > 1e-6)]

# Check rank deficiency before fitting; drop QC covars if needed (keep age_z first)
build_design <- function(cv) {
  form <- as.formula(
    paste("~ group", if (length(cv) > 0) paste("+", paste(cv, collapse = " + ")) else "")
  )
  mm <- model.matrix(form, data = meta)
  list(form = form, mm = mm)
}

res <- build_design(covars)
while (qr(res$mm)$rank < ncol(res$mm) && length(covars) > 0) {
  covars <- covars[-length(covars)]
  res <- build_design(covars)
}
design <- res$mm
colnames(design) <- make.names(colnames(design))

y <- DGEList(counts = counts, group = meta$group)
keep <- filterByExpr(y, group = meta$group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

coef_query <- grep("groupQuery", colnames(design), value = TRUE)
if (length(coef_query) != 1) {
  stop(sprintf("Expected one groupQuery coefficient, got: %s", paste(colnames(design), collapse = ",")))
}

qlf <- glmQLFTest(fit, coef = coef_query)
tab <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))
tab$gene <- rownames(tab)
tab$logFC <- as.numeric(tab$logFC)
tab$logCPM <- as.numeric(tab$logCPM)
tab$F <- as.numeric(tab$F)
tab$PValue <- as.numeric(tab$PValue)
tab$FDR <- as.numeric(tab$FDR)
tab$significant <- (!is.na(tab$FDR)) & (tab$FDR < 0.05) & (abs(tab$logFC) > 0.5)
tab$direction <- ifelse(!tab$significant, "NS", ifelse(tab$logFC > 0, "Up", "Down"))

write.table(tab, file.path(outdir, "deg_full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sig <- tab[tab$significant %in% TRUE, , drop = FALSE]
write.table(sig, file.path(outdir, "deg_sig.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sidecar <- list(
  n_control = as.integer(n_ctrl),
  n_query = as.integer(n_query),
  n_genes_tested = nrow(y),
  n_sig = nrow(sig),
  n_sig_up = sum(sig$logFC > 0, na.rm = TRUE),
  n_sig_down = sum(sig$logFC < 0, na.rm = TRUE),
  formula = paste(deparse(res$form), collapse = ""),
  contrast = "Query(organoid) - Control(primary)",
  covariates = covars
)
write(toJSON(sidecar, auto_unbox = TRUE, pretty = TRUE), file.path(outdir, "design.json"))
cat(sprintf(
  "edgeR-QLF: %d genes tested, %d ctrl + %d query, %d sig (%d up, %d down)\n",
  nrow(y), n_ctrl, n_query, nrow(sig), sidecar$n_sig_up, sidecar$n_sig_down
))
