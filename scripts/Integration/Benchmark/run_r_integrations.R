# R integrations for one subsample iteration (Seurat v3/v4 + CSS).
# Reads method toggles from config.yaml (or BENCH_METHODS env as fallback).

library(Seurat)
library(Matrix)

ITER_ID <- as.integer(Sys.getenv("BENCHMARK_ITER_ID", unset = "1"))
BENCH_DIR <- Sys.getenv(
  "BENCH_DIR",
  unset = "/data-master/pure-workspace/labss/hmami/atlas_shortcut_new_start/Organization/integration/benchmark"
)
CONFIG_YAML <- Sys.getenv("CONFIG_YAML", unset = file.path(BENCH_DIR, "config.yaml"))

# Prefer YAML; fall back to BENCH_METHODS env from the orchestrator
if (requireNamespace("yaml", quietly = TRUE) && file.exists(CONFIG_YAML)) {
  cfg <- yaml::read_yaml(CONFIG_YAML)
  flags <- cfg$methods
  enabled <- function(name) isTRUE(flags[[name]])
  cat("Loaded method flags from", CONFIG_YAML, "\n")
} else {
  BENCH_METHODS <- strsplit(Sys.getenv("BENCH_METHODS", unset = ""), ",", fixed = TRUE)[[1]]
  BENCH_METHODS <- trimws(BENCH_METHODS)
  BENCH_METHODS <- BENCH_METHODS[nzchar(BENCH_METHODS)]
  enabled <- function(name) length(BENCH_METHODS) == 0 || name %in% BENCH_METHODS
  cat("Using BENCH_METHODS env (yaml package or config missing)\n")
}

CHECKPOINT_DIR <- file.path(BENCH_DIR, "checkpoints", paste0("subsample_iter_", ITER_ID))
setwd(CHECKPOINT_DIR)
cat("Working directory:", CHECKPOINT_DIR, "\n")

want_seurat <- enabled("seurat_v3") || enabled("seurat_v4")
want_css <- enabled("css_pearson") || enabled("css_spearman")

if (!want_seurat && !want_css) {
  cat("No R methods enabled — nothing to do.\n")
  quit(save = "no", status = 0)
}

load_seurat <- function(normalize = FALSE) {
  cat("  Reading MTX matrix and metadata...\n")
  counts_matrix <- readMM("heoca_counts.mtx")
  genes <- read.csv("heoca_genes.csv", header = FALSE)$V1
  metadata <- read.csv("heoca_metadata.csv", row.names = 1)

  rownames(counts_matrix) <- genes
  colnames(counts_matrix) <- rownames(metadata)

  seu <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)
  VariableFeatures(seu) <- rownames(seu)

  if (normalize) {
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
  }

  rm(counts_matrix, genes, metadata)
  gc()
  return(seu)
}

if (want_seurat) {
  cat("\nLoading master Seurat object\n")
  seu_master <- load_seurat(normalize = FALSE)
  cat("Master object ready. Cells:", ncol(seu_master), "| Genes:", nrow(seu_master), "\n")
  obj.list <- SplitObject(seu_master, split.by = "sample_id")

  if (enabled("seurat_v3")) {
    cat("\nSeurat v3 (CCA)\n")
    t0 <- Sys.time()
    obj.list.v3 <- lapply(X = obj.list, FUN = function(x) {
      x <- NormalizeData(x, verbose = FALSE)
      VariableFeatures(x) <- rownames(x)
      return(x)
    })
    cat("  Finding CCA Anchors...\n")
    anchors_v3 <- FindIntegrationAnchors(object.list = obj.list.v3, dims = 1:30, reduction = "cca")
    seu_v3 <- IntegrateData(anchorset = anchors_v3, dims = 1:30)
    seu_v3 <- ScaleData(seu_v3, verbose = FALSE)
    seu_v3 <- RunPCA(seu_v3, npcs = 50, verbose = FALSE)
    write.csv(Embeddings(seu_v3, reduction = "pca"), "seurat_v3_embeddings.csv", quote = FALSE)
    cat("Seurat v3 complete | Time:", round(Sys.time() - t0, 2), "\n")
    rm(anchors_v3, seu_v3, obj.list.v3)
    gc()
  } else {
    cat("SKIP seurat_v3 (disabled)\n")
  }

  if (enabled("seurat_v4")) {
    cat("\nSeurat v4 (RPCA)\n")
    t0 <- Sys.time()
    obj.list.v4 <- lapply(X = obj.list, FUN = function(x) {
      x <- NormalizeData(x, verbose = FALSE)
      VariableFeatures(x) <- rownames(x)
      x <- ScaleData(x, verbose = FALSE)
      x <- RunPCA(x, features = rownames(x), verbose = FALSE)
      return(x)
    })
    cat("  Finding RPCA Anchors...\n")
    anchors_v4 <- FindIntegrationAnchors(object.list = obj.list.v4, dims = 1:30, reduction = "rpca")
    seu_v4 <- IntegrateData(anchorset = anchors_v4, dims = 1:30)
    seu_v4 <- ScaleData(seu_v4, verbose = FALSE)
    seu_v4 <- RunPCA(seu_v4, npcs = 50, verbose = FALSE)
    write.csv(Embeddings(seu_v4, reduction = "pca"), "seurat_v4_embeddings.csv", quote = FALSE)
    cat("Seurat v4 complete | Time:", round(Sys.time() - t0, 2), "\n")
    rm(anchors_v4, seu_v4, obj.list.v4)
    gc()
  } else {
    cat("SKIP seurat_v4 (disabled)\n")
  }

  rm(obj.list, seu_master)
  gc()
}

if (want_css) {
  cat("\nReloading master object with PCA for CSS\n")
  seu_master <- load_seurat(normalize = TRUE)

  if (!requireNamespace("simspec", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("quadbio/simspec")
  }
  library(simspec)

  if (enabled("css_pearson")) {
    cat("\nCSS (Pearson)\n")
    t0 <- Sys.time()
    seu_css_pearson <- cluster_sim_spectrum(seu_master, label_tag = "sample_id", corr_method = "pearson")
    write.csv(Embeddings(seu_css_pearson, reduction = "css"), "css_pearson_embeddings.csv", quote = FALSE)
    cat("CSS Pearson complete | Time:", round(Sys.time() - t0, 2), "\n")
    rm(seu_css_pearson)
    gc()
  } else {
    cat("SKIP css_pearson (disabled)\n")
  }

  if (enabled("css_spearman")) {
    cat("\nCSS (Spearman)\n")
    t0 <- Sys.time()
    seu_css_spearman <- cluster_sim_spectrum(seu_master, label_tag = "sample_id", corr_method = "spearman")
    write.csv(Embeddings(seu_css_spearman, reduction = "css"), "css_spearman_embeddings.csv", quote = FALSE)
    cat("CSS Spearman complete | Time:", round(Sys.time() - t0, 2), "\n")
    rm(seu_css_spearman)
    gc()
  } else {
    cat("SKIP css_spearman (disabled)\n")
  }

  rm(seu_master)
  gc()
}

cat("\nR integrations complete. CSVs saved to:", CHECKPOINT_DIR, "\n")
