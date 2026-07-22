# HKOCA Integration Benchmarking - R Integration Methods
# Configuration: benchmark_config.yaml

library(Seurat)
library(Matrix)

if (!requireNamespace("yaml", quietly = TRUE)) {
  install.packages("yaml")
}
library(yaml)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 2 && args[1] == "--config") {
  args[2]
} else {
  Sys.getenv("BENCHMARK_CONFIG", unset = file.path(script_dir, "benchmark_config.yaml"))
}
cfg <- yaml::read_yaml(config_path)

ITER_ID <- 1
ITER_ID <- as.integer(Sys.getenv(cfg$iteration$slurm_env_var, unset = as.character(cfg$iteration$default_iter_id)))

checkpoint_pattern <- cfg$iteration$checkpoint_subdir_pattern
CHECKPOINT_DIR <- file.path(
  cfg$paths$checkpoints_dir,
  gsub("\\{iter_id\\}", as.character(ITER_ID), checkpoint_pattern)
)
setwd(CHECKPOINT_DIR)
cat("Working directory:", CHECKPOINT_DIR, "\n")

BATCH_KEY <- cfg$metadata$batch_key
R_CFG <- cfg$r_integration
SEURAT_DIMS <- 1:R_CFG$seurat_dims
SEURAT_NPCS <- R_CFG$seurat_npcs
EXPORT_CFG <- cfg$export
EMBEDDING_FILES <- R_CFG$embedding_files


load_seurat <- function(normalize = FALSE) {
  cat("  Reading MTX matrix and metadata...\n")
  counts_matrix <- readMM(EXPORT_CFG$mtx_counts)
  genes          <- read.csv(EXPORT_CFG$mtx_genes, header = FALSE)$V1
  metadata       <- read.csv(EXPORT_CFG$mtx_metadata, row.names = 1)

  rownames(counts_matrix) <- genes
  colnames(counts_matrix) <- rownames(metadata)

  seu <- CreateSeuratObject(counts = counts_matrix, meta.data = metadata)
  VariableFeatures(seu) <- rownames(seu)

  if (normalize) {
    seu <- NormalizeData(seu, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    seu <- RunPCA(seu, npcs = SEURAT_NPCS, verbose = FALSE)
  }

  rm(counts_matrix, genes, metadata)
  gc()
  return(seu)
}


cat("\nLoading master Seurat object\n")
seu_master <- load_seurat(normalize = FALSE)
cat("Master object ready. Cells:", ncol(seu_master), "| Genes:", nrow(seu_master), "\n")

obj.list <- SplitObject(seu_master, split.by = BATCH_KEY)


cat("\nCCA integration\n")
t0 <- Sys.time()

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  VariableFeatures(x) <- rownames(x)
  return(x)
})

cat("  Finding CCA anchors...\n")
anchors_cca <- FindIntegrationAnchors(object.list = obj.list, dims = SEURAT_DIMS, reduction = "cca")
seu_cca     <- IntegrateData(anchorset = anchors_cca, dims = SEURAT_DIMS)
seu_cca     <- ScaleData(seu_cca, verbose = FALSE)
seu_cca     <- RunPCA(seu_cca, npcs = SEURAT_NPCS, verbose = FALSE)

emb_cca <- Embeddings(seu_cca, reduction = "pca")
write.csv(emb_cca, EMBEDDING_FILES$cca, quote = FALSE)

cat("CCA complete. Saved:", EMBEDDING_FILES$cca, "| Time:", round(Sys.time() - t0, 2), "\n")
rm(anchors_cca, seu_cca, emb_cca)
gc()


cat("\nRPCA integration\n")
t0 <- Sys.time()

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, features = rownames(x), verbose = FALSE)
  return(x)
})

cat("  Finding RPCA anchors...\n")
anchors_rpca <- FindIntegrationAnchors(object.list = obj.list, dims = SEURAT_DIMS, reduction = "rpca")
seu_rpca     <- IntegrateData(anchorset = anchors_rpca, dims = SEURAT_DIMS)
seu_rpca     <- ScaleData(seu_rpca, verbose = FALSE)
seu_rpca     <- RunPCA(seu_rpca, npcs = SEURAT_NPCS, verbose = FALSE)

emb_rpca <- Embeddings(seu_rpca, reduction = "pca")
write.csv(emb_rpca, EMBEDDING_FILES$rpca, quote = FALSE)

cat("RPCA complete. Saved:", EMBEDDING_FILES$rpca, "| Time:", round(Sys.time() - t0, 2), "\n")
rm(anchors_rpca, seu_rpca, emb_rpca, obj.list, seu_master)
gc()


cat("\nReloading master object with PCA for CSS\n")
seu_master <- load_seurat(normalize = TRUE)
cat("Master object ready for CSS. Cells:", ncol(seu_master), "\n")


if (!requireNamespace("simspec", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("quadbio/simspec")
}
library(simspec)

cat("\nCSS (Pearson)\n")
t0 <- Sys.time()

seu_css_pearson  <- cluster_sim_spectrum(seu_master, label_tag = R_CFG$css_label_tag, corr_method = "pearson")
emb_css_pearson  <- Embeddings(seu_css_pearson, reduction = "css")
write.csv(emb_css_pearson, EMBEDDING_FILES$css_pearson, quote = FALSE)

cat("CSS Pearson complete. Saved:", EMBEDDING_FILES$css_pearson, "| Time:", round(Sys.time() - t0, 2), "\n")
rm(seu_css_pearson, emb_css_pearson)
gc()


cat("\nCSS (Spearman)\n")
t0 <- Sys.time()

seu_css_spearman <- cluster_sim_spectrum(seu_master, label_tag = R_CFG$css_label_tag, corr_method = "spearman")
emb_css_spearman <- Embeddings(seu_css_spearman, reduction = "css")
write.csv(emb_css_spearman, EMBEDDING_FILES$css_spearman, quote = FALSE)

cat("CSS Spearman complete. Saved:", EMBEDDING_FILES$css_spearman, "| Time:", round(Sys.time() - t0, 2), "\n")
rm(seu_css_spearman, emb_css_spearman, seu_master)
gc()


cat("\nAll R integrations complete. CSVs saved to:", CHECKPOINT_DIR, "\n")
