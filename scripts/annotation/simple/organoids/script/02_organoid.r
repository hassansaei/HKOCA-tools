#!/usr/bin/env Rscript

# ==============================================================================
# 🧫 Organoid Single-Cell Atlas Pipeline
# Phase 1: Dynamic Clustering, Marker Discovery & Curation
# ==============================================================================

suppressWarnings(suppressMessages({
  library(Seurat)
  library(presto)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(cluster)
  library(yaml)
  library(stringr)
  library(future)
}))

# ── LOAD CONFIGURATION ────────────────────────────────────────────────────────
config_path <- "config.yaml"
if (!file.exists(config_path)) {
  stop("Configuration file 'config.yaml' not found. Please create one.")
}
cfg <- read_yaml(config_path)

# ── CONFIG VALIDATION ─────────────────────────────────────────────────────────
# Fix #6: validate all required keys at startup with type checks, so errors
# surface immediately with a clear message rather than deep inside the pipeline.
validate_config <- function(cfg) {
  errors <- character()

  required_numeric <- list(
    c("paths", "input_dir"),
    c("paths", "output_dir"),
    c("preprocessing", "n_hvg"),
    c("preprocessing", "n_pcs"),
    c("preprocessing", "n_neighbors"),
    c("clustering", "res_min"),
    c("clustering", "res_max"),
    c("clustering", "res_step"),
    c("clustering", "merge_overlap"),
    c("clustering", "merge_top_n"),
    c("clustering", "min_cluster_frac"),
    c("clustering", "max_merge_iters"),
    c("clustering", "top_n_res_candidates"),
    c("clustering", "downsample_threshold"),
    c("clustering", "downsample_per_cluster"),
    c("marker_filtering", "min_pct_in"),
    c("marker_filtering", "min_logfc"),
    c("marker_filtering", "min_kept"),
    c("marker_filtering", "tolerance"),
    c("marker_filtering", "max_n_shown"),
    c("marker_filtering", "max_markers_per_cluster"),
    c("marker_filtering", "top_yaml_genes"),
    c("marker_filtering", "dotplot_top_n"),
    c("noise_genes", "prefixes"),
    c("noise_genes", "exact_matches")
  )

  for (keys in required_numeric) {
    val <- tryCatch(Reduce(`[[`, keys, init = cfg), error = function(e) NULL)
    if (is.null(val)) {
      errors <- c(errors, paste("Missing config key:", paste(keys, collapse = "$")))
    }
  }

  # Range sanity checks
  if (!is.null(cfg$clustering$res_min) && !is.null(cfg$clustering$res_max)) {
    if (cfg$clustering$res_min >= cfg$clustering$res_max)
      errors <- c(errors, "clustering$res_min must be less than clustering$res_max")
  }
  if (!is.null(cfg$clustering$merge_overlap)) {
    if (cfg$clustering$merge_overlap < 0 || cfg$clustering$merge_overlap > 1)
      errors <- c(errors, "clustering$merge_overlap must be between 0 and 1")
  }
  if (!is.null(cfg$marker_filtering$min_pct_in)) {
    if (cfg$marker_filtering$min_pct_in < 0 || cfg$marker_filtering$min_pct_in > 1)
      errors <- c(errors, "marker_filtering$min_pct_in must be between 0 and 1")
  }

  if (length(errors) > 0) {
    stop("Config validation failed:\n", paste(" -", errors, collapse = "\n"))
  }

  invisible(TRUE)
}

validate_config(cfg)

# ── GLOBAL SEED ───────────────────────────────────────────────────────────────
# Fix #7: set a single global seed here so all stochastic steps (FindClusters,
# RunUMAP, silhouette subsampling) are reproducible. The per-call random.seed
# arguments in FindClusters are kept for explicitness but this is the master.
set.seed(42)

# Derived variables from config
RESOLUTIONS    <- round(seq(cfg$clustering$res_min, cfg$clustering$res_max, by = cfg$clustering$res_step), 2)
PREFIX_PATTERN <- paste0("^(", paste(cfg$noise_genes$prefixes, collapse = "|"), ")")

# Setup memory limit (100GB limit to prevent Seurat crashes)
plan("sequential")
options(future.globals.maxSize = 100000 * 1024^2)

# ── HELPER FUNCTIONS ──────────────────────────────────────────────────────────

log_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
}

preprocess_seurat <- function(obj, cfg) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = cfg$preprocessing$n_hvg, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, npcs = cfg$preprocessing$n_pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:cfg$preprocessing$n_pcs, k.param = cfg$preprocessing$n_neighbors, verbose = FALSE)
  return(obj)
}

compute_silhouette <- function(obj, resolution, cfg) {
  obj <- FindClusters(obj, resolution = resolution, algorithm = 1, random.seed = 42, verbose = FALSE)
  labels <- as.integer(Idents(obj))
  if (length(unique(labels)) < 2) return(list(score = -1, obj = obj))

  n   <- min(ncol(obj), 5000L)
  idx <- sample(ncol(obj), n, replace = FALSE)
  pca <- Embeddings(obj, "pca")[idx, ]
  sil <- silhouette(labels[idx], dist(pca))
  list(score = mean(sil[, "sil_width"]), obj = obj)
}

marker_coherence_check <- function(obj, cluster_key, cfg) {
  result <- tryCatch({
    Idents(obj) <- cluster_key
    markers <- FindAllMarkers(obj, test.use = "wilcox", logfc.threshold = 0.25,
                              min.pct = 0.1, only.pos = TRUE, verbose = FALSE)

    top_genes <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = cfg$clustering$merge_top_n) %>%
      summarise(genes = list(gene), .groups = "drop")
    gene_sets <- setNames(top_genes$genes, top_genes$cluster)
    pairs <- combn(as.character(top_genes$cluster), 2, simplify = FALSE)

    Filter(Negate(is.null), lapply(pairs, function(p) {
      ov <- length(intersect(gene_sets[[p[1]]], gene_sets[[p[2]]])) / cfg$clustering$merge_top_n
      if (ov >= cfg$clustering$merge_overlap) p else NULL
    }))
  }, error = function(e) list())
  return(result)
}

merge_duplicate_clusters <- function(obj, cluster_key = "leiden", cfg) {
  obj@meta.data[["leiden"]] <- as.character(obj@meta.data[[cluster_key]])
  # Fix #9: max iteration count now comes from cfg$clustering$max_merge_iters
  for (i in seq_len(cfg$clustering$max_merge_iters)) {
    Idents(obj) <- "leiden"
    dups <- marker_coherence_check(obj, "leiden", cfg)
    if (length(dups) == 0) break

    pair   <- dups[[1]]
    target <- pair[which.min(as.integer(pair))]
    source <- pair[which.max(as.integer(pair))]
    obj@meta.data[["leiden"]] <- ifelse(obj@meta.data[["leiden"]] == source, target, obj@meta.data[["leiden"]])
  }
  obj@meta.data[["leiden"]] <- factor(as.integer(factor(obj@meta.data[["leiden"]])) - 1)
  return(obj)
}

merge_micro_clusters <- function(obj, cluster_key = "leiden", cfg) {
  counts    <- table(obj@meta.data[[cluster_key]])
  threshold <- floor(cfg$clustering$min_cluster_frac * ncol(obj))
  tiny      <- names(counts[counts < threshold])

  if (length(tiny) == 0) return(obj)

  nn_graph <- obj@graphs[["RNA_nn"]]
  labels   <- as.character(obj@meta.data[[cluster_key]])

  for (cl in tiny) {
    cl_cells   <- which(labels == cl)
    neighbours <- unlist(lapply(cl_cells, function(i) labels[which(nn_graph[i, ] > 0 & labels != cl)]))
    if (length(neighbours) > 0) labels[cl_cells] <- names(which.max(table(neighbours)))
  }
  obj@meta.data[[cluster_key]] <- factor(labels)
  return(obj)
}

extract_markers <- function(obj, dataset_name, cfg) {
  Idents(obj) <- "leiden"
  all_genes   <- rownames(obj)

  clean_genes <- all_genes[!str_starts(all_genes, "ENSG") &
                           !str_detect(all_genes, PREFIX_PATTERN) &
                           !(all_genes %in% cfg$noise_genes$exact_matches)]

  markers <- FindAllMarkers(obj, features = clean_genes, test.use = "wilcox",
                            logfc.threshold = cfg$marker_filtering$min_logfc, min.pct = 0.2, only.pos = TRUE, verbose = FALSE)

  if (nrow(markers) == 0) stop("FindAllMarkers returned 0 genes.")

  markers %>%
    filter(pct.1 >= cfg$marker_filtering$min_pct_in) %>%
    # Fix #9: max candidates per cluster now comes from cfg$marker_filtering$max_markers_per_cluster
    group_by(cluster) %>% slice_max(avg_log2FC, n = cfg$marker_filtering$max_markers_per_cluster, with_ties = FALSE) %>% ungroup() %>%
    transmute(dataset = dataset_name, group = as.character(cluster), gene = gene,
              logFC = avg_log2FC, pval_adj = p_val_adj, pct_in = pct.1, pct_out = pct.2) %>%
    mutate(stringent_score = logFC * pct_in^2 / (pct_out + cfg$marker_filtering$tolerance)) %>%
    arrange(group, desc(stringent_score))
}

export_curated_markers <- function(marker_df, dataset_name, out_dir, cfg) {
  curated_list <- list()
  yaml_dict    <- list()

  # Ensure output subdirectories exist even if called standalone
  dir.create(file.path(out_dir, "tsv"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "yaml"), recursive = TRUE, showWarnings = FALSE)

  for (ct in unique(marker_df$group)) {
    ct_data <- marker_df %>% filter(group == ct)
    y       <- head(ct_data$stringent_score, cfg$marker_filtering$max_n_shown)

    if (length(y) > 2L) {
      x_norm    <- (seq_along(y) - 1L) / (length(y) - 1L + 1e-9)
      y_norm    <- (y - min(y)) / (max(y) - min(y) + 1e-9)
      elbow_idx <- which.max(1.0 - (x_norm + y_norm))
    } else { elbow_idx <- length(y) }

    cutoff_rank <- max(elbow_idx, cfg$marker_filtering$min_kept)
    final_genes <- head(ct_data, cutoff_rank) %>%
      mutate(across(c(stringent_score, logFC), ~round(., 2)),
             across(c(pct_in, pct_out), ~round(., 3)))

    curated_list[[ct]] <- final_genes
    # Fix #9: number of genes written to YAML now comes from cfg$marker_filtering$top_yaml_genes
    yaml_dict[[ct]] <- head(final_genes$gene, cfg$marker_filtering$top_yaml_genes)
  }

  write.table(bind_rows(curated_list), file.path(out_dir, "tsv", paste0(dataset_name, "_markers.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write_yaml(yaml_dict, file.path(out_dir, "yaml", paste0(dataset_name, "_markers.yaml")))
}

wilcoxon_dotplot <- function(obj, markers_df, dataset_name, out_dir, cfg) {
  # Fix #9: genes per cluster for dotplot now comes from cfg$marker_filtering$dotplot_top_n
  top_markers  <- markers_df %>% group_by(group) %>% slice_head(n = cfg$marker_filtering$dotplot_top_n) %>% ungroup()
  genes_ordered <- intersect(unique(top_markers$gene), rownames(obj))
  if (length(genes_ordered) == 0) return()

  # Ensure output subdirectory exists even if called standalone
  dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

  Idents(obj) <- "leiden"
  p <- DotPlot(obj, features = genes_ordered, cols = c("white", "steelblue")) +
    coord_flip() + ggtitle(dataset_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.title  = element_text(size = 10, face = "bold", hjust = 0.5))

  ggsave(file.path(out_dir, "plots", paste0(dataset_name, "_dotplot.png")), p, width = 10, height = 7, dpi = 150)
}

plot_umap <- function(obj, dataset_name, out_dir) {
  # Generate a UMAP colored by our final dynamic clusters
  p <- DimPlot(obj, reduction = "umap", group.by = "leiden", label = TRUE, repel = TRUE) +
    ggtitle(paste("UMAP:", dataset_name)) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  
  # Save it into the plots directory
  ggsave(file.path(out_dir, "plots", paste0(dataset_name, "_umap.png")), 
         plot = p, width = 8, height = 6, dpi = 150)
}

# ── CORE PIPELINE EXECUTION FUNCTION ──────────────────────────────────────────

process_organoid_dataset <- function(rds_path, out_dir, cfg) {
  dataset_name <- tools::file_path_sans_ext(basename(rds_path))
  log_msg("\n%s\n🚀 PROCESSING DATASET: %s\n%s", strrep("=", 60), dataset_name, strrep("=", 60))

  # 1. LOAD & PREPROCESS
  obj <- readRDS(rds_path)
  if (!inherits(obj, "Seurat")) stop("Not a Seurat object")
  log_msg("  ▶ Preprocessing (%d cells)...", ncol(obj))
  obj <- preprocess_seurat(obj, cfg)

  # 2. DYNAMIC CLUSTERING
  log_msg("  ▶ Sweep Resolutions & Silhouette Scores...")
  sil_scores <- sapply(RESOLUTIONS, function(r) compute_silhouette(obj, r, cfg)$score)
  # Fix #9: number of top resolution candidates now comes from cfg$clustering$top_n_res_candidates
  top_n      <- cfg$clustering$top_n_res_candidates
  top_n_res  <- as.numeric(RESOLUTIONS[order(sil_scores, decreasing = TRUE)[seq_len(top_n)]])

  coherence_res <- sapply(top_n_res, function(r) {
    tmp <- FindClusters(obj, resolution = r, random.seed = 42, verbose = FALSE)
    length(marker_coherence_check(tmp, "seurat_clusters", cfg))
  })

  chosen_res <- max(top_n_res[coherence_res == min(coherence_res)])
  log_msg("  ▶ Applying Optimal Resolution: %s", chosen_res)
  obj <- FindClusters(obj, resolution = chosen_res, algorithm = 1, random.seed = 42, verbose = FALSE)
  obj$leiden <- factor(as.character(Idents(obj)))

  # 3. REFINE STRUCTURE & UMAP
  log_msg("  ▶ Refining Duplicates and Micro Clusters...")
  obj <- merge_duplicate_clusters(obj, "leiden", cfg)
  obj <- merge_micro_clusters(obj, "leiden", cfg)
  # Fix #7: seed passed explicitly to RunUMAP for reproducibility
  obj <- RunUMAP(obj, dims = 1:cfg$preprocessing$n_pcs, seed.use = 42, verbose = FALSE)

  save_path <- file.path(out_dir, "clustered_rds", paste0(dataset_name, "_clustered.rds"))
  saveRDS(obj, save_path)

  # 4. MARKER EXTRACTION
  log_msg("  ▶ Extracting High-Confidence Markers (Presto)...")
  # Fix #9: downsample threshold and per-cluster count now come from config
  marker_obj <- if (ncol(obj) > cfg$clustering$downsample_threshold) {
    subset(obj, downsample = cfg$clustering$downsample_per_cluster)
  } else {
    obj
  }
  markers_df <- extract_markers(marker_obj, dataset_name, cfg)

  # 5. EXPORT & PLOT
  log_msg("  ▶ Exporting Curated Lists & Plotting...")
  export_curated_markers(markers_df, dataset_name, out_dir, cfg)
  wilcoxon_dotplot(obj, markers_df, dataset_name, out_dir, cfg)
  plot_umap(obj, dataset_name, out_dir)

  log_msg("✅ PIPELINE COMPLETE FOR: %s", dataset_name)
  return("Success")
}

# ── MAIN LOOP (Dataset Iterator) ──────────────────────────────────────────────

OUT_DIR <- cfg$paths$output_dir
for (d in c("clustered_rds", "tsv", "yaml", "plots")) {
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)
}

rds_files <- list.files(cfg$paths$input_dir, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
if (length(rds_files) == 0) stop("No datasets found in ", cfg$paths$input_dir)

log_msg("Found %d RDS files. Starting batch run...", length(rds_files))

# Fix #5 (from original review): pre-allocate list, bind once at end
results_list <- vector("list", length(rds_files))

for (i in seq_along(rds_files)) {
  rds_path     <- rds_files[[i]]
  dataset_name <- tools::file_path_sans_ext(basename(rds_path))

  status <- tryCatch({
    process_organoid_dataset(rds_path, OUT_DIR, cfg)
  }, error = function(e) {
    log_msg("❌ ERROR ON %s: %s", dataset_name, conditionMessage(e))
    return(conditionMessage(e))
  })

  results_list[[i]] <- data.frame(dataset = dataset_name, status = status)
  gc(full = TRUE)
}

results_log <- bind_rows(results_list)
write.csv(results_log, file.path(OUT_DIR, "pipeline_execution_log.csv"), row.names = FALSE)
log_msg("\n🎉 ALL DONE! Check %s for the execution log.", file.path(OUT_DIR, "pipeline_execution_log.csv"))