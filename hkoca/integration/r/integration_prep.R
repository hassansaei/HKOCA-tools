# =============================================================================
# HKOCA Integration - Prep stage (per-sample SCT + merge + resolution selection)
#
# Loads a QC-filtered Seurat RDS, splits by sample_id, runs SCTransform on each
# sample, merges into one object (layers retained), then PCA / clustree /
# silhouette / non-integrated UMAP. JoinLayers is not applied here; join RNA
# only later for FindAllMarkers if needed.
#
# Usage:
#   Rscript integration_prep.R --input_rds PATH --output_dir PATH [--config PATH]
# =============================================================================

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

.get_script_dir <- function() {
    cmd <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd, value = TRUE)
    if (length(file_arg) > 0)
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
    getwd()
}

.parse_cli_args <- function(args) {
    out <- list(); i <- 1L
    while (i <= length(args)) {
        arg <- args[[i]]
        if (startsWith(arg, "--")) {
            key <- sub("^--", "", arg); val <- "TRUE"
            if (grepl("=", key, fixed = TRUE)) {
                kv <- strsplit(key, "=", fixed = TRUE)[[1]]
                key <- kv[1]; val <- paste(kv[-1], collapse = "=")
            } else if ((i + 1L) <= length(args) && !startsWith(args[[i + 1L]], "--")) {
                val <- args[[i + 1L]]; i <- i + 1L
            }
            out[[key]] <- val
        }
        i <- i + 1L
    }
    out
}

.as_bool <- function(x, default = FALSE) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return(default)
    tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

.as_num <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") return(default)
    v <- suppressWarnings(as.numeric(x))
    if (is.na(v)) default else v
}

.resolve_path <- function(path_value, base_dir = getwd()) {
    if (is.null(path_value) || path_value == "" || path_value == "null") return(NA_character_)
    if (grepl("^~", path_value)) path_value <- path.expand(path_value)
    if (grepl("^(/|[A-Za-z]:)", path_value))
        return(normalizePath(path_value, mustWork = FALSE))
    normalizePath(file.path(base_dir, path_value), mustWork = FALSE)
}

.read_config_dcf <- function(config_path) {
    if (!file.exists(config_path))
        stop(sprintf("Config file not found: %s", config_path))
    lines <- readLines(config_path, warn = FALSE)
    cfg <- list()
    current_key <- NULL
    current_val <- ""
    flush <- function() {
        if (!is.null(current_key)) cfg[[current_key]] <<- current_val
    }
    for (ln in lines) {
        if (grepl("^\\s*#", ln) || !nzchar(trimws(ln))) next
        key_parts <- regmatches(ln, regexec("^\\s*([^:]+):\\s*(.*)$", ln))[[1]]
        if (length(key_parts) == 3) {
            flush()
            current_key <- trimws(key_parts[2])
            current_val <- key_parts[3]
        } else if (!is.null(current_key)) {
            current_val <- paste(current_val, sub("\\s+$", "", ln), sep = "\n")
        }
    }
    flush()
    cfg
}

.parse_resolutions <- function(x) {
    parts <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
    parts <- parts[nzchar(parts)]
    if (!length(parts)) stop("No resolutions configured.")
    as.numeric(parts)
}

.parse_dims <- function(x) {
    txt <- trimws(as.character(x))
    if (grepl(":", txt)) {
        bounds <- as.integer(unlist(strsplit(txt, ":", fixed = TRUE)))
        if (length(bounds) != 2L) stop(sprintf("Invalid neighbor_dims: %s", txt))
        return(seq.int(bounds[1], bounds[2]))
    }
    as.integer(unlist(strsplit(txt, ",", fixed = TRUE)))
}

.log_info <- function(...) cat(sprintf("[%s] INFO  %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

.log_warn <- function(...) cat(sprintf("[%s] WARN  %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

.log_error <- function(...) {
    cat(sprintf("[%s] ERROR %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)), file = stderr())
}

.remove_path <- function(path, label = "") {
    if (is.na(path) || !nzchar(path) || !file.exists(path)) return(invisible(FALSE))
    ok <- file.remove(path)
    if (ok) {
        .log_info("Removed intermediate%s: %s", if (nzchar(label)) sprintf(" (%s)", label) else "", path)
    } else {
        .log_warn("Could not remove intermediate: %s", path)
    }
    invisible(ok)
}

.load_integration_packages <- function() {
    required <- c(
        "Seurat", "ggplot2", "patchwork", "clustree", "glmGamPoi",
        "cluster", "digest", "jsonlite", "yaml", "dplyr", "tidyr",
        "dittoSeq", "reticulate"
    )
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing))
        stop(sprintf(
            "Missing R packages: %s\nCreate conda env: conda env create -f conda/environment_integration.yaml",
            paste(missing, collapse = ", ")
        ))
    suppressPackageStartupMessages({
        library(Seurat)
        library(ggplot2)
        library(patchwork)
        library(clustree)
        library(glmGamPoi)
        library(cluster)
        library(digest)
        library(jsonlite)
        library(yaml)
        library(dplyr)
        library(tidyr)
        library(dittoSeq)
        library(reticulate)
    })
}

.load_hkoca_color_palettes <- function(script_dir, cfg) {
    colors_yaml <- .resolve_path(cfg$celltype_colors_yaml)
    if (is.na(colors_yaml) || !nzchar(colors_yaml)) {
        colors_yaml <- Sys.getenv("HKOCA_CELLTYPE_COLORS", unset = "")
    }
    if (!nzchar(colors_yaml)) {
        colors_yaml <- normalizePath(
            file.path(script_dir, "..", "..", "config", "celltype_colors.yaml"),
            mustWork = FALSE
        )
    }
    colors_r <- file.path(script_dir, "hkoca_colors.R")
    if (!file.exists(colors_r))
        stop(sprintf("Missing HKOCA color helpers: %s", colors_r))
    source(colors_r, local = FALSE)
    palettes <- load_hkoca_celltype_palettes(colors_yaml)
    .log_info("Loaded HKOCA cell-type colors: %s", palettes$path)
    viz_r <- file.path(script_dir, "integration_viz.R")
    if (!file.exists(viz_r))
        stop(sprintf("Missing integration visualization helpers: %s", viz_r))
    source(viz_r, local = FALSE)
    palettes
}

.ensure_mito_percent <- function(obj, mito_col) {
    if (mito_col %in% colnames(obj@meta.data)) return(obj)
    .log_info("Column '%s' missing; computing from MT- genes.", mito_col)
    obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    if (mito_col != "percent.mito") {
        obj[[mito_col]] <- obj[["percent.mito"]]
    }
    obj
}

.strip_barcode_suffix <- function(barcodes) {
    sub("-[0-9]+$", "", as.character(barcodes))
}

.transfer_annotation_from_h5ad <- function(obj, h5ad_path) {
    wanted <- c("Level_1", "Level_2", "Level_3")
    already <- wanted[wanted %in% colnames(obj@meta.data)]
    if (length(already) == length(wanted)) {
        .log_info("Level_1/2/3 already present in Seurat metadata; skipping h5ad transfer.")
        return(obj)
    }
    if (is.na(h5ad_path) || !nzchar(h5ad_path) || !file.exists(h5ad_path)) {
        .log_warn("No annotated h5ad provided; Level proportion plots may be skipped.")
        return(obj)
    }

    .log_info("Transferring Level_1/2/3 from annotated h5ad: %s", h5ad_path)
    ann_py <- Sys.getenv("HKOCA_ANNOTATION_PYTHON", unset = "")
    if (nzchar(ann_py) && file.exists(ann_py)) {
        tryCatch(reticulate::use_python(ann_py, required = FALSE), error = function(e) NULL)
    }
    if (!reticulate::py_module_available("scanpy") && !reticulate::py_module_available("anndata")) {
        .log_warn("Python scanpy/anndata not available via reticulate; cannot transfer annotation.")
        return(obj)
    }

    obs_df <- tryCatch({
        if (reticulate::py_module_available("scanpy")) {
            sc <- reticulate::import("scanpy", convert = FALSE)
            ad <- sc$read_h5ad(h5ad_path)
        } else {
            anndata <- reticulate::import("anndata", convert = FALSE)
            ad <- anndata$read_h5ad(h5ad_path)
        }
        reticulate::py_to_r(ad$obs)
    }, error = function(e) {
        .log_warn("Failed to read annotated h5ad obs: %s", conditionMessage(e))
        NULL
    })
    if (is.null(obs_df) || !nrow(obs_df)) return(obj)

    avail <- intersect(wanted, colnames(obs_df))
    if (!length(avail)) {
        .log_warn("Annotated h5ad has no Level_1/2/3 columns.")
        return(obj)
    }

    seurat_cells <- colnames(obj)
    h5ad_cells <- rownames(obs_df)
    map_from <- seurat_cells
    map_to <- match(seurat_cells, h5ad_cells)

    if (sum(!is.na(map_to)) < 0.5 * length(seurat_cells)) {
        # Try stripped barcode match (Seurat -1 suffixes / h5ad bare barcodes)
        map_to <- match(.strip_barcode_suffix(seurat_cells), .strip_barcode_suffix(h5ad_cells))
    }
    if (sum(!is.na(map_to)) < 0.5 * length(seurat_cells)) {
        # Try matching on trailing barcode token after '_'
        seurat_tail <- sub("^.*_", "", seurat_cells)
        h5ad_tail <- sub("^.*_", "", h5ad_cells)
        map_to <- match(seurat_tail, h5ad_tail)
    }

    n_matched <- sum(!is.na(map_to))
    .log_info("Matched %d / %d cells to annotated h5ad.", n_matched, length(seurat_cells))
    if (n_matched == 0L) {
        .log_warn("No overlapping cell barcodes between RDS and annotated h5ad.")
        return(obj)
    }

    for (col in avail) {
        vals <- rep(NA_character_, length(seurat_cells))
        vals[!is.na(map_to)] <- as.character(obs_df[[col]][map_to[!is.na(map_to)]])
        obj[[col]] <- vals
        .log_info("  Added metadata column '%s' (%d non-missing).", col, sum(!is.na(vals) & nzchar(vals)))
    }
    obj
}

.count_batch_levels <- function(obj, batch_col) {
    vals <- as.character(obj@meta.data[[batch_col]])
    length(unique(vals[nzchar(vals) & !is.na(vals)]))
}

.detect_batch_column <- function(obj, preferred = NULL) {
    if (!is.null(preferred) && nzchar(preferred)) {
        if (!preferred %in% colnames(obj@meta.data)) {
            stop(sprintf("Configured integration_batch_key '%s' not found in metadata.", preferred))
        }
        return(preferred)
    }
    candidates <- c("sample_id", "orig.ident", "study", "batch", "dataset")
    for (col in candidates) {
        if (col %in% colnames(obj@meta.data) && .count_batch_levels(obj, col) >= 2L) {
            return(col)
        }
    }
    for (col in candidates) {
        if (col %in% colnames(obj@meta.data)) return(col)
    }
    stop("No batch column found in metadata (expected sample_id, orig.ident, or study).")
}

.sct_one_sample <- function(x, mito_col, sample_label) {
    .log_info(
        "SCTransform sample '%s' (%s cells).",
        sample_label, format(ncol(x), big.mark = ",")
    )
    DefaultAssay(x) <- "RNA"
    if ("SCT" %in% Assays(x)) {
        DefaultAssay(x) <- "RNA"
        x[["SCT"]] <- NULL
    }
    assay_obj <- tryCatch(x[["RNA"]], error = function(e) NULL)
    if (!is.null(assay_obj) && inherits(assay_obj, "Assay5") && exists("JoinLayers", mode = "function")) {
        n_layers <- tryCatch(length(Layers(assay_obj)), error = function(e) 1L)
        if (n_layers > 1L) x[["RNA"]] <- JoinLayers(assay_obj)
    }
    SCTransform(
        x,
        assay = "RNA",
        method = "glmGamPoi",
        vars.to.regress = mito_col,
        vst.flavor = "v2",
        verbose = FALSE
    )
}

.split_sct_merge <- function(obj, batch_col, mito_col, n_features) {
    n_batches <- .count_batch_levels(obj, batch_col)
    .log_info("Integration batch column: '%s' (%d level(s)).", batch_col, n_batches)

    if (n_batches < 2L) {
        .log_info("Single batch; running one SCTransform (no merge).")
        obj <- .sct_one_sample(obj, mito_col, "all")
        DefaultAssay(obj) <- "SCT"
        obj <- FindVariableFeatures(
            obj, selection.method = "vst", assay = "SCT",
            nfeatures = n_features, verbose = FALSE
        )
        return(obj)
    }

    .log_info("SplitObject by '%s' for per-sample SCTransform.", batch_col)
    obj_list <- SplitObject(obj, split.by = batch_col)
    sample_names <- names(obj_list)
    if (is.null(sample_names) || !length(sample_names)) {
        sample_names <- as.character(seq_along(obj_list))
    }

    obj_list <- lapply(seq_along(obj_list), function(i) {
        .sct_one_sample(obj_list[[i]], mito_col, sample_names[[i]])
    })
    names(obj_list) <- sample_names

    .log_info("Merging %d SCT-normalized samples (layers retained).", length(obj_list))
    obj <- merge(x = obj_list[[1]], y = obj_list[-1], project = "hkoca_integration")
    DefaultAssay(obj) <- "SCT"

    .log_info("Selecting integration variable features (n=%d).", n_features)
    features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = n_features)
    VariableFeatures(obj) <- features
    .log_info(
        "Merged object: %s cells, %d variable features.",
        format(ncol(obj), big.mark = ","), length(features)
    )
    obj
}

.run_nonintegrated_umap_and_plots <- function(obj, out_dir, dims, cluster_col,
                                              markers_yaml_path, palettes = NULL,
                                              resolution = NA_real_) {
    .log_info("Running non-integrated UMAP from PCA.")
    obj <- RunUMAP(
        obj,
        dims = dims,
        reduction = "pca",
        reduction.name = "umap_unintegrated",
        reduction.key = "UMAPUN_",
        verbose = FALSE
    )
    vis <- .save_integration_figures(
        obj, out_dir,
        reduction = "umap_unintegrated",
        cluster_col = cluster_col,
        markers_yaml_path = markers_yaml_path,
        method_label = "Non-integrated",
        feature_assay = "RNA",
        palettes = palettes,
        resolution = resolution
    )
    list(object = obj, outputs = vis)
}

.calc_silhouette <- function(obj, res, prefix, dims, n_subsample = 5000, seed = 42L) {
    Idents(obj) <- paste0(prefix, res)
    set.seed(seed)
    cells_use <- sample(colnames(obj), min(n_subsample, ncol(obj)))
    obj_sub <- subset(obj, cells = cells_use)
    clusters <- as.integer(Idents(obj_sub))
    if (length(unique(clusters)) < 2L) return(NA_real_)
    pca_embed <- Embeddings(obj_sub, "pca")[, dims, drop = FALSE]
    dist_mat <- dist(pca_embed)
    sil <- silhouette(clusters, dist_mat)
    mean(sil[, 3])
}

.print_usage <- function(default_config) {
    cat("Usage:\n")
    cat("  Rscript integration_prep.R --input_rds PATH --output_dir PATH [OPTIONS]\n\n")
    cat("Options:\n")
    cat(sprintf("  --config PATH            Config DCF (default: %s)\n", default_config))
    cat("  --input_rds PATH         QC-filtered Seurat RDS (singlets filtered)\n")
    cat("  --annotated_h5ad PATH    Annotated h5ad from hkoca annotation (optional)\n")
    cat("  --output_dir PATH        Output directory\n")
    cat("  --force_overwrite        Re-run even if outputs exist\n")
    cat("  --remove_intermediate    Delete heavy intermediate files after success\n")
    cat("  --help                   Show this message\n")
}

# =============================================================================
# Main
# =============================================================================

tryCatch({
    script_dir <- .get_script_dir()
    default_config <- file.path(script_dir, "integration.config.dcf")
    cli_args <- .parse_cli_args(commandArgs(trailingOnly = TRUE))

    if (.as_bool(cli_args$help, FALSE)) {
        .print_usage(default_config)
        quit(save = "no", status = 0)
    }

    config_path <- if (!is.null(cli_args$config)) .resolve_path(cli_args$config) else default_config
    cfg <- .read_config_dcf(config_path)
    for (nm in setdiff(names(cli_args), c("config", "help")))
        cfg[[nm]] <- cli_args[[nm]]

    r_max_vsize <- cfg$r_max_vsize %||% "150Gb"
    Sys.setenv(R_MAX_VSIZE = r_max_vsize)
    .log_info("Set R_MAX_VSIZE=%s", r_max_vsize)

    input_rds <- .resolve_path(cfg$input_rds)
    annotated_h5ad <- .resolve_path(cfg$annotated_h5ad)
    output_dir <- .resolve_path(cfg$output_dir %||% "integration_results")
    force_flag <- .as_bool(cfg$force_overwrite, FALSE)

    if (is.na(input_rds) || !nzchar(input_rds) || !file.exists(input_rds))
        stop("Provide --input_rds pointing to an existing QC-filtered Seurat RDS.")

    prep_dir <- file.path(output_dir, "prep")
    nonint_dir <- file.path(output_dir, "nonintegrated")
    tab_dir  <- file.path(output_dir, "tables")
    for (d in c(prep_dir, nonint_dir, tab_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    markers_yaml <- .resolve_path(cfg$markers_yaml)
    if (is.na(markers_yaml) || !nzchar(markers_yaml)) {
        markers_yaml <- Sys.getenv("HKOCA_MARKERS_YAML", unset = "")
    }
    if (!nzchar(markers_yaml)) {
        markers_yaml <- normalizePath(
            file.path(script_dir, "..", "..", "config", "snapseed_markers_v4.yaml"),
            mustWork = FALSE
        )
    }
    if (file.exists(markers_yaml)) {
        .log_info("Using markers YAML: %s", markers_yaml)
    } else {
        .log_info("Markers YAML not found, FeaturePlots will be skipped: %s", markers_yaml)
        markers_yaml <- ""
    }

    log_file <- file.path(output_dir, "integration_prep.log")
    sink(log_file, split = TRUE)
    on.exit(sink(), add = TRUE)

    prepared_rds <- file.path(prep_dir, "sct_prepared.rds")
    meta_json <- file.path(prep_dir, ".sct_prepared.meta.json")

    n_features <- as.integer(cfg$n_features %||% 2500L)
    npcs <- as.integer(cfg$npcs %||% 50L)
    neighbor_dims <- .parse_dims(cfg$neighbor_dims %||% "1:30")
    mito_col <- cfg$mito_regress %||% "percent.mito"
    batch_key <- cfg$integration_batch_key %||% "sample_id"
    resolutions <- .parse_resolutions(cfg$resolutions)
    cluster_prefix <- cfg$cluster_prefix %||% "SCT_snn_res."
    sil_subsample <- as.integer(cfg$silhouette_subsample %||% 5000L)
    sil_seed <- as.integer(cfg$silhouette_seed %||% 42L)
    sil_min_resolution <- as.numeric(cfg$silhouette_min_resolution %||% 0.4)
    sil_min_clusters <- as.integer(cfg$silhouette_min_clusters %||% 5L)

    run_params <- list(
        input_rds = input_rds,
        annotated_h5ad = annotated_h5ad,
        n_features = n_features,
        npcs = npcs,
        neighbor_dims = neighbor_dims,
        mito_regress = mito_col,
        integration_batch_key = batch_key,
        resolutions = resolutions,
        cluster_prefix = cluster_prefix,
        silhouette_subsample = sil_subsample,
        silhouette_min_resolution = sil_min_resolution,
        silhouette_min_clusters = sil_min_clusters,
        prep_workflow = "split_sct_merge"
    )

    if (!force_flag && file.exists(prepared_rds) && file.info(prepared_rds)$size > 0 &&
        file.exists(meta_json)) {
        old <- tryCatch(jsonlite::fromJSON(meta_json), error = function(e) NULL)
        if (!is.null(old) && identical(old$param_hash, digest::digest(run_params))) {
            .log_info("Skipping prep; existing output unchanged: %s", prepared_rds)
            quit(save = "no", status = 0)
        }
    }

    .load_integration_packages()
    hkoca_palettes <- .load_hkoca_color_palettes(script_dir, cfg)

    .log_info("Loading Seurat object: %s", input_rds)
    obj <- readRDS(input_rds)
    if (!inherits(obj, "Seurat")) stop("input_rds must be a Seurat object.")
    .log_info("Loaded %s cells x %s features", format(ncol(obj), big.mark = ","), format(nrow(obj), big.mark = ","))

    if (!is.na(annotated_h5ad) && nzchar(annotated_h5ad)) {
        if (file.exists(annotated_h5ad)) {
            .log_info("Annotated h5ad: %s", annotated_h5ad)
        } else {
            .log_info("Annotated h5ad path not found (continuing prep): %s", annotated_h5ad)
            annotated_h5ad <- NA_character_
        }
    }

    # Transfer Level labels before SplitObject/merge so barcodes still match
    obj <- .transfer_annotation_from_h5ad(obj, annotated_h5ad)
    obj <- .ensure_mito_percent(obj, mito_col)
    batch_col <- .detect_batch_column(obj, batch_key)

    # Per-sample SCT then merge (do not JoinLayers on the merged object here)
    obj <- .split_sct_merge(obj, batch_col, mito_col, n_features)
    DefaultAssay(obj) <- "SCT"

    .log_info("Running PCA (npcs=%d).", npcs)
    obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)

    elbow_path <- file.path(prep_dir, "elbow_plot.png")
    .log_info("Saving elbow plot: %s", elbow_path)
    .save_ggplot_png(elbow_path, ElbowPlot(obj, ndims = npcs), width = 8, height = 5)

    .log_info("Finding neighbors (dims %s).", paste(range(neighbor_dims), collapse = ":"))
    obj <- FindNeighbors(obj, dims = neighbor_dims, reduction = "pca", verbose = FALSE)

    .log_info("Testing resolutions: %s", paste(resolutions, collapse = ", "))
    obj <- FindClusters(obj, resolution = resolutions, verbose = FALSE)

    clustree_path <- file.path(prep_dir, "clustree_sct.png")
    .log_info("Saving clustree plot: %s", clustree_path)
    .save_ggplot_png(clustree_path, clustree(obj, prefix = cluster_prefix), width = 12, height = 10)

    .log_info("Calculating silhouette scores (subsample=%d).", sil_subsample)
    sil_scores <- vapply(resolutions, function(r) {
        .log_info("  Testing resolution: %s", r)
        .calc_silhouette(obj, r, prefix = cluster_prefix, dims = neighbor_dims,
                         n_subsample = sil_subsample, seed = sil_seed)
    }, numeric(1))

    sil_df <- data.frame(resolution = resolutions, silhouette = sil_scores)
    sil_csv <- file.path(tab_dir, "silhouette_scores.csv")
    write.csv(sil_df, sil_csv, row.names = FALSE)
    .log_info("Saved silhouette table: %s", sil_csv)

    best_pick <- .select_best_silhouette_resolution(
        sil_df,
        obj,
        cluster_prefix,
        min_resolution = sil_min_resolution,
        min_clusters = sil_min_clusters
    )
    best_res <- best_pick$resolution
    if (!is.na(best_res)) {
        .log_info(
            "Recommended resolution (silhouette, >= %.2f, min %d clusters): %s (score=%.4f)",
            sil_min_resolution, sil_min_clusters, best_res, best_pick$silhouette
        )
    }

    sil_plot_path <- file.path(prep_dir, "silhouette_sct.png")
    p_sil <- ggplot(sil_df, aes(x = resolution, y = silhouette)) +
        geom_point(size = 3) +
        geom_line() +
        theme_minimal() +
        labs(
            title = sprintf(
                "Silhouette score vs resolution (best >= %.2f)",
                sil_min_resolution
            ),
            x = "Resolution",
            y = "Silhouette"
        )
    if (!is.na(best_res)) {
        p_sil <- p_sil +
            geom_vline(xintercept = best_res, color = "red", linetype = "dashed")
    }
    .save_ggplot_png(sil_plot_path, p_sil, width = 8, height = 5)
    .log_info("Saved silhouette plot: %s", sil_plot_path)

    summary_df <- data.frame(
        best_resolution = best_res,
        best_silhouette = best_pick$silhouette,
        silhouette_min_resolution = sil_min_resolution,
        silhouette_min_clusters = sil_min_clusters,
        npcs = npcs,
        neighbor_dims = paste(range(neighbor_dims), collapse = ":"),
        n_features = n_features,
        stringsAsFactors = FALSE
    )
    summary_csv <- file.path(tab_dir, "resolution_summary.csv")
    write.csv(summary_df, summary_csv, row.names = FALSE)
    .log_info("Recommended resolution (silhouette): %s", best_res)

    cluster_col <- paste0(cluster_prefix, best_res)
    if (!cluster_col %in% colnames(obj@meta.data)) {
        stop(sprintf("Selected cluster column not found in metadata: %s", cluster_col))
    }
    nonint_result <- .run_nonintegrated_umap_and_plots(
        obj, nonint_dir, neighbor_dims, cluster_col, markers_yaml,
        palettes = hkoca_palettes,
        resolution = best_res
    )
    obj <- nonint_result$object
    nonintegrated_outputs <- nonint_result$outputs

    .log_info("Saving prepared object: %s", prepared_rds)
    saveRDS(obj, prepared_rds)

    meta_out <- list(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        param_hash = digest::digest(run_params),
        params = run_params,
        celltype_colors_yaml = hkoca_palettes$path,
        outputs = list(
            prepared_rds = prepared_rds,
            elbow_plot = elbow_path,
            clustree_plot = clustree_path,
            silhouette_plot = sil_plot_path,
            nonintegrated_dir = nonint_dir,
            nonintegrated_umap = nonintegrated_outputs$overall,
            nonintegrated_split_columns = nonintegrated_outputs$split_columns,
            nonintegrated_feature_plots = length(nonintegrated_outputs$feature_plots),
            nonintegrated_proportion_plots = nonintegrated_outputs$proportion_plots,
            silhouette_csv = sil_csv,
            resolution_summary = summary_csv
        ),
        best_resolution = best_res
    )
    jsonlite::write_json(meta_out, meta_json, auto_unbox = TRUE, pretty = TRUE)

    if (.as_bool(cfg$remove_intermediate, FALSE)) {
        input_in_output <- tryCatch(
            startsWith(
                normalizePath(input_rds, mustWork = FALSE),
                normalizePath(output_dir, mustWork = FALSE)
            ),
            error = function(e) FALSE
        )
        if (input_in_output) {
            .remove_path(input_rds, "input RDS copy under output_dir")
        }
    }

    .log_info("Integration prep completed. Outputs under: %s", output_dir)
    quit(save = "no", status = 0)
}, error = function(e) {
    .log_error("%s", conditionMessage(e))
    quit(save = "no", status = 1)
})
