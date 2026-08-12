# =============================================================================
# HKOCA Integration - Prep stage (SCT + resolution selection)
#
# Loads a QC-filtered Seurat RDS (singlets filtered), normalizes RNA, runs
# SCTransform (glmGamPoi), PCA elbow, clustree, and silhouette scoring.
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

.log_error <- function(...) {
    cat(sprintf("[%s] ERROR %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)), file = stderr())
}

.load_integration_packages <- function() {
    required <- c(
        "Seurat", "ggplot2", "patchwork", "clustree", "glmGamPoi",
        "cluster", "digest", "jsonlite", "yaml", "dplyr", "tidyr",
        "ggpubr", "rstatix", "dittoSeq"
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
        library(ggpubr)
        library(rstatix)
        library(dittoSeq)
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

.join_rna_layers_if_needed <- function(obj) {
    assay_name <- "RNA"
    if (!assay_name %in% Assays(obj)) return(obj)
    assay_obj <- tryCatch(obj[[assay_name]], error = function(e) NULL)
    if (!is.null(assay_obj) && inherits(assay_obj, "Assay5") && exists("JoinLayers", mode = "function")) {
        obj <- JoinLayers(obj, assay = assay_name)
    }
    obj
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
    fig_dir <- file.path(output_dir, "figures")
    tab_dir <- file.path(output_dir, "tables")
    log_dir <- file.path(output_dir, "logs")
    for (d in c(prep_dir, fig_dir, tab_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    log_file <- file.path(log_dir, "integration_prep.log")
    sink(log_file, split = TRUE)
    on.exit(sink(), add = TRUE)

    prepared_rds <- file.path(prep_dir, "sct_prepared.rds")
    meta_json <- file.path(prep_dir, ".sct_prepared.meta.json")

    n_features <- as.integer(cfg$n_features %||% 2500L)
    npcs <- as.integer(cfg$npcs %||% 50L)
    neighbor_dims <- .parse_dims(cfg$neighbor_dims %||% "1:30")
    mito_col <- cfg$mito_regress %||% "percent.mito"
    resolutions <- .parse_resolutions(cfg$resolutions)
    cluster_prefix <- cfg$cluster_prefix %||% "SCT_snn_res."
    sil_subsample <- as.integer(cfg$silhouette_subsample %||% 5000L)
    sil_seed <- as.integer(cfg$silhouette_seed %||% 42L)

    run_params <- list(
        input_rds = input_rds,
        annotated_h5ad = annotated_h5ad,
        n_features = n_features,
        npcs = npcs,
        neighbor_dims = neighbor_dims,
        mito_regress = mito_col,
        resolutions = resolutions,
        cluster_prefix = cluster_prefix,
        silhouette_subsample = sil_subsample
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
            .log_info("Annotated h5ad registered for later steps: %s", annotated_h5ad)
        } else {
            .log_info("Annotated h5ad path not found (continuing prep): %s", annotated_h5ad)
        }
    }

    obj <- .join_rna_layers_if_needed(obj)
    obj <- .ensure_mito_percent(obj, mito_col)

    .log_info("Normalizing RNA counts.")
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)

    .log_info("Running SCTransform (glmGamPoi, regress %s).", mito_col)
    obj <- SCTransform(
        obj,
        assay = "RNA",
        method = "glmGamPoi",
        vars.to.regress = mito_col,
        vst.flavor = "v2",
        verbose = FALSE
    )
    DefaultAssay(obj) <- "SCT"

    .log_info("Selecting variable features (n=%d).", n_features)
    obj <- FindVariableFeatures(obj, selection.method = "vst", assay = "SCT", nfeatures = n_features, verbose = FALSE)

    .log_info("Running PCA (npcs=%d).", npcs)
    obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)

    elbow_path <- file.path(fig_dir, "elbow_plot.pdf")
    .log_info("Saving elbow plot: %s", elbow_path)
    pdf(elbow_path, width = 8, height = 5)
    print(ElbowPlot(obj, ndims = npcs))
    dev.off()

    .log_info("Finding neighbors (dims %s).", paste(range(neighbor_dims), collapse = ":"))
    obj <- FindNeighbors(obj, dims = neighbor_dims, reduction = "pca", verbose = FALSE)

    .log_info("Testing resolutions: %s", paste(resolutions, collapse = ", "))
    obj <- FindClusters(obj, resolution = resolutions, verbose = FALSE)

    clustree_path <- file.path(fig_dir, "clustree_sct.pdf")
    .log_info("Saving clustree plot: %s", clustree_path)
    pdf(clustree_path, width = 12, height = 10)
    print(clustree(obj, prefix = cluster_prefix))
    dev.off()

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

    best_idx <- which.max(sil_df$silhouette)
    best_res <- if (length(best_idx)) sil_df$resolution[best_idx[1]] else NA_real_

    sil_plot_path <- file.path(fig_dir, "silhouette_sct.pdf")
    p_sil <- ggplot(sil_df, aes(x = resolution, y = silhouette)) +
        geom_point(size = 3) +
        geom_line() +
        theme_minimal() +
        labs(title = "Silhouette score vs resolution (subsampled)", x = "Resolution", y = "Silhouette")
    if (!is.na(best_res)) {
        p_sil <- p_sil +
            geom_vline(xintercept = best_res, color = "red", linetype = "dashed")
    }
    ggsave(sil_plot_path, plot = p_sil, width = 8, height = 5)
    .log_info("Saved silhouette plot: %s", sil_plot_path)

    summary_df <- data.frame(
        best_resolution = best_res,
        best_silhouette = if (length(best_idx)) sil_df$silhouette[best_idx[1]] else NA_real_,
        npcs = npcs,
        neighbor_dims = paste(range(neighbor_dims), collapse = ":"),
        n_features = n_features,
        stringsAsFactors = FALSE
    )
    summary_csv <- file.path(tab_dir, "resolution_summary.csv")
    write.csv(summary_df, summary_csv, row.names = FALSE)
    .log_info("Recommended resolution (silhouette): %s", best_res)

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
            silhouette_csv = sil_csv,
            resolution_summary = summary_csv
        ),
        best_resolution = best_res
    )
    jsonlite::write_json(meta_out, meta_json, auto_unbox = TRUE, pretty = TRUE)

    .log_info("Integration prep completed. Outputs under: %s", output_dir)
    quit(save = "no", status = 0)
}, error = function(e) {
    .log_error("%s", conditionMessage(e))
    quit(save = "no", status = 1)
})
