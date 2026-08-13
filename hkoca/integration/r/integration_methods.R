# =============================================================================
# HKOCA Integration - Methods stage (Harmony / RPCA / CCA + re-annotation)
#
# Reads the SCT-prepared Seurat RDS from integration_prep.R (per-sample SCT +
# merge + PCA). Runs IntegrateLayers only — no second SCTransform.
#
# Usage:
#   Rscript integration_methods.R \
#       --prepared_rds  PATH/prep/sct_prepared.rds \
#       --output_dir    PATH \
#       [--methods      harmony,rpca,cca] \
#       [--annotated_h5ad PATH] \
#       [--config       PATH/integration.config.dcf]
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

.parse_dims <- function(x) {
    txt <- trimws(as.character(x))
    if (grepl(":", txt)) {
        bounds <- as.integer(unlist(strsplit(txt, ":", fixed = TRUE)))
        return(seq.int(bounds[1], bounds[2]))
    }
    as.integer(unlist(strsplit(txt, ",", fixed = TRUE)))
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
    cfg <- list(); current_key <- NULL; current_val <- ""
    flush <- function() { if (!is.null(current_key)) cfg[[current_key]] <<- current_val }
    for (ln in lines) {
        if (grepl("^\\s*#", ln) || !nzchar(trimws(ln))) next
        key_parts <- regmatches(ln, regexec("^\\s*([^:]+):\\s*(.*)$", ln))[[1]]
        if (length(key_parts) == 3) {
            flush(); current_key <- trimws(key_parts[2]); current_val <- key_parts[3]
        } else if (!is.null(current_key)) {
            current_val <- paste(current_val, sub("\\s+$", "", ln), sep = "\n")
        }
    }
    flush(); cfg
}

.log_info <- function(...) cat(sprintf("[%s] INFO  %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
.log_warn <- function(...) cat(sprintf("[%s] WARN  %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
.log_error <- function(...) cat(sprintf("[%s] ERROR %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)), file = stderr())

.load_packages <- function() {
    required <- c(
        "Seurat", "ggplot2", "patchwork", "harmony",
        "cluster", "digest", "jsonlite", "yaml", "dplyr", "reticulate"
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
        library(harmony)
        library(cluster)
        library(digest)
        library(jsonlite)
        library(yaml)
        library(dplyr)
        library(reticulate)
    })
    script_dir <- .get_script_dir()
    viz_r <- file.path(script_dir, "integration_viz.R")
    if (!file.exists(viz_r))
        stop(sprintf("Missing integration visualization helpers: %s", viz_r))
    source(viz_r, local = FALSE)
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

.propagate_parent_labels_df <- function(assignments) {
    if (is.null(assignments) || !ncol(assignments)) return(assignments)
    label_cols <- grep("_score$|_auc$|_expr$", colnames(assignments), invert = TRUE, value = TRUE)
    if (length(label_cols) < 2L) return(assignments)
    out <- assignments
    for (i in seq_len(length(label_cols) - 1L)) {
        parent <- label_cols[i]
        child <- label_cols[i + 1L]
        parent_vals <- as.character(out[[parent]])
        child_vals <- as.character(out[[child]])
        missing <- !nzchar(child_vals) | is.na(child_vals) |
            tolower(child_vals) %in% c("unassigned", "na", "nan")
        out[[child]][missing] <- parent_vals[missing]
    }
    out
}

# =============================================================================
# Snapseed re-annotation (Level 3 -> celltype_final)
# =============================================================================

.run_snapseed_reannotation <- function(obj, cluster_col, markers_yaml_path,
                                        ann_dir, method_label) {
    if (!nzchar(markers_yaml_path) || !file.exists(markers_yaml_path)) {
        .log_warn("[%s] Markers YAML not found, skipping re-annotation.", method_label)
        return(obj)
    }
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        .log_warn("[%s] Package 'reticulate' not available, skipping re-annotation.", method_label)
        return(obj)
    }

    py <- Sys.getenv("HKOCA_ANNOTATION_PYTHON", unset = "")
    if (!nzchar(py)) py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
    if (nzchar(py)) {
        tryCatch(reticulate::use_python(py, required = TRUE), error = function(e) {
            .log_warn("[%s] Could not use annotation Python (%s): %s", method_label, py, conditionMessage(e))
        })
    }
    if (!reticulate::py_module_available("snapseed")) {
        .log_warn(
            "[%s] Python snapseed not available (set HKOCA_ANNOTATION_PYTHON to hkoca_harmonize python). Skipping re-annotation.",
            method_label
        )
        return(obj)
    }
    if (!reticulate::py_module_available("scanpy")) {
        .log_warn("[%s] Python scanpy not available, skipping re-annotation.", method_label)
        return(obj)
    }

    dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)
    .log_info("[%s] Running Snapseed re-annotation (Level 3 -> celltype_final).", method_label)

    marker_dict <- tryCatch(yaml::read_yaml(markers_yaml_path), error = function(e) NULL)
    if (is.null(marker_dict)) {
        .log_warn("[%s] Could not parse markers YAML, skipping re-annotation.", method_label)
        return(obj)
    }

    Idents(obj) <- cluster_col
    DefaultAssay(obj) <- "RNA"
    obj <- .ensure_joined_normalized_rna(obj)

    data_mat <- tryCatch(
        GetAssayData(obj, assay = "RNA", layer = "data"),
        error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
    )
    snapseed <- reticulate::import("snapseed")
    sc <- reticulate::import("scanpy")

    X <- reticulate::r_to_py(Matrix::t(data_mat))
    obs <- reticulate::r_to_py(obj@meta.data, convert = TRUE)
    var <- reticulate::r_to_py(data.frame(index = rownames(obj)), convert = TRUE)
    adata <- sc$AnnData(X = X, obs = obs, var = var)
    adata$obs[[cluster_col]] <- reticulate::r_to_py(as.character(Idents(obj)))

    marker_py <- reticulate::r_to_py(marker_dict, convert = TRUE)
    results <- tryCatch(
        snapseed$annotate_hierarchy(adata, marker_py, group_name = cluster_col),
        error = function(e) {
            .log_warn("[%s] Snapseed failed: %s", method_label, conditionMessage(e))
            NULL
        }
    )
    if (is.null(results)) return(obj)

    assignments <- reticulate::py_to_r(results[["assignments"]])
    metrics <- tryCatch(reticulate::py_to_r(results[["metrics"]]), error = function(e) NULL)
    if (is.null(assignments) || !nrow(assignments)) {
        .log_warn("[%s] Snapseed returned no assignments.", method_label)
        return(obj)
    }

    if (!"cluster" %in% colnames(assignments) && !is.null(rownames(assignments))) {
        assignments$cluster <- rownames(assignments)
    }
    assignments <- .propagate_parent_labels_df(assignments)

    label_cols <- grep("_score$|_auc$|_expr$", colnames(assignments), invert = TRUE, value = TRUE)
    level3_col <- tail(label_cols, 1)
    cluster_ids <- if ("cluster" %in% colnames(assignments)) {
        as.character(assignments$cluster)
    } else {
        rownames(assignments)
    }
    cluster_to_label <- setNames(as.character(assignments[[level3_col]]), cluster_ids)
    obj$celltype_final <- unname(cluster_to_label[as.character(Idents(obj))])

    score_col <- paste0(level3_col, "_score")
    if (score_col %in% colnames(assignments)) {
        score_map <- setNames(as.numeric(assignments[[score_col]]), cluster_ids)
        obj$celltype_final_score <- unname(score_map[as.character(Idents(obj))])
    }

    ann_csv <- file.path(ann_dir, "celltype_final_per_cluster.csv")
    write.csv(assignments, ann_csv, row.names = FALSE)
    .log_info("[%s] Saved re-annotation table: %s", method_label, ann_csv)

    if (!is.null(metrics)) {
        metrics_dir <- file.path(ann_dir, "snapseed_metrics")
        dir.create(metrics_dir, recursive = TRUE, showWarnings = FALSE)
        if (is.list(metrics)) {
            for (nm in names(metrics)) {
                metric_df <- tryCatch(reticulate::py_to_r(metrics[[nm]]), error = function(e) NULL)
                if (!is.null(metric_df) && nrow(metric_df) > 0) {
                    write.csv(metric_df, file.path(metrics_dir, paste0(nm, ".csv")), row.names = FALSE)
                }
            }
        }
    }

    per_cell <- data.frame(
        cell = colnames(obj),
        cluster = as.character(Idents(obj)),
        celltype_final = obj$celltype_final,
        celltype_final_score = if ("celltype_final_score" %in% colnames(obj@meta.data)) obj$celltype_final_score else NA_real_,
        stringsAsFactors = FALSE
    )
    write.csv(per_cell, file.path(ann_dir, "celltype_final_per_cell.csv"), row.names = FALSE)
    .log_info("[%s] Saved per-cell re-annotation: %s", method_label, file.path(ann_dir, "celltype_final_per_cell.csv"))

    obj
}

# =============================================================================
# Validate prep object for IntegrateLayers (no second SCTransform)
# =============================================================================

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

.count_assay5_batch_layers <- function(obj, assay = "RNA") {
    if (!assay %in% Assays(obj)) return(0L)
    if (!exists("Layers", mode = "function")) return(1L)
    layers <- tryCatch(Layers(obj[[assay]]), error = function(e) character(0))
    if (!length(layers)) return(0L)
    data_layers <- grep("^data\\.", layers, value = TRUE)
    if (!length(data_layers)) data_layers <- grep("^counts\\.", layers, value = TRUE)
    length(data_layers)
}

.count_sct_models <- function(obj) {
    if (!"SCT" %in% Assays(obj)) return(0L)
    assay_obj <- obj[["SCT"]]
    if (inherits(assay_obj, "SCTAssay") && "SCTModel.list" %in% slotNames(assay_obj)) {
        return(length(slot(assay_obj, "SCTModel.list")))
    }
    .count_assay5_batch_layers(obj, "SCT")
}

.validate_prepared_for_integration <- function(obj, cfg) {
    batch_col <- .detect_batch_column(obj, cfg$integration_batch_key %||% "sample_id")
    n_batches <- .count_batch_levels(obj, batch_col)
    n_sct_models <- .count_sct_models(obj)
    n_rna_layers <- .count_assay5_batch_layers(obj, "RNA")

    .log_info(
        "Prepared object check: batch='%s' (%d), SCT models=%d, RNA layers=%d.",
        batch_col, n_batches, n_sct_models, n_rna_layers
    )

    if (!"SCT" %in% Assays(obj)) {
        stop("Prepared object has no SCT assay. Re-run: hkoca integration prep")
    }
    if (!"pca" %in% Reductions(obj)) {
        stop("Prepared object has no PCA reduction. Re-run: hkoca integration prep")
    }

    if (n_batches < 2L) {
        return(list(
            object = obj,
            batch_col = batch_col,
            n_batches = n_batches,
            integratable = FALSE,
            skip_reason = "single_batch"
        ))
    }

    if (n_sct_models < 2L && n_rna_layers < 2L) {
        return(list(
            object = obj,
            batch_col = batch_col,
            n_batches = n_batches,
            integratable = FALSE,
            skip_reason = sprintf(
                "need_split_sct_prep (SCT models=%d, RNA layers=%d); re-run integration prep",
                n_sct_models, n_rna_layers
            )
        ))
    }

    DefaultAssay(obj) <- "SCT"
    list(
        object = obj,
        batch_col = batch_col,
        n_batches = n_batches,
        integratable = TRUE,
        skip_reason = NA_character_
    )
}

# =============================================================================
# Per-method integration runners
# =============================================================================

.integrate_harmony <- function(obj, cfg, method_dir, markers_yaml, ann_dir) {
    dims <- .parse_dims(cfg$integration_dims %||% "1:30")
    res  <- as.numeric(cfg$harmony_resolution %||% "0.6")

    .log_info("[harmony] Running IntegrateLayers (HarmonyIntegration).")
    obj <- IntegrateLayers(
        object = obj, method = HarmonyIntegration,
        orig.reduction = "pca", new.reduction = "harmony",
        assay = "SCT", normalization.method = "SCT", verbose = TRUE
    )
    obj <- FindNeighbors(obj, reduction = "harmony", dims = dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = res, cluster.name = "harmony_clusters", verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "harmony", dims = dims,
                   reduction.name = "umap.harmony", verbose = FALSE)

    obj <- .run_snapseed_reannotation(obj, "harmony_clusters",
                                      markers_yaml, ann_dir, "harmony")
    .save_integration_figures(obj, method_dir, "umap.harmony", "harmony_clusters",
                         markers_yaml, "Harmony", feature_assay = "RNA")
    obj
}

.integrate_rpca <- function(obj, cfg, method_dir, markers_yaml, ann_dir) {
    dims <- .parse_dims(cfg$integration_dims %||% "1:30")
    res  <- as.numeric(cfg$rpca_resolution %||% "0.5")

    .log_info("[rpca] Running IntegrateLayers (RPCAIntegration).")
    obj <- IntegrateLayers(
        object = obj, method = RPCAIntegration,
        orig.reduction = "pca", new.reduction = "integrated.rpca",
        assay = "SCT", normalization.method = "SCT", verbose = TRUE
    )
    obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = res, cluster.name = "rpca_clusters", verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = dims,
                   reduction.name = "umap.rpca", verbose = FALSE)

    obj <- .run_snapseed_reannotation(obj, "rpca_clusters",
                                      markers_yaml, ann_dir, "rpca")
    .save_integration_figures(obj, method_dir, "umap.rpca", "rpca_clusters",
                         markers_yaml, "RPCA", feature_assay = "RNA")
    obj
}

.integrate_cca <- function(obj, cfg, method_dir, markers_yaml, ann_dir) {
    dims <- .parse_dims(cfg$integration_dims %||% "1:30")
    res  <- as.numeric(cfg$cca_resolution %||% "0.5")

    .log_info("[cca] Running IntegrateLayers (CCAIntegration).")
    obj <- IntegrateLayers(
        object = obj, method = CCAIntegration,
        orig.reduction = "pca", new.reduction = "integrated.cca",
        assay = "SCT", normalization.method = "SCT", verbose = TRUE
    )
    obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = res, cluster.name = "cca_clusters", verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "integrated.cca", dims = dims,
                   reduction.name = "umap.cca", verbose = FALSE)

    obj <- .run_snapseed_reannotation(obj, "cca_clusters",
                                      markers_yaml, ann_dir, "cca")
    .save_integration_figures(obj, method_dir, "umap.cca", "cca_clusters",
                         markers_yaml, "CCA", feature_assay = "RNA")
    obj
}

# =============================================================================
# Main
# =============================================================================

tryCatch({
    script_dir   <- .get_script_dir()
    default_cfg  <- file.path(script_dir, "integration.config.dcf")
    cli_args     <- .parse_cli_args(commandArgs(trailingOnly = TRUE))

    if (.as_bool(cli_args$help, FALSE)) {
        cat("Usage:\n  Rscript integration_methods.R --prepared_rds PATH --output_dir PATH [--methods harmony,rpca,cca] [OPTIONS]\n")
        quit(save = "no", status = 0)
    }

    config_path <- if (!is.null(cli_args$config)) .resolve_path(cli_args$config) else default_cfg
    cfg <- .read_config_dcf(config_path)
    for (nm in setdiff(names(cli_args), c("config", "help")))
        cfg[[nm]] <- cli_args[[nm]]

    r_max_vsize <- cfg$r_max_vsize %||% "150Gb"
    Sys.setenv(R_MAX_VSIZE = r_max_vsize)
    .log_info("Set R_MAX_VSIZE=%s", r_max_vsize)

    prepared_rds  <- .resolve_path(cfg$prepared_rds)
    output_dir    <- .resolve_path(cfg$output_dir %||% "integration_results")
    force_flag    <- .as_bool(cfg$force_overwrite, FALSE)
    remove_flag   <- .as_bool(cfg$remove_intermediate, FALSE)
    methods_raw   <- tolower(trimws(cfg$methods %||% "harmony,rpca,cca"))
    methods       <- trimws(unlist(strsplit(methods_raw, ",", fixed = TRUE)))
    methods       <- methods[nzchar(methods)]

    if (is.na(prepared_rds) || !nzchar(prepared_rds) || !file.exists(prepared_rds))
        stop("Provide --prepared_rds pointing to the SCT-prepared RDS from integration_prep.")

    # Resolve markers YAML
    markers_yaml <- .resolve_path(cfg$markers_yaml)
    if (is.na(markers_yaml) || !nzchar(markers_yaml))
        markers_yaml <- Sys.getenv("HKOCA_MARKERS_YAML", unset = "")
    if (!nzchar(markers_yaml))
        markers_yaml <- normalizePath(
            file.path(script_dir, "..", "..", "config", "snapseed_markers_v4.yaml"),
            mustWork = FALSE
        )
    if (!file.exists(markers_yaml)) {
        .log_warn("Markers YAML not found, FeaturePlots and re-annotation will be skipped: %s", markers_yaml)
        markers_yaml <- ""
    } else {
        .log_info("Using markers YAML: %s", markers_yaml)
    }

    tab_dir  <- file.path(output_dir, "tables")
    obj_dir  <- file.path(output_dir, "objects")
    for (d in c(tab_dir, obj_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    log_file <- file.path(output_dir, "integration_methods.log")
    sink(log_file, split = TRUE)
    on.exit(sink(), add = TRUE)

    .load_packages()

    .log_info("Loading prepared Seurat object: %s", prepared_rds)
    obj_check <- readRDS(prepared_rds)
    if (!inherits(obj_check, "Seurat"))
        stop("prepared_rds must be a Seurat object.")
    .log_info("Loaded %s cells x %s features", format(ncol(obj_check), big.mark = ","), format(nrow(obj_check), big.mark = ","))

    prep <- .validate_prepared_for_integration(obj_check, cfg)
    rm(obj_check)
    gc(verbose = FALSE)

    if (!prep$integratable) {
        reason <- prep$skip_reason %||% "unknown"
        if (identical(reason, "single_batch")) {
            .log_warn(
                "Harmony/RPCA/CCA require >= 2 batches in '%s' (found %d). Skipping integration run.",
                prep$batch_col, prep$n_batches
            )
            .log_warn(
                "For single-sample data use prep/nonintegrated outputs. To integrate, provide an RDS with multiple samples (sample_id levels)."
            )
            summary_df <- data.frame(
                method = methods,
                rds = NA_character_,
                status = "skipped_single_batch",
                stringsAsFactors = FALSE
            )
            summary_csv <- file.path(tab_dir, "integration_methods_summary.csv")
            write.csv(summary_df, summary_csv, row.names = FALSE)
            .log_info("Integration methods summary: %s", summary_csv)
            quit(save = "no", status = 0)
        }
        stop(sprintf(
            "Prepared object not ready for IntegrateLayers (%s). Re-run: hkoca integration prep --force-overwrite",
            reason
        ))
    }

    .log_info(
        "Using prep SCT+PCA object for IntegrateLayers (batch='%s', %d levels).",
        prep$batch_col, prep$n_batches
    )

    results <- list()

    for (method in methods) {
        .log_info("Starting integration method: %s", method)
        method_fig_dir <- file.path(output_dir, method)
        method_ann_dir <- file.path(output_dir, "annotation", method)
        out_rds <- file.path(obj_dir, sprintf("integrated_%s.rds", method))

        if (!force_flag && file.exists(out_rds) && file.info(out_rds)$size > 0) {
            .log_info("[%s] Skipping, output already exists: %s", method, out_rds)
            results[[method]] <- list(rds = out_rds, status = "skipped")
            next
        }

        obj_method <- tryCatch({
            obj_work <- readRDS(prepared_rds)
            DefaultAssay(obj_work) <- "SCT"
            switch(method,
                harmony = .integrate_harmony(obj_work, cfg, method_fig_dir, markers_yaml, method_ann_dir),
                rpca    = .integrate_rpca(obj_work, cfg, method_fig_dir, markers_yaml, method_ann_dir),
                cca     = .integrate_cca(obj_work, cfg, method_fig_dir, markers_yaml, method_ann_dir),
                {
                    .log_warn("Unknown method '%s', skipping.", method)
                    NULL
                }
            )
        }, error = function(e) {
            .log_error("[%s] Integration failed: %s", method, conditionMessage(e))
            NULL
        })

        if (is.null(obj_method)) {
            results[[method]] <- list(rds = NA_character_, status = "failed")
            next
        }

        .log_info("[%s] Saving integrated object: %s", method, out_rds)
        saveRDS(obj_method, out_rds)
        bench_csv <- file.path(output_dir, "benchmark", sprintf("%s_embeddings.csv", method))
        dir.create(dirname(bench_csv), recursive = TRUE, showWarnings = FALSE)
        red <- switch(method, harmony = "harmony", rpca = "integrated.rpca", cca = "integrated.cca", NA_character_)
        if (!is.na(red) && red %in% Reductions(obj_method)) {
            emb <- as.matrix(Embeddings(obj_method, red))
            npcs_keep <- min(30L, ncol(emb))
            write.csv(emb[, seq_len(npcs_keep), drop = FALSE], bench_csv, row.names = TRUE)
            .log_info("[%s] Wrote embeddings: %s", method, bench_csv)
        }
        results[[method]] <- list(rds = out_rds, status = "done")
        .log_info("[%s] Completed.", method)
        rm(obj_method)
        gc(verbose = FALSE)
    }

    summary_df <- data.frame(
        method = names(results),
        rds    = vapply(results, `[[`, character(1), "rds"),
        status = vapply(results, `[[`, character(1), "status"),
        stringsAsFactors = FALSE
    )
    summary_csv <- file.path(tab_dir, "integration_methods_summary.csv")
    write.csv(summary_df, summary_csv, row.names = FALSE)
    .log_info("Integration methods summary: %s", summary_csv)

    failed <- sum(summary_df$status == "failed")
    if (failed > 0) {
        .log_warn("%d method(s) failed. Check logs.", failed)
        quit(save = "no", status = 1)
    }

    if (remove_flag) {
        prep_dir <- dirname(prepared_rds)
        prep_meta <- file.path(prep_dir, ".sct_prepared.meta.json")
        split_cache <- file.path(prep_dir, "integration_split_prepared.rds")
        .remove_path(prepared_rds, "SCT prepared RDS")
        .remove_path(split_cache, "legacy split integration RDS")
        .remove_path(prep_meta, "prep metadata")
    }

    .log_info("All integration methods completed. Outputs under: %s", output_dir)
    quit(save = "no", status = 0)

}, error = function(e) {
    .log_error("%s", conditionMessage(e))
    quit(save = "no", status = 1)
})
