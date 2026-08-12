# =============================================================================
# HKOCA Integration - Methods stage (Harmony / RPCA / CCA + re-annotation)
#
# Reads the SCT-prepared Seurat RDS produced by integration_prep.R, runs one
# or more integration methods, generates per-method UMAP + FeaturePlot figures
# in dedicated subdirectories, optionally re-annotates with Snapseed (Level 3,
# stored as celltype_final), and saves one RDS per method.
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
# Marker helpers (same logic as integration_prep.R)
# =============================================================================

.categorical_levels <- function(x) {
    vals <- unique(trimws(as.character(x)))
    sort(vals[nzchar(vals) & !is.na(vals)])
}

.metadata_split_columns <- function(meta) {
    preferred <- c(
        "sample_id", "model", "sort_podxl", "study", "source",
        "diff_protocol", "sc_protocol", "sequencing", "genome_build",
        "Age", "type", "transduction", "MOI", "AAV"
    )
    is_cat <- vapply(meta, function(col) {
        if (!(is.character(col) || is.factor(col) || is.logical(col))) return(FALSE)
        n <- length(.categorical_levels(col))
        n >= 2L && n <= 24L
    }, logical(1))
    detected <- names(meta)[is_cat]
    c(intersect(preferred, detected), setdiff(detected, preferred))
}

.split_plot_layout <- function(n_panels) {
    if (n_panels <= 1L) return(list(num_columns = 1L, width = 8, height = 6))
    num_columns <- if (n_panels <= 2L) n_panels else min(4L, ceiling(sqrt(n_panels)))
    n_rows <- ceiling(n_panels / num_columns)
    list(num_columns = num_columns, width = max(8, 4.5 * num_columns), height = max(6, 3.8 * n_rows))
}

.feature_plot_layout <- function(n_genes) {
    nc <- min(3L, n_genes)
    list(ncol = nc, width = max(9, 3.5 * nc), height = max(3.5, 3.5 * ceiling(n_genes / nc)))
}

.extract_level3_marker_groups <- function(yaml_path) {
    if (is.na(yaml_path) || !nzchar(yaml_path) || !file.exists(yaml_path)) return(list())
    raw <- tryCatch(yaml::read_yaml(yaml_path), error = function(e) NULL)
    if (is.null(raw)) return(list())
    groups <- list()
    for (l1_name in names(raw)) {
        l1 <- raw[[l1_name]]; if (!is.list(l1)) next
        subs <- l1$subtypes; if (is.null(subs)) next
        for (l2_name in names(subs)) {
            l2 <- subs[[l2_name]]; if (!is.list(l2)) next
            l3_subs <- l2$subtypes
            if (!is.null(l3_subs)) {
                for (l3_name in names(l3_subs)) {
                    l3 <- l3_subs[[l3_name]]
                    genes <- if (is.character(l3$marker_genes)) l3$marker_genes else character(0)
                    if (length(genes) > 0) groups[[l3_name]] <- list(label = l3_name, genes = genes)
                }
            } else {
                genes <- if (is.character(l2$marker_genes)) l2$marker_genes else character(0)
                if (length(genes) > 0) groups[[l2_name]] <- list(label = l2_name, genes = genes)
            }
        }
    }
    groups
}

# =============================================================================
# Shared visualization: UMAP (cluster + split) + FeaturePlots per method
# =============================================================================

.save_method_figures <- function(obj, out_dir, reduction, cluster_col,
                                  markers_yaml_path, method_label) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    Idents(obj) <- cluster_col
    overall_path <- file.path(out_dir, "umap.pdf")
    overall_plot <- DimPlot(obj, reduction = reduction, label = TRUE, repel = TRUE) +
        ggtitle(sprintf("%s UMAP", method_label))
    ggsave(overall_path, plot = overall_plot, width = 8, height = 6)
    .log_info("[%s] Saved UMAP: %s", method_label, overall_path)

    split_cols <- .metadata_split_columns(obj@meta.data)
    for (col_name in split_cols) {
        levels_n <- length(.categorical_levels(obj@meta.data[[col_name]]))
        layout <- .split_plot_layout(levels_n)
        out_path <- file.path(out_dir, sprintf("umap_split_%s.pdf", col_name))
        split_plot <- DimPlot(obj, reduction = reduction, split.by = col_name,
                              ncol = layout$num_columns, label = FALSE) +
            ggtitle(sprintf("%s UMAP split by %s", method_label, col_name))
        ggsave(out_path, plot = split_plot, width = layout$width, height = layout$height)
        .log_info("[%s] Saved split UMAP: %s (column=%s, levels=%d)", method_label, out_path, col_name, levels_n)
    }

    # FeaturePlots use RNA assay for unscaled expression
    DefaultAssay(obj) <- "RNA"
    marker_groups <- .extract_level3_marker_groups(markers_yaml_path)
    for (grp in marker_groups) {
        present <- intersect(grp$genes, rownames(obj))
        if (length(present) == 0) {
            .log_info("[%s] FeaturePlot: skipping %s (no genes in object)", method_label, grp$label)
            next
        }
        safe_label <- gsub("[^A-Za-z0-9_-]", "_", grp$label)
        out_path <- file.path(out_dir, sprintf("featureplot_%s.png", safe_label))
        layout <- .feature_plot_layout(length(present))
        fp <- FeaturePlot(obj, features = present, reduction = reduction,
                          ncol = layout$ncol, pt.size = 0.5, order = TRUE) &
            theme(plot.title = element_text(size = 9))
        fp_titled <- fp + patchwork::plot_annotation(title = grp$label)
        ggsave(out_path, plot = fp_titled, width = layout$width, height = layout$height, dpi = 150)
        .log_info("[%s] Saved FeaturePlot: %s (%d genes)", method_label, out_path, length(present))
    }

    # If celltype_final exists, add re-annotation UMAP
    if ("celltype_final" %in% colnames(obj@meta.data)) {
        Idents(obj) <- "celltype_final"
        ann_path <- file.path(out_dir, "umap_celltype_final.pdf")
        ann_plot <- DimPlot(obj, reduction = reduction, label = TRUE, repel = TRUE,
                            group.by = "celltype_final") +
            ggtitle(sprintf("%s UMAP (celltype_final)", method_label))
        ggsave(ann_path, plot = ann_plot, width = 10, height = 7)
        .log_info("[%s] Saved re-annotated UMAP: %s", method_label, ann_path)
    }

    # Reset to SCT for downstream steps
    DefaultAssay(obj) <- "SCT"
    invisible(NULL)
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
    needs_norm <- TRUE
    if (exists("Layers", mode = "function")) {
        layers <- tryCatch(Layers(obj[["RNA"]]), error = function(e) character(0))
        needs_norm <- !("data" %in% layers)
    }
    if (needs_norm && exists("NormalizeData", mode = "function")) {
        obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
    }

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
# Per-method integration runners
# =============================================================================

.integrate_harmony <- function(obj, cfg, method_dir, markers_yaml, ann_dir) {
    dims <- .parse_dims(cfg$integration_dims %||% "1:30")
    res  <- as.numeric(cfg$harmony_resolution %||% "0.6")

    .log_info("[harmony] Running IntegrateLayers (HarmonyIntegration).")
    obj <- IntegrateLayers(
        object = obj, method = HarmonyIntegration,
        orig.reduction = "pca", new.reduction = "harmony",
        normalization.method = "SCT", verbose = TRUE
    )
    obj <- FindNeighbors(obj, reduction = "harmony", dims = dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = res, cluster.name = "harmony_clusters", verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "harmony", dims = dims,
                   reduction.name = "umap.harmony", verbose = FALSE)

    obj <- .run_snapseed_reannotation(obj, "harmony_clusters",
                                      markers_yaml, ann_dir, "harmony")
    .save_method_figures(obj, method_dir, "umap.harmony", "harmony_clusters",
                         markers_yaml, "Harmony")
    obj
}

.integrate_rpca <- function(obj, cfg, method_dir, markers_yaml, ann_dir) {
    dims <- .parse_dims(cfg$integration_dims %||% "1:30")
    res  <- as.numeric(cfg$rpca_resolution %||% "0.5")

    .log_info("[rpca] Running IntegrateLayers (RPCAIntegration).")
    obj <- IntegrateLayers(
        object = obj, method = RPCAIntegration,
        orig.reduction = "pca", new.reduction = "integrated.rpca",
        normalization.method = "SCT", verbose = TRUE
    )
    obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = res, cluster.name = "rpca_clusters", verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = dims,
                   reduction.name = "umap.rpca", verbose = FALSE)

    obj <- .run_snapseed_reannotation(obj, "rpca_clusters",
                                      markers_yaml, ann_dir, "rpca")
    .save_method_figures(obj, method_dir, "umap.rpca", "rpca_clusters",
                         markers_yaml, "RPCA")
    obj
}

.integrate_cca <- function(obj, cfg, method_dir, markers_yaml, ann_dir) {
    dims <- .parse_dims(cfg$integration_dims %||% "1:30")
    res  <- as.numeric(cfg$cca_resolution %||% "0.5")

    .log_info("[cca] Running IntegrateLayers (CCAIntegration).")
    obj <- IntegrateLayers(
        object = obj, method = CCAIntegration,
        orig.reduction = "pca", new.reduction = "integrated.cca",
        normalization.method = "SCT", verbose = TRUE
    )
    obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = res, cluster.name = "cca_clusters", verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "integrated.cca", dims = dims,
                   reduction.name = "umap.cca", verbose = FALSE)

    obj <- .run_snapseed_reannotation(obj, "cca_clusters",
                                      markers_yaml, ann_dir, "cca")
    .save_method_figures(obj, method_dir, "umap.cca", "cca_clusters",
                         markers_yaml, "CCA")
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

    fig_dir  <- file.path(output_dir, "figures")
    tab_dir  <- file.path(output_dir, "tables")
    obj_dir  <- file.path(output_dir, "objects")
    log_dir  <- file.path(output_dir, "logs")
    for (d in c(fig_dir, tab_dir, obj_dir, log_dir))
        dir.create(d, recursive = TRUE, showWarnings = FALSE)

    log_file <- file.path(log_dir, "integration_methods.log")
    sink(log_file, split = TRUE)
    on.exit(sink(), add = TRUE)

    .load_packages()

    .log_info("Loading prepared Seurat object: %s", prepared_rds)
    obj_base <- readRDS(prepared_rds)
    if (!inherits(obj_base, "Seurat"))
        stop("prepared_rds must be a Seurat object.")
    .log_info("Loaded %s cells x %s features", format(ncol(obj_base), big.mark = ","), format(nrow(obj_base), big.mark = ","))

    results <- list()

    for (method in methods) {
        .log_info("Starting integration method: %s", method)
        method_fig_dir <- file.path(fig_dir, method)
        method_ann_dir <- file.path(output_dir, "annotation", method)
        out_rds <- file.path(obj_dir, sprintf("integrated_%s.rds", method))

        if (!force_flag && file.exists(out_rds) && file.info(out_rds)$size > 0) {
            .log_info("[%s] Skipping, output already exists: %s", method, out_rds)
            results[[method]] <- list(rds = out_rds, status = "skipped")
            next
        }

        obj_method <- tryCatch({
            obj_work <- obj_base
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
        results[[method]] <- list(rds = out_rds, status = "done")
        .log_info("[%s] Completed.", method)
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
        prep_meta <- file.path(dirname(prepared_rds), ".sct_prepared.meta.json")
        .remove_path(prepared_rds, "SCT prepared RDS")
        .remove_path(prep_meta, "prep metadata")
    }

    .log_info("All integration methods completed. Outputs under: %s", output_dir)
    quit(save = "no", status = 0)

}, error = function(e) {
    .log_error("%s", conditionMessage(e))
    quit(save = "no", status = 1)
})
