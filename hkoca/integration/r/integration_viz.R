# Shared integration figure helpers (UMAP, FeaturePlot, silhouette selection).

HKOCA_PLOT_DPI <- 300L
HKOCA_UMAP_PT_SIZE <- 1.2
HKOCA_UMAP_ALPHA <- 1.0
HKOCA_FEATURE_PT_SIZE <- 0.8
HKOCA_FEATURE_COLS <- c("grey80", "black")

.categorical_levels <- function(x) {
    vals <- unique(trimws(as.character(x)))
    sort(vals[nzchar(vals) & !is.na(vals)])
}

.is_excluded_split_column <- function(col_name) {
    col_name <- as.character(col_name)
    if (grepl("_res[0-9]", col_name)) return(TRUE)
    if (grepl("_res\\.", col_name)) return(TRUE)
    if (grepl("^SCT_snn_res", col_name)) return(TRUE)
    if (grepl("^RNA_snn_res", col_name)) return(TRUE)
    if (grepl("^integrated\\.", col_name)) return(TRUE)
    if (grepl("^leiden_res", col_name, ignore.case = TRUE)) return(TRUE)
    if (grepl("^Level_", col_name)) return(TRUE)
    if (col_name %in% c("celltype_final", "celltype")) return(TRUE)
    FALSE
}

.metadata_split_columns <- function(meta) {
    preferred <- c(
        "sample_id", "model", "sort_podxl", "study", "source",
        "diff_protocol", "sc_protocol", "sequencing", "genome_build",
        "Age", "type", "transduction", "MOI", "AAV"
    )
    is_categorical <- vapply(names(meta), function(col_name) {
        col <- meta[[col_name]]
        if (.is_excluded_split_column(col_name)) return(FALSE)
        if (!(is.character(col) || is.factor(col) || is.logical(col))) return(FALSE)
        n_levels <- length(.categorical_levels(col))
        n_levels >= 2L && n_levels <= 24L
    }, logical(1))
    detected <- names(meta)[is_categorical]
    c(intersect(preferred, detected), setdiff(detected, preferred))
}

.split_plot_layout <- function(n_panels) {
    if (n_panels <= 1L) {
        return(list(num_columns = 1L, width = 8, height = 6))
    }
    num_columns <- if (n_panels <= 2L) n_panels else min(4L, ceiling(sqrt(n_panels)))
    n_rows <- ceiling(n_panels / num_columns)
    width <- max(8, 4.5 * num_columns)
    height <- max(6, 3.8 * n_rows)
    list(num_columns = num_columns, width = width, height = height)
}

.feature_plot_layout <- function(n_genes) {
    ncol <- min(3L, n_genes)
    nrow <- ceiling(n_genes / ncol)
    list(ncol = ncol, width = max(9, 3.5 * ncol), height = max(3.5, 3.5 * nrow))
}

.save_ggplot_png <- function(path, plot, width, height, dpi = HKOCA_PLOT_DPI) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    ggsave(path, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
}

.annotation_level_columns <- function(meta) {
    preferred <- c("Level_1", "Level_2", "Level_3")
    present <- preferred[preferred %in% colnames(meta)]
    if (length(present)) return(present)
    # Fallback: resolution-suffixed aliases
    for (prefix in c("Level_1", "Level_2", "Level_3")) {
        hits <- grep(paste0("^", prefix, "(_|$)"), colnames(meta), value = TRUE)
        hits <- hits[!grepl("_score$", hits)]
        if (length(hits)) present <- c(present, hits[[1]])
    }
    unique(present)
}

.save_celltype_proportion_plots <- function(obj, out_dir, group.by = "sample_id",
                                              palettes = NULL, method_label = "Non-integrated") {
    if (!requireNamespace("dittoSeq", quietly = TRUE)) {
        .log_warn("[%s] Package dittoSeq not available; skipping proportion bar plots.", method_label)
        return(character(0))
    }
    if (is.null(palettes) || !is.list(palettes)) {
        .log_warn("[%s] No HKOCA palettes loaded; skipping proportion bar plots.", method_label)
        return(character(0))
    }

    level_cols <- .annotation_level_columns(obj@meta.data)
    if (!length(level_cols)) {
        .log_warn(
            "[%s] No Level_1/2/3 columns in metadata; skipping proportion plots (pass --annotated-h5ad).",
            method_label
        )
        return(character(0))
    }

    if (!group.by %in% colnames(obj@meta.data)) {
        alt <- if ("orig.ident" %in% colnames(obj@meta.data)) "orig.ident" else NA_character_
        if (is.na(alt)) {
            .log_warn("[%s] group.by '%s' missing; skipping proportion plots.", method_label, group.by)
            return(character(0))
        }
        .log_warn("[%s] '%s' missing; using '%s' for proportion plots.", method_label, group.by, alt)
        group.by <- alt
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saved <- character(0)
    for (var_col in level_cols) {
        colors <- tryCatch(
            ditto_colors_for_meta(obj, var_col, palettes = palettes),
            error = function(e) {
                .log_warn("[%s] Color lookup failed for %s: %s", method_label, var_col, conditionMessage(e))
                NULL
            }
        )
        if (is.null(colors) || !length(colors)) next

        title <- sprintf("%s proportions per %s", gsub("_", " ", var_col, fixed = TRUE), group.by)
        out_path <- file.path(out_dir, sprintf("celltype_proportion_%s.png", var_col))
        p <- tryCatch({
            dittoSeq::dittoBarPlot(
                obj,
                var = var_col,
                group.by = group.by,
                main = title
            ) +
                ggplot2::scale_fill_manual(
                    values = colors,
                    na.value = palettes$fallback %||% "#BBBBBB"
                ) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        }, error = function(e) {
            .log_warn("[%s] dittoBarPlot failed for %s: %s", method_label, var_col, conditionMessage(e))
            NULL
        })
        if (is.null(p)) next

        .save_ggplot_png(out_path, p, width = 8, height = 6)
        .log_info("[%s] Saved proportion plot: %s", method_label, out_path)
        saved <- c(saved, out_path)
    }
    saved
}

.hkoca_dim_plot <- function(obj, reduction, ..., group.by = NULL) {
    args <- list(
        object = obj,
        reduction = reduction,
        pt.size = HKOCA_UMAP_PT_SIZE,
        alpha = HKOCA_UMAP_ALPHA,
        raster = FALSE,
        ...
    )
    if (!is.null(group.by)) args$group.by <- group.by
    do.call(DimPlot, args)
}

.ensure_joined_normalized_rna <- function(obj) {
    if (!"RNA" %in% Assays(obj)) return(obj)
    assay_obj <- obj[["RNA"]]
    if (inherits(assay_obj, "Assay5") && exists("JoinLayers", mode = "function")) {
        obj[["RNA"]] <- JoinLayers(assay_obj)
    }
    has_data <- TRUE
    if (exists("Layers", mode = "function")) {
        ly <- tryCatch(Layers(obj[["RNA"]]), error = function(e) character(0))
        has_data <- any(ly == "data" | startsWith(ly, "data."))
    }
    if (!has_data) {
        obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
    }
    obj
}

.split_rna_by_batch <- function(obj) {
    if (!"RNA" %in% Assays(obj)) return(obj)
    batch_col <- if ("sample_id" %in% colnames(obj@meta.data)) {
        "sample_id"
    } else if ("orig.ident" %in% colnames(obj@meta.data)) {
        "orig.ident"
    } else {
        return(obj)
    }
    vals <- as.character(obj@meta.data[[batch_col]])
    n_levels <- length(unique(vals[nzchar(vals) & !is.na(vals)]))
    if (n_levels < 2L) return(obj)
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[batch_col]])
    obj
}

.hkoca_feature_plot <- function(obj, features, reduction, ncol) {
    FeaturePlot(
        obj,
        features = features,
        reduction = reduction,
        ncol = ncol,
        pt.size = HKOCA_FEATURE_PT_SIZE,
        alpha = HKOCA_UMAP_ALPHA,
        order = TRUE,
        cols = HKOCA_FEATURE_COLS,
        raster = FALSE
    ) & theme(plot.title = element_text(size = 9))
}

.extract_level3_marker_groups <- function(yaml_path) {
    if (is.na(yaml_path) || !nzchar(yaml_path) || !file.exists(yaml_path))
        return(list())
    raw <- tryCatch(yaml::read_yaml(yaml_path), error = function(e) NULL)
    if (is.null(raw)) return(list())
    groups <- list()
    for (l1_name in names(raw)) {
        l1 <- raw[[l1_name]]
        if (!is.list(l1)) next
        subs <- l1$subtypes
        if (is.null(subs)) next
        for (l2_name in names(subs)) {
            l2 <- subs[[l2_name]]
            if (!is.list(l2)) next
            l3_subs <- l2$subtypes
            if (!is.null(l3_subs)) {
                for (l3_name in names(l3_subs)) {
                    l3 <- l3_subs[[l3_name]]
                    genes <- if (is.character(l3$marker_genes)) l3$marker_genes else character(0)
                    if (length(genes) > 0)
                        groups[[l3_name]] <- list(label = l3_name, genes = genes)
                }
            } else {
                genes <- if (is.character(l2$marker_genes)) l2$marker_genes else character(0)
                if (length(genes) > 0)
                    groups[[l2_name]] <- list(label = l2_name, genes = genes)
            }
        }
    }
    groups
}

.save_integration_figures <- function(obj, out_dir, reduction, cluster_col,
                                        markers_yaml_path, method_label,
                                        feature_assay = "RNA") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    Idents(obj) <- cluster_col
    overall_path <- file.path(out_dir, "umap.png")
    overall_plot <- .hkoca_dim_plot(
        obj,
        reduction = reduction,
        label = TRUE,
        repel = TRUE
    ) + ggtitle(sprintf("%s UMAP", method_label))
    .save_ggplot_png(overall_path, overall_plot, width = 8, height = 6)
    .log_info("[%s] Saved UMAP: %s", method_label, overall_path)

    split_cols <- .metadata_split_columns(obj@meta.data)
    split_outputs <- list()
    for (col_name in split_cols) {
        levels_n <- length(.categorical_levels(obj@meta.data[[col_name]]))
        layout <- .split_plot_layout(levels_n)
        out_path <- file.path(out_dir, sprintf("umap_split_%s.png", col_name))
        split_plot <- .hkoca_dim_plot(
            obj,
            reduction = reduction,
            split.by = col_name,
            ncol = layout$num_columns,
            label = FALSE
        ) + ggtitle(sprintf("%s UMAP split by %s", method_label, col_name))
        .save_ggplot_png(out_path, split_plot, width = layout$width, height = layout$height)
        .log_info(
            "[%s] Saved split UMAP: %s (column=%s, levels=%d)",
            method_label, out_path, col_name, levels_n
        )
        split_outputs[[col_name]] <- list(
            path = out_path, n_levels = levels_n, ncol = layout$num_columns,
            width = layout$width, height = layout$height
        )
    }

    marker_groups <- .extract_level3_marker_groups(markers_yaml_path)
    feature_outputs <- list()
    if (length(marker_groups) > 0 && feature_assay %in% Assays(obj)) {
        if (identical(feature_assay, "RNA")) {
            obj <- .ensure_joined_normalized_rna(obj)
        }
        DefaultAssay(obj) <- feature_assay
        for (grp in marker_groups) {
            present <- intersect(grp$genes, rownames(obj))
            if (length(present) == 0) {
                .log_info("[%s] FeaturePlot: skipping %s (no genes found in object)", method_label, grp$label)
                next
            }
            safe_label <- gsub("[^A-Za-z0-9_-]", "_", grp$label)
            out_path <- file.path(out_dir, sprintf("featureplot_%s.png", safe_label))
            layout <- .feature_plot_layout(length(present))
            fp <- tryCatch(
                .hkoca_feature_plot(obj, present, reduction, layout$ncol),
                error = function(e) {
                    .log_warn(
                        "[%s] FeaturePlot failed for %s: %s",
                        method_label, grp$label, conditionMessage(e)
                    )
                    NULL
                }
            )
            if (is.null(fp)) next
            fp_titled <- fp + patchwork::plot_annotation(title = grp$label)
            .save_ggplot_png(out_path, fp_titled, width = layout$width, height = layout$height)
            .log_info(
                "[%s] Saved FeaturePlot: %s (genes=%d, ncol=%d)",
                method_label, out_path, length(present), layout$ncol
            )
            feature_outputs[[grp$label]] <- list(
                path = out_path, genes = present, ncol = layout$ncol
            )
        }
        if (identical(feature_assay, "RNA")) {
            obj <- .split_rna_by_batch(obj)
        }
    }

    if ("celltype_final" %in% colnames(obj@meta.data)) {
        Idents(obj) <- "celltype_final"
        ann_path <- file.path(out_dir, "umap_celltype_final.png")
        ann_plot <- .hkoca_dim_plot(
            obj,
            reduction = reduction,
            label = TRUE,
            repel = TRUE,
            group.by = "celltype_final"
        ) + ggtitle(sprintf("%s UMAP (celltype_final)", method_label))
        .save_ggplot_png(ann_path, ann_plot, width = 10, height = 7)
        .log_info("[%s] Saved re-annotated UMAP: %s", method_label, ann_path)
    }

    if ("SCT" %in% Assays(obj)) {
        DefaultAssay(obj) <- "SCT"
    }

    list(
        overall = overall_path,
        split = split_outputs,
        feature_plots = feature_outputs,
        split_columns = split_cols
    )
}

.select_best_silhouette_resolution <- function(sil_df, obj, cluster_prefix,
                                               min_resolution = 0.4,
                                               min_clusters = 5L) {
    eligible <- sil_df[is.finite(sil_df$silhouette), , drop = FALSE]
    if (!nrow(eligible)) {
        return(list(resolution = NA_real_, silhouette = NA_real_, eligible = eligible))
    }

    filtered <- eligible[eligible$resolution >= min_resolution, , drop = FALSE]
    if (!nrow(filtered)) {
        .log_warn(
            "No silhouette scores at resolution >= %s; using full resolution range.",
            min_resolution
        )
        filtered <- eligible
    }

    if (!is.null(obj) && nrow(filtered) > 0) {
        n_clusters <- vapply(filtered$resolution, function(r) {
            col <- paste0(cluster_prefix, r)
            if (!col %in% colnames(obj@meta.data)) return(0L)
            length(unique(as.character(obj@meta.data[[col]])))
        }, integer(1))
        with_min_clusters <- filtered[n_clusters >= min_clusters, , drop = FALSE]
        if (nrow(with_min_clusters) > 0) {
            filtered <- with_min_clusters
        } else {
            .log_warn(
                "No resolutions with >= %d clusters; ignoring cluster-count filter.",
                min_clusters
            )
        }
    }

    best_idx <- which.max(filtered$silhouette)
    list(
        resolution = filtered$resolution[best_idx[1]],
        silhouette = filtered$silhouette[best_idx[1]],
        eligible = filtered
    )
}
