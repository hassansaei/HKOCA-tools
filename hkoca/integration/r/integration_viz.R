# Shared integration figure helpers (UMAP, FeaturePlot, silhouette selection).

HKOCA_PLOT_DPI <- 300L
HKOCA_UMAP_PT_SIZE <- 1.5
HKOCA_UMAP_ALPHA <- 1.0
HKOCA_FEATURE_PT_SIZE <- 1.2
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
    if (grepl("_clusters$", col_name)) return(TRUE)
    if (col_name %in% c("celltype_final", "celltype", "seurat_clusters", "ident")) return(TRUE)
    FALSE
}

.metadata_split_columns <- function(meta) {
    preferred <- c(
        "sample_id", "model", "sort_podxl", "study", "source",
        "diff_protocol", "sc_protocol", "sequencing", "genome_build",
        "Age", "type", "transduction", "MOI", "AAV",
        "condition", "genotype", "reporter", "transgene_status"
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

.transgene_split_column <- function(meta) {
    for (col in c("sample_id", "orig.ident")) {
        if (!col %in% colnames(meta)) next
        if (length(.categorical_levels(meta[[col]])) >= 2L) return(col)
    }
    NA_character_
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

.resolution_label <- function(resolution) {
    if (is.null(resolution) || length(resolution) == 0 || is.na(resolution[1])) {
        return(NA_character_)
    }
    format(as.numeric(resolution[1]), trim = TRUE, scientific = FALSE)
}

.apply_point_style <- function(plot, size = HKOCA_UMAP_PT_SIZE, alpha = HKOCA_UMAP_ALPHA) {
    style_ggplot <- function(p) {
        if (!inherits(p, "ggplot")) return(p)
        for (i in seq_along(p$layers)) {
            geom_class <- class(p$layers[[i]]$geom)
            if (!any(grepl("Point|Scattermore", geom_class))) next
            if (any(grepl("Scattermore", geom_class))) {
                p$layers[[i]]$aes_params$pointsize <- max(size * 3, 6)
                p$layers[[i]]$aes_params$pixels <- 1200
            } else {
                p$layers[[i]]$aes_params$size <- size
            }
            p$layers[[i]]$aes_params$alpha <- alpha
            p$layers[[i]]$aes_params$stroke <- 0
        }
        p
    }
    if (inherits(plot, "patchwork")) {
        n <- tryCatch(length(plot), error = function(e) 0L)
        if (is.finite(n) && n > 0L) {
            for (i in seq_len(n)) {
                if (inherits(plot[[i]], "ggplot")) plot[[i]] <- style_ggplot(plot[[i]])
            }
        }
        return(plot)
    }
    style_ggplot(plot)
}

.hkoca_dim_plot <- function(obj, reduction, ..., group.by = NULL) {
    extra <- list(...)
    args <- c(
        list(
            object = obj,
            reduction = reduction,
            pt.size = HKOCA_UMAP_PT_SIZE,
            raster = FALSE
        ),
        extra
    )
    if (!is.null(group.by)) args$group.by <- group.by
    plot <- do.call(DimPlot, args)
    .apply_point_style(plot, size = HKOCA_UMAP_PT_SIZE, alpha = HKOCA_UMAP_ALPHA)
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

.hkoca_feature_plot <- function(obj, features, reduction, ncol, split.by = NULL) {
    args <- list(
        object = obj,
        features = features,
        reduction = reduction,
        ncol = ncol,
        pt.size = HKOCA_FEATURE_PT_SIZE,
        alpha = HKOCA_UMAP_ALPHA,
        order = TRUE,
        cols = HKOCA_FEATURE_COLS,
        raster = FALSE
    )
    if (!is.null(split.by) && nzchar(split.by)) args$split.by <- split.by
    fp <- do.call(FeaturePlot, args) & theme(plot.title = element_text(size = 9))
    .apply_point_style(fp, size = HKOCA_FEATURE_PT_SIZE, alpha = HKOCA_UMAP_ALPHA)
}

.parse_name_list <- function(x) {
    parts <- trimws(unlist(strsplit(as.character(x %||% ""), ",", fixed = TRUE)))
    parts[nzchar(parts) & !tolower(parts) %in% c("null", "none", "na")]
}

.resolve_transgenes <- function(cfg = NULL) {
    from_cfg <- .parse_name_list(if (!is.null(cfg)) cfg$transgenes else NULL)
    if (length(from_cfg)) return(unique(from_cfg))
    from_env <- .parse_name_list(Sys.getenv("HKOCA_TRANSGENES", unset = ""))
    if (length(from_env)) return(unique(from_env))
    c("AAV", "EGFP", "mCherry", "GFP", "eGFP", "LK03_eGFP")
}

.match_features_in_object <- function(obj, names) {
    wanted <- unique(.parse_name_list(names))
    if (!length(wanted)) return(character(0))
    feats <- rownames(obj)
    key <- stats::setNames(feats, toupper(feats))
    found <- vapply(wanted, function(g) {
        hit <- unname(key[toupper(g)])
        hit <- hit[!is.na(hit) & nzchar(hit)]
        if (length(hit)) hit[[1]] else NA_character_
    }, character(1))
    unique(found[!is.na(found)])
}

.assay_matrix <- function(obj, assay = DefaultAssay(obj), layers = c("counts", "data")) {
    for (layer in layers) {
        mat <- tryCatch(
            GetAssayData(obj, assay = assay, layer = layer),
            error = function(e) tryCatch(
                GetAssayData(obj, assay = assay, slot = layer),
                error = function(e2) NULL
            )
        )
        if (!is.null(mat) && nrow(mat) > 0) return(mat)
    }
    NULL
}

.prepare_rna_for_split_feature_plot <- function(obj, split_col) {
    if (!"RNA" %in% Assays(obj)) return(obj)
    obj <- .ensure_joined_normalized_rna(obj)
    if (!split_col %in% colnames(obj@meta.data)) return(obj)
    if (length(.categorical_levels(obj@meta.data[[split_col]])) < 2L) return(obj)
    DefaultAssay(obj) <- "RNA"
    if (exists("split", mode = "function")) {
        obj[["RNA"]] <- tryCatch(
            split(obj[["RNA"]], f = obj@meta.data[[split_col]]),
            error = function(e) {
                .log_warn("Could not split RNA by '%s': %s", split_col, conditionMessage(e))
                obj[["RNA"]]
            }
        )
    }
    obj
}

.features_with_expression <- function(obj, features) {
    features <- unique(features[features %in% rownames(obj)])
    if (!length(features)) return(character(0))
    mat <- .assay_matrix(obj)
    if (is.null(mat)) return(features)
    keep <- vapply(features, function(g) {
        if (!g %in% rownames(mat)) return(FALSE)
        as.numeric(Matrix::nnzero(mat[g, , drop = FALSE])) > 0
    }, logical(1))
    features[keep]
}

.save_one_feature_plot <- function(obj, features, reduction, out_dir, label, method_label) {
    if (!length(features)) return(NULL)
    safe_label <- gsub("[^A-Za-z0-9_-]", "_", label)
    out_path <- file.path(out_dir, sprintf("featureplot_%s.png", safe_label))
    layout <- .feature_plot_layout(length(features))
    fp <- tryCatch(
        .hkoca_feature_plot(obj, features, reduction, layout$ncol),
        error = function(e) {
            .log_warn(
                "[%s] FeaturePlot failed for %s: %s",
                method_label, label, conditionMessage(e)
            )
            NULL
        }
    )
    if (is.null(fp)) return(NULL)
    fp_titled <- fp + patchwork::plot_annotation(title = label)
    .save_ggplot_png(out_path, fp_titled, width = layout$width, height = layout$height)
    .log_info(
        "[%s] Saved FeaturePlot: %s (genes=%d, ncol=%d)",
        method_label, out_path, length(features), layout$ncol
    )
    list(path = out_path, genes = features, ncol = layout$ncol)
}

.save_split_feature_plots <- function(obj, features, reduction, out_dir, split_cols,
                                        label, method_label) {
    if (!length(features)) return(list())
    if (!length(split_cols)) {
        .log_info(
            "[%s] Split FeaturePlot skipped for %s: no metadata column with 2-24 levels.",
            method_label, label
        )
        return(list())
    }
    saved <- list()
    for (col_name in split_cols) {
        levels_n <- length(.categorical_levels(obj@meta.data[[col_name]]))
        if (levels_n < 2L) next
        layout <- .split_plot_layout(levels_n)
        obj_split <- .prepare_rna_for_split_feature_plot(obj, col_name)
        panels <- lapply(features, function(g) {
            tryCatch(
                .hkoca_feature_plot(
                    obj_split, g, reduction, layout$num_columns, split.by = col_name
                ),
                error = function(e) {
                    .log_warn(
                        "[%s] Split FeaturePlot failed for %s / %s: %s",
                        method_label, g, col_name, conditionMessage(e)
                    )
                    NULL
                }
            )
        })
        panels <- Filter(Negate(is.null), panels)
        if (!length(panels)) next
        combined <- if (length(panels) == 1L) {
            panels[[1]]
        } else {
            patchwork::wrap_plots(panels, ncol = 1L)
        }
        combined <- combined + patchwork::plot_annotation(
            title = sprintf("%s %s split by %s", method_label, label, col_name)
        )
        safe_label <- gsub("[^A-Za-z0-9_-]", "_", label)
        out_path <- file.path(out_dir, sprintf("featureplot_%s_split_%s.png", safe_label, col_name))
        height <- layout$height * max(1, length(panels))
        .save_ggplot_png(out_path, combined, width = layout$width, height = height)
        .log_info(
            "[%s] Saved split FeaturePlot: %s (column=%s, genes=%d, levels=%d)",
            method_label, out_path, col_name, length(panels), levels_n
        )
        saved[[col_name]] <- list(
            path = out_path, genes = features, n_levels = levels_n
        )
    }
    saved
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

.save_level_annotation_umaps <- function(obj, out_dir, reduction, palettes, method_label) {
    level_cols <- .annotation_level_columns(obj@meta.data)
    if (!length(level_cols)) {
        .log_warn(
            "[%s] No Level_1/2/3 columns; skipping annotation UMAPs (pass --annotated-h5ad at prep).",
            method_label
        )
        return(character(0))
    }
    saved <- character(0)
    fb <- if (!is.null(palettes)) palettes$fallback %||% "#BBBBBB" else "#BBBBBB"
    for (col_name in level_cols) {
        colors <- tryCatch(
            ditto_colors_for_meta(obj, col_name, palettes = palettes),
            error = function(e) {
                .log_warn("[%s] Color lookup failed for %s: %s", method_label, col_name, conditionMessage(e))
                NULL
            }
        )
        plot_args <- list(
            obj = obj,
            reduction = reduction,
            group.by = col_name,
            label = TRUE,
            repel = TRUE
        )
        # Named vector from celltype_colors.yaml; DimPlot matches identity names.
        if (!is.null(colors) && length(colors)) {
            colors[is.na(colors) | !nzchar(colors)] <- fb
            plot_args$cols <- colors
        }
        p <- do.call(.hkoca_dim_plot, plot_args) +
            ggtitle(sprintf("%s UMAP (%s)", method_label, gsub("_", " ", col_name, fixed = TRUE)))
        out_path <- file.path(out_dir, sprintf("umap_%s.png", col_name))
        .save_ggplot_png(out_path, p, width = 10, height = 7)
        .log_info(
            "[%s] Saved annotation UMAP: %s (palette=%s, labels=%d)",
            method_label, out_path, infer_palette_level(col_name), length(colors)
        )
        saved <- c(saved, out_path)
    }
    saved
}

.save_integration_figures <- function(obj, out_dir, reduction, cluster_col,
                                        markers_yaml_path, method_label,
                                        feature_assay = "RNA", palettes = NULL,
                                        resolution = NA_real_) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    res_lab <- .resolution_label(resolution)
    Idents(obj) <- cluster_col
    if (!is.na(res_lab)) {
        overall_path <- file.path(out_dir, sprintf("umap_res%s.png", res_lab))
        overall_title <- sprintf("%s UMAP (resolution %s)", method_label, res_lab)
    } else {
        overall_path <- file.path(out_dir, "umap.png")
        overall_title <- sprintf("%s UMAP", method_label)
    }
    overall_plot <- .apply_point_style(
        .hkoca_dim_plot(
            obj,
            reduction = reduction,
            label = TRUE,
            repel = TRUE
        ) + ggtitle(overall_title)
    )
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

    level_umaps <- .save_level_annotation_umaps(obj, out_dir, reduction, palettes, method_label)
    prop_paths <- .save_celltype_proportion_plots(
        obj, out_dir, group.by = "sample_id", palettes = palettes, method_label = method_label
    )

    marker_groups <- .extract_level3_marker_groups(markers_yaml_path)
    transgene_names <- .resolve_transgenes()
    feature_outputs <- list()
    need_features <- feature_assay %in% Assays(obj) &&
        (length(marker_groups) > 0 || length(transgene_names) > 0)
    if (need_features) {
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
            saved <- .save_one_feature_plot(
                obj, present, reduction, out_dir, grp$label, method_label
            )
            if (!is.null(saved)) feature_outputs[[grp$label]] <- saved
        }
        transgene_feats <- .match_features_in_object(obj, transgene_names)
        transgene_feats <- transgene_feats[transgene_feats %in% rownames(obj)]
        if (length(transgene_feats)) {
            saved <- .save_one_feature_plot(
                obj, transgene_feats, reduction, out_dir, "transgenes", method_label
            )
            if (!is.null(saved)) feature_outputs[["transgenes"]] <- saved
            split_col <- .transgene_split_column(obj@meta.data)
            if (is.na(split_col)) {
                .log_warn(
                    "[%s] Transgenes detected (%s) but no sample_id/orig.ident with >=2 levels for split FeaturePlot.",
                    method_label, paste(transgene_feats, collapse = ", ")
                )
            } else {
                split_fp <- .save_split_feature_plots(
                    obj, transgene_feats, reduction, out_dir, split_col,
                    "transgenes", method_label
                )
                if (length(split_fp)) {
                    feature_outputs[["transgenes_split"]] <- split_fp
                } else {
                    .log_warn(
                        "[%s] No split transgene FeaturePlots saved (split column: %s).",
                        method_label, split_col
                    )
                }
            }
        } else if (length(transgene_names)) {
            .log_info(
                "[%s] FeaturePlot: no transgene expression in object (looked for: %s)",
                method_label, paste(transgene_names, collapse = ", ")
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
        .log_info("[%s] Saved celltype_final UMAP: %s", method_label, ann_path)
    }

    if ("SCT" %in% Assays(obj)) {
        DefaultAssay(obj) <- "SCT"
    }

    list(
        overall = overall_path,
        split = split_outputs,
        feature_plots = feature_outputs,
        split_columns = split_cols,
        level_umaps = level_umaps,
        proportion_plots = prop_paths
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
