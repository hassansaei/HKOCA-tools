# Shared HKOCA cell-type color palettes (same YAML as hkoca annotation).
#
# Usage:
#   palettes <- load_hkoca_celltype_palettes()
#   cols <- colors_for_labels(c("Podocytes", "PT"), level = "level3", palettes = palettes)
#   p + scale_color_hkoca_celltypes(level = "level3")

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

.load_yaml_safe <- function(path) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
        if (!requireNamespace("jsonlite", quietly = TRUE))
            stop("Install yaml or jsonlite to read HKOCA color palettes.")
    }
    if (requireNamespace("yaml", quietly = TRUE)) {
        return(yaml::read_yaml(path))
    }
    stop("Package 'yaml' is required to load HKOCA cell-type colors.")
}

.default_colors_yaml <- function() {
    file.path(getwd(), "celltype_colors.yaml")
}

load_hkoca_celltype_palettes <- function(yaml_path = NULL) {
    path <- yaml_path %||% Sys.getenv("HKOCA_CELLTYPE_COLORS", unset = "")
    if (!nzchar(path)) {
        script_dir <- tryCatch({
            cmd <- commandArgs(trailingOnly = FALSE)
            fa <- grep("^--file=", cmd, value = TRUE)
            if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)) else NULL
        }, error = function(e) NULL)
        if (!is.null(script_dir)) {
            candidate <- file.path(script_dir, "..", "..", "config", "celltype_colors.yaml")
            if (file.exists(candidate)) path <- candidate
        }
    }
    if (!nzchar(path) || !file.exists(path)) {
        stop(sprintf(
            "HKOCA cell-type colors YAML not found (set HKOCA_CELLTYPE_COLORS or pass yaml_path). Tried: %s",
            path
        ))
    }
    raw <- .load_yaml_safe(path)
    list(
        fallback = raw$fallback %||% "#BBBBBB",
        level1 = unlist(raw$level1, use.names = TRUE),
        level2 = unlist(raw$level2, use.names = TRUE),
        level3 = unlist(raw$level3, use.names = TRUE),
        combined = c(
            unlist(raw$level1, use.names = TRUE),
            unlist(raw$level2, use.names = TRUE),
            unlist(raw$level3, use.names = TRUE)
        ),
        path = path
    )
}

.normalize_label <- function(x) trimws(as.character(x))

.lookup_celltype_color <- function(label, palette, fallback = "#BBBBBB") {
    key <- .normalize_label(label)
    if (key %in% names(palette)) return(unname(palette[[key]]))
    under <- gsub(" ", "_", key, fixed = TRUE)
    if (under %in% names(palette)) return(unname(palette[[under]]))
    spaced <- gsub("_", " ", key, fixed = TRUE)
    if (spaced %in% names(palette)) return(unname(palette[[spaced]]))
    fallback
}

palette_for_level <- function(palettes, level = c("combined", "level1", "level2", "level3")) {
    level <- match.arg(level)
    palettes[[level]]
}

colors_for_labels <- function(labels, level = "combined", palettes = NULL, yaml_path = NULL) {
    if (is.null(palettes)) palettes <- load_hkoca_celltype_palettes(yaml_path)
    pal <- palette_for_level(palettes, level)
    fb <- palettes$fallback %||% "#BBBBBB"
    stats::setNames(
        vapply(labels, .lookup_celltype_color, character(1), palette = pal, fallback = fb),
        labels
    )
}

infer_palette_level <- function(metadata_col) {
    col <- as.character(metadata_col)
    if (grepl("^Level_1", col) || col == "Level_1") return("level1")
    if (grepl("^Level_2", col) || col == "Level_2") return("level2")
    if (grepl("^Level_3", col) || col == "Level_3") return("level3")
    if (grepl("^Level_latest", col) || col %in% c("Level_latest", "Level_3_latest")) return("combined")
    "combined"
}

scale_color_hkoca_celltypes <- function(level = "combined", palettes = NULL, yaml_path = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("ggplot2 is required for scale_color_hkoca_celltypes().")
    if (is.null(palettes)) palettes <- load_hkoca_celltype_palettes(yaml_path)
    pal <- palette_for_level(palettes, level)
    ggplot2::scale_color_manual(values = pal, ...)
}

scale_fill_hkoca_celltypes <- function(level = "combined", palettes = NULL, yaml_path = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
        stop("ggplot2 is required for scale_fill_hkoca_celltypes().")
    if (is.null(palettes)) palettes <- load_hkoca_celltype_palettes(yaml_path)
    pal <- palette_for_level(palettes, level)
    ggplot2::scale_fill_manual(values = pal, ...)
}

# dittoSeq: named vector for color.panel
ditto_colors_for_meta <- function(seurat_obj, metadata_col, palettes = NULL, yaml_path = NULL) {
    level <- infer_palette_level(metadata_col)
    if (is.null(palettes)) palettes <- load_hkoca_celltype_palettes(yaml_path)
    labels <- sort(unique(as.character(seurat_obj[[metadata_col, drop = TRUE]])))
    colors_for_labels(labels, level = level, palettes = palettes)
}
