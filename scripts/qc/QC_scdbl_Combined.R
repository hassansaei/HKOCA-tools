# =============================================================================
# Single-Cell QC Pipeline — Full Integration
# QC Filtering  →  Doublet Detection (scDblFinder)
#
# Description:
#   Stage 1 (QC): Loads raw Seurat RDS files, applies CSV-driven thresholds
#     (nFeature, nCount, percent.mito), saves filtered RDS + audit PDF +
#     full 7-panel dashboard (bar, intensity, inflation, line plots, ridge plots).
#   Stage 2 (Doublet): Reads QC-filtered RDS, runs scDblFinder (batch-aware,
#     optional per-sample DBR via YAML), saves singlet-only RDS + audit PDF +
#     summary CSV.
#   Stage 3 (Summary): Integrated QC+Doublet end-to-end summary CSV + plots.
#   Both stages share a single config file, a single logger, and the same CLI.
#
# Config file (DCF format — qc_config.dcf):
#   rds_dir             = /path/to/raw_rds
#   decisions_csv       = /path/to/QC_Decisions.csv
#   output_dir          = /path/to/output
#   summary_subdir      = Summary
#   filtered_subdir     = qc_filtered_rds
#   doublet_subdir      = doublet_filtered_rds
#   log_file            = logs/qc_pipeline.log
#   rds_pattern         = \.rds$
#   recursive_discovery = FALSE
#   skip_noisy_datasets = TRUE
#   noisy_datasets      = file1.rds,file2.rds
#   dbl_batch_col       = sample_id
#   dbl_min_count       = 100
#   dbl_min_feature     = 50
#   dbl_umap_dims       = 20
#   dbl_default_platform = 10x
#   sample.<name>.platform = 10x | dropseq | indrops | ...
#   sample.<name>.dbr      = 0.08   (omit for 10x auto-estimation)
#
# Usage:
#   Rscript qc_pipeline.R [--config PATH] [--stage all|qc|doublet] [OPTIONS]
#
# CLI options:
#   --config PATH           Config DCF file (default: qc_config.dcf next to script)
#   --stage STR             all | qc | doublet  (default: all)
#   --rds_dir PATH          Override raw RDS input directory
#   --decisions_csv PATH    Override QC decisions CSV
#   --output_dir PATH       Override output base directory
#   --summary_subdir NAME   Override QC summary subdirectory name
#   --filtered_subdir NAME  Override filtered RDS subdirectory name
#   --doublet_subdir NAME   Override doublet output subdirectory name
#   --log_file PATH         Override log file path
#   --recursive_discovery   Enable recursive RDS file search
#   --rds_pattern REGEX     File name regex for RDS discovery
#   --force_overwrite       Force reprocessing of existing outputs
#   --help                  Print this message and exit
# =============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 0 — Package Loading
# ══════════════════════════════════════════════════════════════════════════════

# NOTE: Dependencies must be installed externally via Conda or renv.
# See environment.yml and README.md.

load_packages_safely <- function() {
    required_packages <- c(
        "Seurat", "MASS", "ggplot2", "patchwork", "cowplot", 
        "scales", "reshape2", "ggridges", "SingleCellExperiment", 
        "scDblFinder", "yaml", "dplyr", "jsonlite", "digest"
    )
    
    missing <- character(0)
    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            missing <- c(missing, pkg)
        }
    }
    
    if (length(missing) > 0) {
        stop(sprintf("Missing required packages: %s\nPlease construct the environment using 'conda env create -f environment.yml' and activate it.", paste(missing, collapse = ", ")))
    }
    
    suppressPackageStartupMessages({
        library(Seurat)
        library(MASS)
        library(ggplot2)
        library(patchwork)
        library(cowplot)
        library(scales)
        library(reshape2)
        library(ggridges)
        library(SingleCellExperiment)
        library(scDblFinder)
        library(yaml)
        library(dplyr)
        library(jsonlite)
        library(digest)
    })
}

tryCatch({
    load_packages_safely()
}, error = function(e) {
    message(e$message)
    quit(save = "no", status = 1)
})

set.seed(1234)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Shared Runtime Utilities
# ══════════════════════════════════════════════════════════════════════════════

# ── 1.1 Script directory ──────────────────────────────────────────────────────
.get_script_dir <- function() {
    cmd      <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd, value = TRUE)
    if (length(file_arg) > 0)
        return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
    getwd()
}

# ── 1.2 CLI argument parser ───────────────────────────────────────────────────
.parse_cli_args <- function(args) {
    out <- list(); i <- 1L
    while (i <= length(args)) {
        arg <- args[[i]]
        if (startsWith(arg, "--")) {
            key <- sub("^--", "", arg); val <- "TRUE"
            if (grepl("=", key, fixed = TRUE)) {
                kv  <- strsplit(key, "=", fixed = TRUE)[[1]]
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

# ── 1.3 Type coercions ────────────────────────────────────────────────────────
.as_bool <- function(x, default = FALSE) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x == "") return(default)
    tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

.parse_csv_values <- function(x) {
    if (is.null(x) || is.na(x) || trimws(x) == "") return(character(0))
    trimws(unlist(strsplit(x, ",", fixed = TRUE)))
}

# ── 1.4 Path resolution ───────────────────────────────────────────────────────
.resolve_path <- function(path_value, base_dir = getwd()) {
    if (is.null(path_value) || path_value == "") return(path_value)
    if (grepl("^~",             path_value)) path_value <- path.expand(path_value)
    if (grepl("^(/|[A-Za-z]:)", path_value))
        return(normalizePath(path_value, mustWork = FALSE))
    normalizePath(file.path(base_dir, path_value), mustWork = FALSE)
}

# ── 1.5 DCF config reader ─────────────────────────────────────────────────────
.read_config_dcf <- function(config_path) {
    if (!file.exists(config_path))
        stop(sprintf("Config file not found: %s", config_path))
    
    # readLines to avoid encoding/formatting issues with read.dcf
    # Strip comments and ensure we have clean key:value pairs
    lines <- readLines(config_path, warn = FALSE)
    clean_lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
    
    if (length(clean_lines) == 0)
        stop(sprintf("Config file is empty or contains only comments: %s", config_path))
        
    # Create a clean temporary DCF file for read.dcf
    tmp_dcf <- tempfile(fileext = ".dcf")
    writeLines(clean_lines, tmp_dcf)
    on.exit(unlink(tmp_dcf))
    
    as.list(read.dcf(tmp_dcf, all = TRUE)[1, , drop = FALSE])
}

# ── 1.6 Null-coalesce operator ────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── 1.7 Logger ────────────────────────────────────────────────────────────────
.LOG_FILE <- NULL

.init_logger <- function(log_file_path) {
    .LOG_FILE <<- log_file_path
    log_dir    <- dirname(.LOG_FILE)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    cat("", file = .LOG_FILE, append = TRUE)
}

.log_message <- function(level, ...) {
    msg  <- paste0(...)
    line <- sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
    cat(line, "\n")
    if (!is.null(.LOG_FILE)) cat(line, "\n", file = .LOG_FILE, append = TRUE)
}

log_info  <- function(...) .log_message("INFO",  ...)
log_warn  <- function(...) .log_message("WARN",  ...)
log_error <- function(...) .log_message("ERROR", ...)

# ── 1.8 RDS file discovery ────────────────────────────────────────────────────
discover_rds_files <- function(root_dir, pattern = "\\.rds$", recursive = FALSE) {
    if (!dir.exists(root_dir))
        stop(sprintf("RDS directory does not exist: %s", root_dir))

    paths <- list.files(root_dir, pattern = pattern, full.names = TRUE,
                        recursive = recursive, ignore.case = TRUE)

    if (length(paths) == 0)
        stop(sprintf("No RDS files found in %s with pattern '%s'", root_dir, pattern))

    mapping <- list()
    for (fp in sort(paths)) {
        key <- sub(pattern, "", basename(fp), ignore.case = TRUE)
        if (!is.null(mapping[[key]])) {
            log_warn(sprintf("Duplicate key '%s' — keeping first, skipping: %s", key, fp))
            next
        }
        mapping[[key]] <- normalizePath(fp, mustWork = FALSE)
    }
    mapping
}

# ── 1.9 Usage printer ─────────────────────────────────────────────────────────
.print_usage <- function(default_config) {
    cat("Usage:\n")
    cat("  Rscript qc_pipeline.R [--config PATH] [--stage all|qc|doublet] [OPTIONS]\n\n")
    cat("Options:\n")
    cat(sprintf("  --config PATH            Config DCF file (default: %s)\n", default_config))
    cat("  --stage STR              all | qc | doublet  (default: all)\n")
    cat("  --rds_dir PATH           Override raw RDS input directory\n")
    cat("  --decisions_csv PATH     Override QC decisions CSV\n")
    cat("  --output_dir PATH        Override output base directory\n")
    cat("  --summary_subdir NAME    Override QC summary subdirectory\n")
    cat("  --filtered_subdir NAME   Override QC filtered RDS subdirectory\n")
    cat("  --doublet_subdir NAME    Override doublet output subdirectory\n")
    cat("  --log_file PATH          Override log file path\n")
    cat("  --recursive_discovery    Enable recursive RDS file search\n")
    cat("  --rds_pattern REGEX      RDS file name regex\n")
    cat("  --force_overwrite        Force reprocessing of existing outputs\n")
    cat("  --help                   Show this message\n")
}

# ── 1.10 Parameter Metadata Checkpointing ────────────────────────────────────
# Saves a hidden .<file>.json metadata file to track thresholds.
# If thresholds in CSV/CLI change, the script automatically re-runs.
.get_meta_path <- function(target_rds) {
    file.path(dirname(target_rds), paste0(".", basename(target_rds), ".meta.json"))
}

.should_skip_run <- function(target_rds, current_params, force_flag) {
    if (.as_bool(force_flag, FALSE)) return(FALSE)
    if (!file.exists(target_rds) || file.info(target_rds)$size == 0) return(FALSE)
    
    meta_path <- .get_meta_path(target_rds)
    if (!file.exists(meta_path)) return(FALSE)
    
    tryCatch({
        old_meta <- jsonlite::fromJSON(meta_path)
        # Compare current params vs stored params
        # We use digest to handle complex list comparisons reliably
        return(digest::digest(current_params) == old_meta$param_hash)
    }, error = function(e) return(FALSE))
}

.save_run_metadata <- function(target_rds, current_params) {
    meta_path <- .get_meta_path(target_rds)
    meta_data <- list(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        param_hash = digest::digest(current_params),
        params = current_params
    )
    jsonlite::write_json(meta_data, meta_path, auto_unbox = TRUE, pretty = TRUE)
}



# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Configuration & Path Setup
# ══════════════════════════════════════════════════════════════════════════════

# Run full pipeline only when executed as a script.
if (sys.nframe() == 0L || !exists(".read_config_dcf")) {

script_dir     <- .get_script_dir()
default_config <- file.path(script_dir, "qc_config.dcf")
cli_args       <- .parse_cli_args(commandArgs(trailingOnly = TRUE))

if (.as_bool(cli_args$help, FALSE)) { .print_usage(default_config); quit(save = "no", status = 0) }

config_path <- if (!is.null(cli_args$config)) .resolve_path(cli_args$config, getwd()) else default_config
cfg         <- .read_config_dcf(config_path)

# CLI overrides win over config file
for (nm in setdiff(names(cli_args), c("config", "help", "stage"))) cfg[[nm]] <- cli_args[[nm]]

# Stage selection
RUN_STAGE <- tolower(cli_args$stage %||% "all")
if (!RUN_STAGE %in% c("all", "qc", "doublet"))
    stop(sprintf("Invalid --stage '%s'. Choose: all | qc | doublet", RUN_STAGE))

RUN_QC      <- RUN_STAGE %in% c("all", "qc")
RUN_DOUBLET <- RUN_STAGE %in% c("all", "doublet")

# ── Resolve all paths ─────────────────────────────────────────────────────────
rds_dir         <- .resolve_path(cfg$rds_dir       %||% NULL,             getwd())
decisions_csv   <- .resolve_path(cfg$decisions_csv %||% NULL, getwd())
BASE_OUT_DIR    <- .resolve_path(cfg$output_dir    %||% NULL,                  getwd())

summary_subdir  <- if (!is.null(cfg$summary_subdir)  && nzchar(cfg$summary_subdir))  cfg$summary_subdir  else "Summary"
filtered_subdir <- if (!is.null(cfg$filtered_subdir) && nzchar(cfg$filtered_subdir)) cfg$filtered_subdir else "qc_filtered_rds"
doublet_subdir  <- if (!is.null(cfg$doublet_subdir)  && nzchar(cfg$doublet_subdir))  cfg$doublet_subdir  else "doublet_filtered_rds"

QC_DIR              <- file.path(BASE_OUT_DIR, summary_subdir)
FILTERED_DIR        <- file.path(BASE_OUT_DIR, filtered_subdir)
DOUBLET_DIR         <- file.path(BASE_OUT_DIR, doublet_subdir)
DOUBLET_SUMMARY_DIR <- file.path(DOUBLET_DIR, "Summary")
DOUBLET_SINGLET_DIR <- file.path(DOUBLET_DIR, "singlets")
DOUBLET_WITH_DOUBLETS_DIR <- file.path(DOUBLET_DIR, "with doublets")

default_log   <- file.path(BASE_OUT_DIR, "logs", "qc_pipeline.log")
log_file_path <- .resolve_path(
    if (!is.null(cfg$log_file) && nzchar(cfg$log_file)) cfg$log_file else default_log,
    getwd()
)

# ── scDblFinder tuning parameters ────────────────────────────────────────────
DBL_BATCH_COL      <- cfg$dbl_batch_col      %||% "sample_id"
DBL_MIN_COUNT      <- as.integer(cfg$dbl_min_count   %||% 100)
DBL_MIN_FEATURE    <- as.integer(cfg$dbl_min_feature %||% 50)
DBL_UMAP_DIMS      <- 1:as.integer(cfg$dbl_umap_dims %||% 20)
DBL_DEFAULT_PLATFORM <- tolower(cfg$dbl_default_platform %||% "10x")

# ── Per-sample DBR table parsed from DCF keys ─────────────────────────────────
# Keys follow the convention:
#   sample.<sample_name>.platform = 10x | dropseq | indrops | ...
#   sample.<sample_name>.dbr      = 0.08   (omit for 10x auto-estimation)
#
# Example in qc_config.dcf:
#   sample.my_sample_A.platform: 10x
#   sample.my_sample_A.dbr: 0.08
#   sample.my_sample_B.platform: dropseq
#   sample.my_sample_B.dbr: 0.05
#
# The parser scans all cfg keys that start with "sample." and builds a lookup
# list: DBL_SAMPLE_CFG[["my_sample_A"]] = list(platform="10x", dbr=0.08)
DBL_SAMPLE_CFG <- list()
sample_keys    <- grep("^sample\\.", names(cfg), value = TRUE)
for (sk in sample_keys) {
    parts <- strsplit(sk, "\\.")[[1]]   # ["sample", "<name>", "<field>"]
    if (length(parts) != 3) next
    sname <- parts[2]; field <- parts[3]
    if (is.null(DBL_SAMPLE_CFG[[sname]])) DBL_SAMPLE_CFG[[sname]] <- list()
    DBL_SAMPLE_CFG[[sname]][[field]] <- cfg[[sk]]
}

# ── Discovery settings ────────────────────────────────────────────────────────
rds_pattern         <- cfg$rds_pattern         %||% "\\.rds$"
recursive_discovery <- .as_bool(cfg$recursive_discovery %||% "FALSE", FALSE)

# ── Create output directories ─────────────────────────────────────────────────
for (d in c(QC_DIR, FILTERED_DIR, DOUBLET_DIR, DOUBLET_SUMMARY_DIR,
            DOUBLET_SINGLET_DIR, DOUBLET_WITH_DOUBLETS_DIR))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Initialise shared logger ──────────────────────────────────────────────────
.init_logger(log_file_path)
run_ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

log_info("══════════════════════════════════════════════════")
log_info(sprintf("QC PIPELINE START  [stage: %s]", RUN_STAGE))
log_info(sprintf("  Config       : %s", config_path))
log_info(sprintf("  Raw RDS dir  : %s", rds_dir))
log_info(sprintf("  Decisions CSV: %s", decisions_csv))
log_info(sprintf("  Output base  : %s", BASE_OUT_DIR))
log_info(sprintf("  Log file     : %s", log_file_path))
log_info("══════════════════════════════════════════════════")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Threshold Algorithm Implementations
# ══════════════════════════════════════════════════════════════════════════════

# ── Triangle threshold ────────────────────────────────────────────────────────
# BUG NOTE (older versions): which.max(distances) searched ALL bins, not just
# the region between peak and endpoint, producing absurd thresholds.
# FIX: restrict search_range to bins between peak and endpoint only.
.thresh_triangle <- function(x, sub = "higher", bins = 256) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(0)
    if (length(unique(x)) < 5) return(if (sub == "higher") max(x) else min(x))

    h        <- hist(x, breaks = bins, plot = FALSE)
    counts   <- h$counts; mids <- h$mids
    peak_idx <- which.max(counts)
    end_idx  <- if (sub == "higher") length(counts) else 1L

    if (peak_idx == end_idx) return(mids[peak_idx])

    x1 <- peak_idx; y1 <- counts[peak_idx]
    x2 <- end_idx;  y2 <- counts[end_idx]

    distances <- abs(
        (y2 - y1) * seq_along(counts) -
        (x2 - x1) * counts +
         x2 * y1  -  y2 * x1
    ) / sqrt((y2 - y1)^2 + (x2 - x1)^2)

    search_range <- if (sub == "higher") seq(peak_idx, length(counts)) else seq(1L, peak_idx)
    mids[search_range[which.max(distances[search_range])]]
}

# ── Rényi entropy threshold ───────────────────────────────────────────────────
.thresh_renyi <- function(x, bins = 256) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(0)
    if (length(unique(x)) < 5) return(max(x))

    h <- hist(x, breaks = bins, plot = FALSE)
    counts <- h$counts; mids <- h$mids
    if (sum(counts) == 0) return(max(x))
    probs <- counts / sum(counts)

    best_thresh <- mids[which.max(counts)]; best_score <- -Inf
    if (length(probs) < 3) return(best_thresh)

    for (i in seq(2L, length(probs) - 1L)) {
        p1 <- probs[1:i]; p2 <- probs[(i + 1):length(probs)]
        s1 <- sum(p1, na.rm = TRUE); s2 <- sum(p2, na.rm = TRUE)
        if (is.na(s1) || is.na(s2) || s1 == 0 || s2 == 0) next
        h1 <- -log(sum((p1 / s1)^2) + 1e-12)
        h2 <- -log(sum((p2 / s2)^2) + 1e-12)
        if ((h1 + h2) > best_score) { best_score <- h1 + h2; best_thresh <- mids[i] }
    }
    best_thresh
}

# ── KneePoint (Kneedle) threshold ────────────────────────────────────────────
.thresh_knee <- function(x, use_log = TRUE) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(0)
    if (length(unique(x)) < 5) return(median(x))

    x_sorted <- sort(x[x > 0], decreasing = TRUE)
    if (length(x_sorted) < 2) return(median(x))
    ranks <- seq_along(x_sorted)

    curve_x <- if (use_log) log10(ranks)        else as.numeric(ranks)
    curve_y <- if (use_log) log10(x_sorted + 1) else x_sorted

    x_norm <- (curve_x - min(curve_x)) / (max(curve_x) - min(curve_x))
    y_norm <- (curve_y - min(curve_y)) / (max(curve_y) - min(curve_y))

    x_sorted[which.max((x_norm + y_norm - 1) / sqrt(2))]
}

# ── Quantile threshold ────────────────────────────────────────────────────────
.thresh_quantile <- function(x, q) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(0)
    as.numeric(quantile(x, probs = q, na.rm = TRUE))
}

# ── Dispatcher: CSV method string → numeric threshold ─────────────────────────
get_thresh <- function(method_str, values, bound = "lower") {
    method_str <- trimws(as.character(method_str %||% ""))
    values <- values[is.finite(values)]
    if (method_str == "lower_tri") return(.thresh_triangle(values, sub = "lower"))
    if (method_str == "upper_tri") return(.thresh_triangle(values, sub = "higher"))
    if (method_str == "renyi")     return(.thresh_renyi(values))
    if (method_str == "knee")      return(.thresh_knee(values))
    if (method_str == "mad3") {
        m <- median(values, na.rm = TRUE); d <- mad(values, na.rm = TRUE)
        return(if (bound == "upper") m + 3 * d else m - 3 * d)
    }
    if (method_str == "none") return(if (bound == "upper") Inf else -Inf)
    if (grepl("^manual_", method_str)) return(as.numeric(gsub("manual_", "", method_str)))
    log_warn(sprintf("Unrecognised threshold method '%s' — defaulting to no filter", method_str))
    ifelse(bound == "upper", Inf, -Inf)
}

r.shorten_label <- function(x) {
    x <- as.character(x)
    x <- gsub("d10_1016_j_",              "", x)
    x <- gsub("d10_1126_sciadv_",         "", x)
    x <- gsub("d10_1038_",                "", x)
    x <- gsub("dno_doi_kidney_organoid_", "kidney_", x)
    x
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — QC Plot Builders
# ══════════════════════════════════════════════════════════════════════════════

# ── 2-D kernel density (scatter colour) ──────────────────────────────────────
point_density <- function(x, y, n = 200) {
    ok  <- is.finite(x) & is.finite(y)
    out <- rep(NA_real_, length(x))
    if (sum(ok) < 10 || var(x[ok]) == 0 || var(y[ok]) == 0) return(out)
    kd      <- MASS::kde2d(x[ok], y[ok], n = n)
    ix      <- findInterval(x[ok], kd$x, all.inside = TRUE)
    iy      <- findInterval(y[ok], kd$y, all.inside = TRUE)
    out[ok] <- kd$z[cbind(ix, iy)]
    out
}

# ── Density scatter ───────────────────────────────────────────────────────────
make_scatter <- function(x, y, x_label, y_label, title,
                         max_x, max_y,
                         cutoff_l = NULL, cutoff_h = NULL,
                         show_cutoffs = TRUE) {
    dens        <- point_density(x, y)
    df          <- data.frame(X1 = x, X2 = y, dens = dens)
    finite_dens <- dens[is.finite(dens)]
    dens_limit  <- if (length(finite_dens) == 0 || max(finite_dens) == 0) 1 else {
        dl <- quantile(finite_dens, 0.99)
        if (!is.finite(dl) || dl <= 0) max(finite_dens) else dl
    }

    use_raster <- nrow(df) > 50000 && requireNamespace("ggrastr", quietly = TRUE)
    geom_pts   <- if (use_raster) ggrastr::geom_point_rast(size = 0.3, raster.dpi = 300)
                  else            geom_point(size = 0.3)

    p <- ggplot(df, aes(x = X1, y = X2, color = dens)) +
        geom_pts +
        scale_x_continuous(limits = c(0, max_x), oob = scales::squish) +
        scale_y_continuous(limits = c(0, max_y), oob = scales::squish) +
        scale_color_viridis_c(option = "magma", na.value = "grey70",
                              limits = c(0, dens_limit), oob = scales::squish) +
        labs(x = x_label, y = y_label, color = "Density") +
        ggtitle(title) +
        theme_bw() +
        theme(
            plot.background  = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = "grey80"),
            panel.grid.major = element_line(color = "grey92"),
            panel.grid.minor = element_blank(),
            plot.title       = element_text(size = 9, face = "bold"),
            axis.title       = element_text(size = 8),
            axis.text        = element_text(size = 7),
            legend.position  = "bottom",
            legend.title     = element_text(size = 7, face = "bold"),
            legend.text      = element_text(size = 6)
        ) +
        guides(color = guide_colourbar(barwidth = 8, barheight = 0.5))

    if (show_cutoffs && !is.null(cutoff_l) && is.finite(cutoff_l))
        p <- p + geom_hline(yintercept = cutoff_l, linewidth = 0.6, color = "red", linetype = "dashed")
    if (show_cutoffs && !is.null(cutoff_h) && is.finite(cutoff_h))
        p <- p + geom_hline(yintercept = cutoff_h, linewidth = 0.6, color = "red", linetype = "dashed")
    p
}

# ── Violin (avoids Seurat v5 S4 slot bug) ────────────────────────────────────
make_violin <- function(obj, feat, title, max_y, cutoff_l = NULL, cutoff_h = NULL) {
    df <- data.frame(x = "cells", y = as.numeric(obj@meta.data[[feat]]))
    p  <- ggplot(df, aes(x = x, y = y)) +
        geom_violin(fill = "#4C8BE0", color = NA, alpha = 0.8, scale = "width") +
        geom_jitter(width = 0.2, size = 0.1, alpha = 0.15, color = "grey30") +
        scale_y_continuous(limits = c(0, max_y), oob = scales::squish) +
        labs(x = NULL, y = feat, title = title) +
        theme_bw() +
        theme(
            plot.background  = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = "grey80"),
            axis.text.x      = element_text(hjust = 0.5, angle = 0),
            legend.position  = "none",
            plot.title       = element_text(size = 9, face = "bold")
        )
    if (!is.null(cutoff_l) && is.finite(cutoff_l))
        p <- p + geom_hline(yintercept = cutoff_l, linewidth = 0.5, color = "red")
    if (!is.null(cutoff_h) && is.finite(cutoff_h))
        p <- p + geom_hline(yintercept = cutoff_h, linewidth = 0.5, color = "red")
    p
}

# ── Audit density (PDF before/after pages) ───────────────────────────────────
plot_dist <- function(meta, col_name, title, lower_val, upper_val,
                      is_log = TRUE, fill_color = "steelblue") {
    p <- ggplot(meta, aes(x = if (is_log) log1p(.data[[col_name]]) else .data[[col_name]])) +
        geom_density(fill = fill_color, alpha = 0.5) +
        labs(title = title,
             x = ifelse(is_log, paste0("log1p(", col_name, ")"), col_name),
             y = "Density") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold", size = 10))
    if (is.finite(lower_val))
        p <- p + geom_vline(xintercept = if (is_log) log1p(lower_val) else lower_val,
                            color = "blue", linetype = "dashed", linewidth = 1)
    if (is.finite(upper_val))
        p <- p + geom_vline(xintercept = if (is_log) log1p(upper_val) else upper_val,
                            color = "red",  linetype = "dashed", linewidth = 1)
    p
}

# ── Summary stat helpers ──────────────────────────────────────────────────────
s_median <- function(x) round(median(as.numeric(x), na.rm = TRUE), 2)
s_mean   <- function(x) round(mean(as.numeric(x),   na.rm = TRUE), 2)
s_max    <- function(x) round(max(as.numeric(x),    na.rm = TRUE), 2)
s_min    <- function(x) round(min(as.numeric(x),    na.rm = TRUE), 2)

make_summary_row <- function(obj, fobj, nm) {
    cells_raw      <- ncol(obj)
    cells_filtered <- ncol(fobj)
    cells_removed  <- cells_raw - cells_filtered
    pct_removed    <- if (cells_raw > 0) round(100 * cells_removed / cells_raw, 1) else NA_real_

    data.frame(
        sample            = nm,
        cells_raw         = cells_raw,
        cells_filtered    = cells_filtered,
        cells_removed     = cells_removed,
        pct_cells_removed = pct_removed,
        stringsAsFactors  = FALSE
    )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — STAGE 1: QC Filtering
# ══════════════════════════════════════════════════════════════════════════════

if (RUN_QC) {

    log_info("──────────────────────────────────────────────────")
    log_info("STAGE 1: QC FILTERING")
    log_info("──────────────────────────────────────────────────")

    # ── 5.1 Discover raw RDS files ────────────────────────────────────────────
    rds_file_map <- discover_rds_files(rds_dir, pattern = rds_pattern, recursive = recursive_discovery)
    sample_names <- sort(names(rds_file_map))
    log_info(sprintf("Discovered %d RDS files.", length(sample_names)))
    log_info(sprintf("Datasets: %s", paste(sample_names, collapse = ", ")))

    # ── 5.2 Load QC decisions CSV ─────────────────────────────────────────────
    if (!file.exists(decisions_csv))
        stop(sprintf("QC decisions CSV not found: %s", decisions_csv))

    qc_decisions <- read.csv(decisions_csv, stringsAsFactors = FALSE,colClasses = "character", strip.white = TRUE)
    required_cols <- c(
        "Dataset_Name", "Lower_Feature_Method", "Upper_Feature_Method",
        "Lower_Count_Method", "Upper_Count_Method", "Upper_Mito_Method"
    )
    missing_cols <- setdiff(required_cols, colnames(qc_decisions))
    if (length(missing_cols) > 0)
        stop(sprintf("QC decisions CSV missing columns: %s", paste(missing_cols, collapse = ", ")))

    # ── 5.3 Dependency guard ──────────────────────────────────────────────────
    required_vars <- c("rds_dir", "decisions_csv", "rds_file_map", "QC_DIR", "FILTERED_DIR",
                       "s_median", "s_mean", "s_max")
    missing_vars  <- required_vars[!sapply(required_vars, exists)]
    if (length(missing_vars) > 0)
        stop(sprintf("Missing required variables: %s", paste(missing_vars, collapse = ", ")))

    # ── 5.4 Noisy dataset exclusions ─────────────────────────────────────────
    # Hardcoded fallback list (matches original Cell 4 defaults)
    default_noisy_datasets <- c()

    skip_noisy <- .as_bool(cfg$skip_noisy_datasets %||% "TRUE", TRUE)
    noisy_list <- if (!is.null(cfg$noisy_datasets) && nzchar(cfg$noisy_datasets))
                      .parse_csv_values(cfg$noisy_datasets)
                  else
                      default_noisy_datasets

    if (skip_noisy && length(noisy_list) > 0) {
        qc_decisions <- qc_decisions[!basename(qc_decisions$Dataset_Name) %in% noisy_list, ]
        log_info(sprintf("Skipping %d noisy dataset(s): %s",
                         length(noisy_list), paste(noisy_list, collapse = ", ")))
    }

    configured_datasets <- intersect(
        gsub("\\.rds$", "", basename(qc_decisions$Dataset_Name), ignore.case = TRUE),
        sample_names
    )

    if (length(configured_datasets) == 0)
        log_warn("No configured datasets overlap with discovered RDS files.")

    log_info(sprintf("Loaded %d QC decision rows. Running %d dataset(s).",
                     nrow(qc_decisions), length(configured_datasets)))

    # ── 5.5 Open audit PDF (only if at least one dataset is runnable) ───────
    qc_pdf_path <- file.path(QC_DIR, "QC_Before_After_Report.pdf")
    qc_pdf_open <- FALSE
    if (length(configured_datasets) > 0) {
        pdf(qc_pdf_path, width = 14, height = 6)
        qc_pdf_open <- TRUE
    } else {
        log_warn("No configured datasets to process — QC audit PDF will not be generated.")
    }

    summary_rows <- list()

    # ── 5.6 Per-dataset QC loop ───────────────────────────────────────────────
    for (i in seq_len(nrow(qc_decisions))) {

        nm       <- gsub("\\.rds$", "", basename(qc_decisions$Dataset_Name[i]), ignore.case = TRUE)
        rds_path <- rds_file_map[[nm]] # Retrieve full path from map
        row_cfg  <- as.list(qc_decisions[i, , drop = FALSE])
        row_cfg  <- lapply(row_cfg, function(x) trimws(as.character(x)))

        if (is.null(rds_path) || !file.exists(rds_path)) {
            log_warn(sprintf("[MISSING] %s — not found in discovered RDS files", nm)); next
        }

        # Checkpoint / Resume logic (Parameter-Aware)
        filtered_path <- file.path(FILTERED_DIR, paste0(nm, "_filtered.rds"))
        current_qc_params <- as.list(row_cfg)

        if (.should_skip_run(filtered_path, current_qc_params, cfg$force_overwrite)) {
            log_info(sprintf("[QC SKIPPED] %s already processed with identical thresholds.", nm))
            skip_row <- tryCatch({
                fobj <- readRDS(filtered_path)
                obj  <- readRDS(rds_path) # Need raw to calculate 'cells_removed', etc.
                out <- make_summary_row(obj, fobj, nm)
                rm(obj, fobj); gc()
                out
            }, error = function(e) {
                log_warn(sprintf("Could not restore summary stats for skipped %s: %s", nm, e$message))
                NULL
            })

            if (!is.null(skip_row)) summary_rows[[length(summary_rows) + 1L]] <- skip_row
            next
        }

        row_result <- tryCatch({

            log_info(sprintf("[QC] Processing: %s", nm))
            dir.create(file.path(QC_DIR, nm), recursive = TRUE, showWarnings = FALSE)

            # Load and slim
            # Load and slim
            # ── Version-safe object reconstruction ──────────────────────────
            # Load raw without triggering SeuratObject v5 update dispatch
            obj_raw <- readRDS(rds_path)
            log_info(sprintf("  Raw class: %s | Assay class: %s",
                             class(obj_raw),
                             class(obj_raw[["RNA"]])))

            # Extract everything we need directly from slots,
            # bypassing all S4 method dispatch
            cm       <- slot(obj_raw[["RNA"]], "counts")
            meta_in  <- as.data.frame(obj_raw@meta.data)
            barcodes <- colnames(cm)
            genes    <- rownames(cm)
            log_info(sprintf("  Direct slot access OK — %d cells, %d genes",
                             length(barcodes), length(genes)))

            # Build a brand-new minimal Seurat object from scratch
            # This completely avoids UpdateSeuratObject()
            obj <- CreateSeuratObject(counts = cm, meta.data = meta_in)
            rm(obj_raw, cm); gc()
            log_info("  Fresh Seurat object created successfully")

            # Ensure QC metrics exist on the fresh object
            cm <- GetAssayData(obj, layer = "counts")
            if (is.null(obj$nCount_RNA)   || all(obj$nCount_RNA   == 0))
                obj$nCount_RNA   <- Matrix::colSums(cm)
            if (is.null(obj$nFeature_RNA) || all(obj$nFeature_RNA == 0))
                obj$nFeature_RNA <- Matrix::colSums(cm > 0)

            # Manual percent calculation
            gene_names       <- rownames(obj)
            mito_genes       <- grep("^MT-",       gene_names, value = TRUE)
            ribo_genes       <- grep("^(RPL|RPS)", gene_names, value = TRUE, perl = TRUE)
            total_counts     <- Matrix::colSums(cm)
            total_counts[total_counts == 0] <- NA_real_
            obj$percent.mito <- if (length(mito_genes) > 0)
                Matrix::colSums(cm[mito_genes, , drop = FALSE]) / total_counts * 100
                else rep(0, ncol(obj))
            obj$percent.ribo <- if (length(ribo_genes) > 0)
                Matrix::colSums(cm[ribo_genes, , drop = FALSE]) / total_counts * 100
                else rep(0, ncol(obj))

            
            

            # Remove empty barcodes
           # Remove empty barcodes
            nCount_vec   <- as.numeric(obj@meta.data[["nCount_RNA"]])
            nFeature_vec <- as.numeric(obj@meta.data[["nFeature_RNA"]])
            nonzero_keep <- colnames(obj)[nCount_vec > 0 & nFeature_vec > 0]
            if (length(nonzero_keep) < ncol(obj)) obj <- subset(obj, cells = nonzero_keep)

            meta_raw <- obj[[]]

            # Thresholds from CSV
            feat_lower  <- get_thresh(row_cfg$Lower_Feature_Method, meta_raw$nFeature_RNA, "lower")
            feat_upper  <- get_thresh(row_cfg$Upper_Feature_Method, meta_raw$nFeature_RNA, "upper")
            count_lower <- get_thresh(row_cfg$Lower_Count_Method,   meta_raw$nCount_RNA,   "lower")
            count_upper <- get_thresh(row_cfg$Upper_Count_Method,   meta_raw$nCount_RNA,   "upper")
            mito_upper  <- get_thresh(row_cfg$Upper_Mito_Method,    meta_raw$percent.mito, "upper")

            ribo_vals  <- meta_raw$percent.ribo
            ribo_lower <- if (median(ribo_vals, na.rm = TRUE) >= 10) quantile(ribo_vals, 0.05, na.rm = TRUE) else 0
            ribo_upper <- max(ribo_vals, na.rm = TRUE)

            # Filter
            cells_raw      <- ncol(obj); genes_raw <- nrow(obj)
            keep_cells <- colnames(obj)[
                as.numeric(obj@meta.data[["nFeature_RNA"]]) > feat_lower  &
                as.numeric(obj@meta.data[["nFeature_RNA"]]) < feat_upper  &
                as.numeric(obj@meta.data[["nCount_RNA"]])   > count_lower &
                as.numeric(obj@meta.data[["nCount_RNA"]])   < count_upper &
                as.numeric(obj@meta.data[["percent.mito"]]) < mito_upper
            ]
            fobj <- subset(obj, cells = keep_cells)
            cells_filtered <- ncol(fobj)
            cells_removed  <- cells_raw - cells_filtered
            pct_removed    <- round(100 * cells_removed / cells_raw, 1)

            # Axis maxima — raw drives scale so filtered plots are comparable
            max.ct <- max(meta_raw$nCount_RNA,   na.rm = TRUE) * 1.1
            max.ft <- max(meta_raw$nFeature_RNA, na.rm = TRUE) * 1.1
            max.mt <- max(meta_raw$percent.mito, na.rm = TRUE) * 1.1
            max.rb <- max(meta_raw$percent.ribo, na.rm = TRUE) * 1.1

            feat_list <- c("nFeature_RNA", "percent.mito", "percent.ribo", "nCount_RNA")
            cutoffs   <- list(
                nFeature_RNA = c(lower = feat_lower,  upper = feat_upper),
                percent.mito = c(lower = 0,           upper = mito_upper),
                percent.ribo = c(lower = ribo_lower,  upper = ribo_upper),
                nCount_RNA   = c(lower = count_lower, upper = count_upper)
            )
            max_map <- c(nFeature_RNA = max.ft, percent.mito = max.mt,
                         percent.ribo = max.rb, nCount_RNA   = max.ct)

            # 16 PNGs per dataset (raw + filtered × violin + scatter × 4 features)
            for (feat in feat_list) {
                cutoff.l <- cutoffs[[feat]]["lower"]
                cutoff.h <- cutoffs[[feat]]["upper"]
                max.y    <- max_map[feat]

                if (feat != "nCount_RNA") {
                    x_vals <- as.numeric(obj$nCount_RNA);    y_vals <- as.numeric(obj@meta.data[[feat]]); x_lbl <- "nCount_RNA"
                    xf     <- as.numeric(fobj$nCount_RNA);   yf     <- as.numeric(fobj@meta.data[[feat]])
                } else {
                    x_vals <- as.numeric(obj$nFeature_RNA);  y_vals <- as.numeric(obj$nCount_RNA);        x_lbl <- "nFeature_RNA"
                    xf     <- as.numeric(fobj$nFeature_RNA); yf     <- as.numeric(fobj$nCount_RNA)
                }

                p.v   <- make_violin(obj,  feat, nm, max.y, cutoff.l, cutoff.h)
                p.s   <- make_scatter(x_vals, y_vals, x_lbl, feat, nm, max.ct, max.y, cutoff.l, cutoff.h, show_cutoffs = TRUE)
                p.v.f <- make_violin(fobj, feat, nm, max.y)
                p.s.f <- make_scatter(xf, yf, x_lbl, feat, nm, max.ct, max.y, NULL, NULL, show_cutoffs = FALSE)

                ggsave(file.path(QC_DIR, nm, paste0("Violin_qc.raw.",      feat, ".png")), p.v,   width = 6, height = 5, dpi = 150, bg = "white")
                ggsave(file.path(QC_DIR, nm, paste0("Plot_qc.raw.",        feat, ".png")), p.s,   width = 6, height = 5, dpi = 150, bg = "white")
                ggsave(file.path(QC_DIR, nm, paste0("Violin_qc.filtered.", feat, ".png")), p.v.f, width = 6, height = 5, dpi = 150, bg = "white")
                ggsave(file.path(QC_DIR, nm, paste0("Plot_qc.filtered.",   feat, ".png")), p.s.f, width = 6, height = 5, dpi = 150, bg = "white")
            }
            log_info(sprintf("  PNGs: 16 files → %s", file.path(QC_DIR, nm)))

            # PDF audit — BEFORE page
            p1_b <- plot_dist(meta_raw, "nFeature_RNA", paste("Raw nFeature:", nm), feat_lower,  feat_upper,  TRUE,  "lightblue")
            p2_b <- plot_dist(meta_raw, "nCount_RNA",   paste("Raw nCount:",   nm), count_lower, count_upper, TRUE,  "lightgreen")
            p3_b <- plot_dist(meta_raw, "percent.mito", paste("Raw Mito %:",   nm), -Inf,        mito_upper,  FALSE, "lightcoral")
            print((p1_b | p2_b | p3_b) + plot_annotation(
                title    = paste("BEFORE QC:", nm),
                subtitle = sprintf("Lower Feat: %s | Upper Feat: %s | Upper Mito: %s",
                                   row_cfg$Lower_Feature_Method,
                                   row_cfg$Upper_Feature_Method,
                                   row_cfg$Upper_Mito_Method)))

            # PDF audit — AFTER page
            meta_filt <- fobj[[]]
            p1_a <- plot_dist(meta_filt, "nFeature_RNA", paste("Filtered nFeature:", nm), -Inf, Inf, TRUE,  "lightblue")
            p2_a <- plot_dist(meta_filt, "nCount_RNA",   paste("Filtered nCount:",   nm), -Inf, Inf, TRUE,  "lightgreen")
            p3_a <- plot_dist(meta_filt, "percent.mito", paste("Filtered Mito %:",   nm), -Inf, Inf, FALSE, "lightcoral")
            print((p1_a | p2_a | p3_a) + plot_annotation(
                title = paste("AFTER QC:", nm, sprintf("(%.1f%% cells removed)", pct_removed))))

            # Save filtered RDS
            filtered_path <- file.path(FILTERED_DIR, paste0(nm, "_filtered.rds"))
            saveRDS(fobj, filtered_path)
            .save_run_metadata(filtered_path, current_qc_params)

            # Inflation diagnostic
            inflation      <- round(s_median(fobj$nCount_RNA) / s_median(obj$nCount_RNA), 2)
            inflation_flag <- if (inflation > 1.5)   " \u26a0 OVER-FILTERING"
                              else if (inflation < 1.05) " \u26a0 CHECK: very light filtering"
                              else ""
            log_info(sprintf("  Saved  : %s", filtered_path))
            log_info(sprintf("  Cells  : %d \u2192 %d  (%.1f%% removed)", cells_raw, cells_filtered, pct_removed))
            log_info(sprintf("  nCount : median %g \u2192 %g  [inflation %.2fx]%s",
                             s_median(obj$nCount_RNA), s_median(fobj$nCount_RNA),
                             inflation, inflation_flag))

            out <- make_summary_row(obj, fobj, nm)

            rm(obj, fobj, meta_raw, meta_filt); gc()
            out
        }, error = function(e) {
            log_error(sprintf("[FAILED] %s: %s", nm, conditionMessage(e)))
            log_error(sprintf("[CALLS] %s",
                paste(sapply(sys.calls(), function(x) deparse(x[[1]])[1]), collapse = " > ")))
            NULL
        })

        if (!is.null(row_result)) summary_rows[[length(summary_rows) + 1L]] <- row_result
    }

    summary_table <- if (length(summary_rows) > 0) {
        do.call(rbind, summary_rows)
    } else {
        data.frame(
            sample = character(0),
            cells_raw = numeric(0),
            cells_filtered = numeric(0),
            cells_removed = numeric(0),
            pct_cells_removed = numeric(0),
            stringsAsFactors = FALSE
        )
    }

    # TOTALS row
    if (nrow(summary_table) > 0) {
        total_raw      <- sum(summary_table$cells_raw,      na.rm = TRUE)
        total_filtered <- sum(summary_table$cells_filtered, na.rm = TRUE)
        total_removed  <- sum(summary_table$cells_removed,  na.rm = TRUE)
        totals_row <- data.frame(
            sample = "TOTAL",
            cells_raw = total_raw,
            cells_filtered = total_filtered,
            cells_removed = total_removed,
            pct_cells_removed = if (total_raw > 0) round(100 * total_removed / total_raw, 1) else NA_real_,
            stringsAsFactors = FALSE
        )
        summary_table                <- rbind(summary_table, totals_row)
    }

    if (qc_pdf_open) {
        tryCatch(dev.off(), error = function(e) log_warn(paste("PDF device close:", conditionMessage(e))))
        log_info(sprintf("QC audit PDF saved: %s", qc_pdf_path))
    }

    # Print summary and save CSV
    log_info("=== QC Summary ===")
    if (nrow(summary_table) > 0) {
        print_cols <- c("sample", "cells_raw", "cells_filtered", "cells_removed", "pct_cells_removed")
        for (ln in capture.output(print(summary_table[, print_cols], row.names = FALSE))) log_info(ln)
    } else {
        log_warn("No datasets processed.")
    }

    qc_summary_csv <- file.path(QC_DIR, "qc_summary_detailed.csv")
    write.csv(summary_table, qc_summary_csv, row.names = FALSE)
    log_info(sprintf("Full summary (%d rows, %d columns) saved: %s",
                     nrow(summary_table), ncol(summary_table), qc_summary_csv))

    # ── 5.7 QC Dashboard (5-column summary schema) ───────────────────────────
    if (nrow(summary_table) == 0) {
        log_warn("Summary table empty — skipping dashboard.")
    } else {

        plot_data <- summary_table[summary_table$sample != "TOTAL", ]

        # Shorten labels for readability
        plot_data$sample_short <- r.shorten_label(plot_data$sample)

        # P1: Cell retention waterfall
        df_cells      <- melt(plot_data[, c("sample_short", "cells_raw", "cells_filtered")], id.vars = "sample_short")
        p1 <- ggplot(df_cells, aes(x = sample_short, y = value, fill = variable)) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
            scale_fill_manual(values = c("cells_raw" = "#bcbcbc", "cells_filtered" = "#2c7fb8"),
                              labels = c("Before", "After")) +
            labs(title = "Cell Retention Overview", y = "Number of Cells", x = "", fill = "Stage") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")

        # P2: Filtering intensity
        p2 <- ggplot(plot_data, aes(x = sample_short, y = pct_cells_removed, fill = sample_short)) +
            geom_col(alpha = 0.8) +
            geom_hline(yintercept = median(plot_data$pct_cells_removed),
                       linetype = "dashed", color = "black") +
            geom_hline(yintercept = 35, linetype = "dotted", color = "red", linewidth = 0.4) +
            annotate("text", x = 1, y = 36, label = "35% warning", size = 2.5, color = "red", hjust = 0) +
            labs(title = "QC Filtering Intensity", y = "% Cells Removed", x = "") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "none")

        final_layout <- (p1 | p2) +
            plot_annotation(
                title = "QC — Summary Dashboard",
                subtitle = "Cell retention and filtering intensity"
            )

        options(repr.plot.width = 14, repr.plot.height = 6)
        print(final_layout)

        ggsave(file.path(QC_DIR, "summary_qc_full_dashboard.png"),
               plot = final_layout, width = 14, height = 6, dpi = 300, bg = "white")
        log_info(sprintf("Full dashboard saved: %s", file.path(QC_DIR, "summary_qc_full_dashboard.png")))
    }

} # end RUN_QC


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — STAGE 2: Doublet Detection (scDblFinder)
# Input: QC-filtered RDS files in FILTERED_DIR
# ══════════════════════════════════════════════════════════════════════════════

if (RUN_DOUBLET) {

    log_info("──────────────────────────────────────────────────")
    log_info("STAGE 2: DOUBLET DETECTION")
    log_info(sprintf("  Input  : %s", FILTERED_DIR))
    log_info(sprintf("  Output : %s", DOUBLET_DIR))
    log_info("──────────────────────────────────────────────────")

    # ── 6.1 Discover QC-filtered RDS files ───────────────────────────────────
    filtered_files <- list.files(FILTERED_DIR, pattern = "\\.rds$", full.names = TRUE)
    if (length(filtered_files) == 0)
        stop(sprintf("No filtered RDS files found in: %s\n  Run Stage 1 first, or check --filtered_subdir.", FILTERED_DIR))

    log_info(sprintf("Found %d filtered RDS files.", length(filtered_files)))

    # ── 6.2 Per-sample DBR config (read from DCF in Section 2) ──────────────
    # DBL_SAMPLE_CFG is already populated at startup from qc_config.dcf keys:
    #   sample.<sample_name>.platform  =  10x | dropseq | indrops | ...
    #   sample.<sample_name>.dbr       =  0.08   (omit for 10x auto-estimation)
    if (length(DBL_SAMPLE_CFG) > 0) {
        log_info(sprintf("Per-sample DBR overrides loaded for %d sample(s): %s",
                         length(DBL_SAMPLE_CFG),
                         paste(names(DBL_SAMPLE_CFG), collapse = ", ")))
    } else {
        log_info(sprintf("No per-sample DBR overrides in config. Default platform: %s.",
                         DBL_DEFAULT_PLATFORM))
    }

    # ── 6.3 Open doublet audit PDF ────────────────────────────────────────────
    dbl_pdf <- file.path(DOUBLET_SUMMARY_DIR, paste0("Doublet_Audit_Report_", run_ts, ".pdf"))
    pdf(dbl_pdf, width = 16, height = 10)

    dbl_summary_stats <- list()

    # ── 6.4 Per-sample doublet detection loop ─────────────────────────────────
    for (file_path in filtered_files) {

        sample_nm <- gsub(
            "_harmonized_filtered\\.rds$|_filtered\\.rds$|\\.rds$",
            "",
            basename(file_path)
        )
        log_info(sprintf("[DOUBLET] Processing: %s", sample_nm))

        # Checkpoint / Resume logic
        singlets_path <- file.path(DOUBLET_SINGLET_DIR, paste0(sample_nm, "_singlets.rds"))
        sample_cfg <- DBL_SAMPLE_CFG[[sample_nm]]   # NULL if sample not listed
        current_dbl_params <- list(
            stage = "doublet",
            sample = sample_nm,
            dbl_batch_col = DBL_BATCH_COL,
            dbl_min_count = DBL_MIN_COUNT,
            dbl_min_feature = DBL_MIN_FEATURE,
            dbl_umap_dims = as.integer(DBL_UMAP_DIMS),
            dbl_default_platform = DBL_DEFAULT_PLATFORM,
            sample_platform = if (!is.null(sample_cfg$platform) && nzchar(sample_cfg$platform))
                                  tolower(trimws(as.character(sample_cfg$platform)))
                              else NA_character_,
            sample_dbr = if (!is.null(sample_cfg$dbr) && nzchar(sample_cfg$dbr))
                             suppressWarnings(as.numeric(sample_cfg$dbr))
                         else NA_real_
        )

        if (.should_skip_run(singlets_path, current_dbl_params, cfg$force_overwrite)) {
            log_info(sprintf("[DOUBLET SKIPPED] %s already processed (output exists).", sample_nm))
            # Just push NA stats to track it in output
            dbl_summary_stats[[sample_nm]] <- data.frame(
                Sample           = sample_nm,
                Total_Cells      = NA,
                Doublets_Found   = NA,
                Percent_Doublet  = NA,
                DBR_Source       = "SKIPPED_RESUME",
                stringsAsFactors = FALSE
            )
            next 
        }

        # Wrap entire sample processing in a tryCatch to prevent crashing pipeline
        dbl_ok <- tryCatch({
            obj <- readRDS(file_path)

        # Ghost cell filter — very low-count cells (CellBender artifacts) break
        # scDblFinder size-factor estimation
        n_before <- ncol(obj)
        obj      <- subset(obj, subset = nCount_RNA > DBL_MIN_COUNT & nFeature_RNA > DBL_MIN_FEATURE)
        log_info(sprintf("  Ghost-cell filter: %d \u2192 %d cells", n_before, ncol(obj)))

        # S4SXP fix: force meta.data to base data.frame (Seurat v5 compatibility)
        obj@meta.data <- as.data.frame(obj@meta.data)

        # Resolve DBR for this sample from DCF config (DBL_SAMPLE_CFG)
        # Fallback rule: if platform is missing or unrecognised, treat as 10x
        # and let scDblFinder auto-estimate. The pipeline never stops due to a
        # missing platform or missing DBR — it always degrades gracefully.
        dbr        <- NULL
        dbr_source <- "AUTO_10X"
        platform   <- DBL_DEFAULT_PLATFORM   # default: "10x"

        if (!is.null(sample_cfg)) {

            # Platform — default to "10x" if key absent or blank
            platform <- if (!is.null(sample_cfg$platform) && nzchar(sample_cfg$platform))
                            tolower(trimws(sample_cfg$platform))
                        else {
                            log_warn(sprintf(
                                "  No platform specified for '%s' — defaulting to 10x (AUTO)",
                                sample_nm))
                            "10x"
                        }

            dbr_user <- if (!is.null(sample_cfg$dbr) && nzchar(sample_cfg$dbr))
                            suppressWarnings(as.numeric(sample_cfg$dbr))
                        else NA

            if (platform == "10x") {
                if (!is.na(dbr_user)) {
                    if (dbr_user <= 0 || dbr_user > 0.3) {
                        log_warn(sprintf(
                            "  Invalid DBR %.4f for '%s' (must be in (0, 0.30]) — falling back to AUTO",
                            dbr_user, sample_nm))
                    } else {
                        dbr <- dbr_user; dbr_source <- "USER_OVERRIDE_10X"
                    }
                }
                # dbr stays NULL → scDblFinder auto-estimates for 10x

            } else {
                # Non-10x platform
                if (is.na(dbr_user)) {
                    log_warn(sprintf(
                        "  Platform '%s' for '%s' has no dbr — defaulting to 10x AUTO estimation",
                        platform, sample_nm))
                    platform   <- "10x"
                    dbr_source <- "AUTO_10X_FALLBACK"
                } else if (dbr_user <= 0 || dbr_user > 0.3) {
                    log_warn(sprintf(
                        "  Invalid DBR %.4f for '%s' (must be in (0, 0.30]) — falling back to 10x AUTO",
                        dbr_user, sample_nm))
                    platform   <- "10x"
                    dbr_source <- "AUTO_10X_FALLBACK"
                } else {
                    dbr <- dbr_user; dbr_source <- "USER_REQUIRED"
                }
            }

        } else {
            # Sample not listed in config at all — fully silent, use default
            dbr_source <- if (DBL_DEFAULT_PLATFORM == "10x") "AUTO_10X" else "DEFAULT_PLATFORM"
        }

        log_info(sprintf("  Platform: %s | DBR: %s [%s]",
                         platform,
                         ifelse(is.null(dbr), "AUTO", sprintf("%.4f", dbr)),
                         dbr_source))

        # Normalise + reduce (UMAP needed for visualisation)
        log_info("  Normalizing | PCA | UMAP ...")
        obj <- NormalizeData(obj,                 verbose = FALSE)
        obj <- FindVariableFeatures(obj,          verbose = FALSE)
        obj <- ScaleData(obj,                     verbose = FALSE)
        obj <- RunPCA(obj,                        verbose = FALSE)
        obj <- RunUMAP(obj, dims = DBL_UMAP_DIMS, verbose = FALSE)

        # Convert to SCE
        sce <- as.SingleCellExperiment(obj)

        # Batch-aware detection — uses DBL_BATCH_COL to prevent the ~30% doublet
        # artifact that appears when multi-sample merged objects are run as one
        has_batch <- DBL_BATCH_COL %in% colnames(colData(sce)) &&
                     length(unique(sce[[DBL_BATCH_COL]])) > 1

        if (has_batch) {
            n_sub <- length(unique(sce[[DBL_BATCH_COL]]))
            log_info(sprintf("  [BATCH MODE] %d internal samples — processing independently.", n_sub))
            sce <- if (is.null(dbr)) scDblFinder(sce, samples = DBL_BATCH_COL)
                   else              scDblFinder(sce, samples = DBL_BATCH_COL, dbr = dbr)
        } else {
            log_info("  [SINGLE MODE] No internal batching — running standard mode.")
            sce <- if (is.null(dbr)) scDblFinder(sce)
                   else              scDblFinder(sce, dbr = dbr)
        }

        # Map results back to Seurat object
        class_results <- as.character(sce$scDblFinder.class)
        score_results <- as.numeric(sce$scDblFinder.score)
        names(class_results) <- names(score_results) <- colnames(sce)
        obj$scDblFinder.class <- class_results[colnames(obj)]
        obj$scDblFinder.score <- score_results[colnames(obj)]

        n_total     <- ncol(obj)
        n_doublets  <- sum(obj$scDblFinder.class == "doublet", na.rm = TRUE)
        pct_doublet <- round((n_doublets / n_total) * 100, 2)
        log_info(sprintf("  Total: %d | Doublets: %d (%.2f%%)", n_total, n_doublets, pct_doublet))

        # Doublet audit plots (4-panel per sample)
        plot_df        <- as.data.frame(obj@meta.data)
        umap_coords    <- as.data.frame(Embeddings(obj, "umap"))
        plot_df$UMAP_1 <- umap_coords[, 1]
        plot_df$UMAP_2 <- umap_coords[, 2]
        dbl_colors     <- c("singlet" = "grey85", "doublet" = "firebrick3")

        p_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = scDblFinder.class)) +
            geom_point(size = 0.5, alpha = 0.8) +
            scale_color_manual(values = dbl_colors) +
            labs(title    = paste("Doublet Classification \u2014", sample_nm),
                 subtitle = sprintf("Doublets: %d (%.2f%%)", n_doublets, pct_doublet),
                 color = NULL) +
            theme_minimal()

        p_vln <- ggplot(plot_df, aes(scDblFinder.class, nFeature_RNA, fill = scDblFinder.class)) +
            geom_violin(alpha = 0.7, trim = FALSE) +
            scale_fill_manual(values = dbl_colors) +
            labs(title = "Feature Count Distribution", x = NULL, y = "nFeature_RNA") +
            theme_minimal() + theme(legend.position = "none")

        p_sct <- ggplot(plot_df, aes(nCount_RNA, nFeature_RNA, color = scDblFinder.class)) +
            geom_point(size = 0.5, alpha = 0.5) +
            scale_color_manual(values = dbl_colors) +
            labs(title = "Complexity Scatter", color = NULL) +
            theme_minimal()

        p_dns <- ggplot(plot_df, aes(scDblFinder.score, fill = scDblFinder.class)) +
            geom_density(alpha = 0.5) +
            scale_fill_manual(values = dbl_colors) +
            labs(title = "Score Distribution", fill = NULL) +
            theme_minimal()

        print((p_umap + p_sct) / (p_vln + p_dns))

        # Save: full object with calls + singlets-only
        full_calls_path <- file.path(DOUBLET_WITH_DOUBLETS_DIR, paste0(sample_nm, "_with_doublet_calls.rds"))
        saveRDS(obj, full_calls_path)
        obj_singlets <- subset(obj, subset = scDblFinder.class == "singlet")
        singlets_path <- file.path(DOUBLET_SINGLET_DIR, paste0(sample_nm, "_singlets.rds"))
        saveRDS(obj_singlets, singlets_path)
        .save_run_metadata(singlets_path, current_dbl_params)

        log_info(sprintf("  Saved: %s | %s", full_calls_path, singlets_path))

        dbl_summary_stats[[sample_nm]] <- data.frame(
            Sample           = sample_nm,
            Total_Cells      = n_total,
            Doublets_Found   = n_doublets,
            Percent_Doublet  = pct_doublet,
            DBR_Source       = dbr_source,
            stringsAsFactors = FALSE
        )

        rm(obj, obj_singlets, sce, plot_df, umap_coords, p_umap, p_vln, p_sct, p_dns)
        gc()
        TRUE
        
        }, error = function(e) {
            log_error(sprintf("[FAILED DOUBLET] %s: %s", sample_nm, conditionMessage(e)))
            
            # Log failure in summary stats
            dbl_summary_stats[[sample_nm]] <<- data.frame(
                Sample           = sample_nm,
                Total_Cells      = NA,
                Doublets_Found   = NA,
                Percent_Doublet  = NA,
                DBR_Source       = "ERROR",
                stringsAsFactors = FALSE
            )
            FALSE
        })
    }

    tryCatch(dev.off(), error = function(e) log_warn(paste("PDF device close:", conditionMessage(e))))
    log_info(sprintf("Doublet audit PDF saved: %s", dbl_pdf))

    # ── 6.5 Doublet summary CSV ────────────────────────────────────────────────
    if (length(dbl_summary_stats) > 0) {
        dbl_summary_df  <- do.call(rbind, dbl_summary_stats)
        dbl_summary_csv <- file.path(DOUBLET_SUMMARY_DIR, paste0("doublet_summary_", run_ts, ".csv"))
        write.csv(dbl_summary_df, dbl_summary_csv, row.names = FALSE)
        log_info(sprintf("Doublet summary CSV saved: %s", dbl_summary_csv))
        print(dbl_summary_df)
    }

} # end RUN_DOUBLET


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — Pipeline Complete
# ══════════════════════════════════════════════════════════════════════════════

log_info("══════════════════════════════════════════════════")
log_info("PIPELINE COMPLETE")
log_info(sprintf("  Stage run    : %s", RUN_STAGE))
log_info(sprintf("  Output base  : %s", BASE_OUT_DIR))
if (RUN_QC)      log_info(sprintf("  QC audit PDF : %s", file.path(QC_DIR,              "QC_Before_After_Report.pdf")))
if (RUN_QC)      log_info(sprintf("  QC dashboard : %s", file.path(QC_DIR,              "summary_qc_full_dashboard.png")))
if (RUN_QC)      log_info(sprintf("  QC summary   : %s", file.path(QC_DIR,              "qc_summary_detailed.csv")))
if (RUN_DOUBLET) log_info(sprintf("  Doublet PDF  : %s", file.path(DOUBLET_SUMMARY_DIR, paste0("Doublet_Audit_Report_", run_ts, ".pdf"))))
if (RUN_DOUBLET) log_info(sprintf("  Doublet CSV  : %s", file.path(DOUBLET_SUMMARY_DIR, paste0("doublet_summary_",      run_ts, ".csv"))))
log_info(sprintf("  Log file     : %s", log_file_path))
log_info("══════════════════════════════════════════════════")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — Integrated QC + Doublet End-to-End Summary
#
# Runs when both stages have been executed (stage = "all") or when outputs
# from both stages are already present on disk.
#
# Outputs:
#   integrated_qc_doublet_summary.csv  — one row per sample with:
#       cells_raw | cells_post_qc | cells_post_doublet |
#       pct_removed_by_qc | pct_removed_by_doublet | pct_removed_total
#   integrated_summary_dashboard.png   — 9-panel figure covering:
#       (1) 3-stage waterfall bar chart
#       (2) stacked % loss breakdown
#       (3) cumulative funnel (line)
#       (4) log1p(nCount_RNA) density overlay: raw → post-QC → post-doublet
#       (5) log1p(nFeature_RNA) density overlay
#       (6) % mito overlay (raw vs post-QC)
#       (7) doublet score distribution per sample (ridge plot)
#       (8) doublet % vs QC removal % scatter
#       (9) total cells lost per stage (pie / donut)
# ══════════════════════════════════════════════════════════════════════════════

log_info("──────────────────────────────────────────────────")
log_info("SECTION 8: INTEGRATED QC + DOUBLET SUMMARY")
log_info("──────────────────────────────────────────────────")

# ── 8.0 Guard: need at least one stage's outputs to proceed ──────────────────
qc_csv_path  <- file.path(QC_DIR, "qc_summary_detailed.csv")
dbl_csv_glob <- list.files(DOUBLET_SUMMARY_DIR, pattern = "^doublet_summary_.*\\.csv$", full.names = TRUE)
dbl_csv_path <- if (length(dbl_csv_glob) > 0) tail(sort(dbl_csv_glob), 1) else NULL

has_qc_summary  <- file.exists(qc_csv_path)
has_dbl_summary <- !is.null(dbl_csv_path) && file.exists(dbl_csv_path)

if (!has_qc_summary && !has_dbl_summary) {
    log_warn("Section 8: Neither QC nor doublet summary CSVs found on disk — skipping integrated summary.")
} else {

    # ── 8.1 Load per-stage summaries ─────────────────────────────────────────

    # QC summary
    if (has_qc_summary) {
        qc_sum <- read.csv(qc_csv_path, stringsAsFactors = FALSE)
        qc_sum <- qc_sum[qc_sum$sample != "TOTAL", ]   # drop totals row
    } else {
        qc_sum <- data.frame(sample = character(0), cells_raw = integer(0), cells_filtered = integer(0))
        log_warn("Section 8: QC summary CSV missing — cells_raw and cells_post_qc will be NA.")
    }

    # Doublet summary
    if (has_dbl_summary) {
        dbl_sum <- read.csv(dbl_csv_path, stringsAsFactors = FALSE)
        dbl_sum <- dbl_sum[dbl_sum$DBR_Source != "ERROR", ]   # drop failed samples
    } else {
        dbl_sum <- data.frame(Sample = character(0), Total_Cells = integer(0), Doublets_Found = integer(0))
        log_warn("Section 8: Doublet summary CSV missing — doublet columns will be NA.")
    }

    # ── 8.2 Build unified table ───────────────────────────────────────────────
    # QC uses study names ending with "_harmonized" while doublet summaries
    # use non-harmonized names. Merge on a normalized key to avoid NA splits.
    .normalize_study_key <- function(x) {
        x <- as.character(x)
        x <- sub("_harmonized$", "", x)
        x
    }
    .pick_study_label <- function(x) {
        x <- as.character(x)
        harmonized <- x[grepl("_harmonized$", x)]
        if (length(harmonized) > 0) harmonized[1] else x[1]
    }

    if (has_qc_summary && nrow(qc_sum) > 0) {
        qc_tbl <- qc_sum[, c("sample", "cells_raw", "cells_filtered")]
        qc_tbl$study_key <- .normalize_study_key(qc_tbl$sample)
        qc_tbl$cells_post_qc <- qc_tbl$cells_filtered

        qc_labels <- aggregate(sample ~ study_key, data = qc_tbl, FUN = .pick_study_label)
        qc_counts <- aggregate(cbind(cells_raw, cells_post_qc) ~ study_key, data = qc_tbl, FUN = sum, na.rm = TRUE)
        qc_agg <- merge(qc_labels, qc_counts, by = "study_key", all = TRUE)
    } else {
        qc_agg <- data.frame(study_key = character(0), sample = character(0),
                             cells_raw = numeric(0), cells_post_qc = numeric(0),
                             stringsAsFactors = FALSE)
    }

    if (has_dbl_summary && nrow(dbl_sum) > 0) {
        dbl_tbl <- dbl_sum[, c("Sample", "Total_Cells", "Doublets_Found", "Percent_Doublet")]
        dbl_tbl$study_key <- .normalize_study_key(dbl_tbl$Sample)
        dbl_tbl$cells_post_qc_dbl <- dbl_tbl$Total_Cells
        dbl_tbl$doublets_found <- dbl_tbl$Doublets_Found
        dbl_tbl$pct_doublet <- dbl_tbl$Percent_Doublet

        dbl_counts <- aggregate(cbind(cells_post_qc_dbl, doublets_found) ~ study_key,
                                data = dbl_tbl, FUN = sum, na.rm = TRUE)
        dbl_pct <- aggregate(pct_doublet ~ study_key, data = dbl_tbl,
                             FUN = function(v) {
                                 if (all(is.na(v))) NA_real_ else round(mean(v, na.rm = TRUE), 2)
                             })
        dbl_agg <- merge(dbl_counts, dbl_pct, by = "study_key", all = TRUE)
    } else {
        dbl_agg <- data.frame(study_key = character(0), cells_post_qc_dbl = numeric(0),
                              doublets_found = numeric(0), pct_doublet = numeric(0),
                              stringsAsFactors = FALSE)
    }

    int_df <- merge(qc_agg, dbl_agg, by = "study_key", all = TRUE)
    if (!"sample" %in% colnames(int_df)) int_df$sample <- NA_character_
    int_df$sample <- ifelse(!is.na(int_df$sample) & nzchar(int_df$sample),
                            int_df$sample, int_df$study_key)

    # Prefer QC post-filtering count when available.
    int_df$cells_post_qc <- ifelse(
        !is.na(int_df$cells_post_qc), int_df$cells_post_qc, int_df$cells_post_qc_dbl)
    int_df$cells_post_doublet <- ifelse(
        !is.na(int_df$cells_post_qc) & !is.na(int_df$doublets_found),
        int_df$cells_post_qc - int_df$doublets_found,
        NA_real_)
    int_df$cells_post_qc_dbl <- NULL
    int_df$study_key <- NULL

    if (nrow(int_df) > 0) int_df <- int_df[order(int_df$sample), , drop = FALSE]

    # Derived percentages
    int_df$cells_removed_by_qc      <- int_df$cells_raw      - int_df$cells_post_qc
    int_df$cells_removed_by_doublet <- int_df$cells_post_qc  - int_df$cells_post_doublet
    int_df$cells_removed_total       <- int_df$cells_raw      - int_df$cells_post_doublet

    int_df$pct_removed_by_qc      <- ifelse(
        !is.na(int_df$cells_raw) & int_df$cells_raw > 0,
        round(int_df$cells_removed_by_qc / int_df$cells_raw * 100, 2),
        NA_real_)
    int_df$pct_removed_by_doublet <- ifelse(
        !is.na(int_df$cells_raw) & int_df$cells_raw > 0,
        round(int_df$cells_removed_by_doublet / int_df$cells_raw * 100, 2),
        NA_real_)
    int_df$pct_removed_total       <- ifelse(
        !is.na(int_df$cells_raw) & int_df$cells_raw > 0,
        round(int_df$cells_removed_total / int_df$cells_raw * 100, 2),
        NA_real_)
    int_df$pct_retained_final      <- 100 - int_df$pct_removed_total

    # Shorten sample names for plots
    int_df$sample_short <- r.shorten_label(int_df$sample)

    # Add TOTAL row
    mean_pct_doublet <- if (nrow(int_df) > 0 && !all(is.na(int_df$pct_doublet)))
        round(mean(int_df$pct_doublet, na.rm = TRUE), 2) else NA_real_

    total_row <- data.frame(
        sample                   = "TOTAL",
        cells_raw                = sum(int_df$cells_raw,                na.rm = TRUE),
        cells_post_qc            = sum(int_df$cells_post_qc,           na.rm = TRUE),
        cells_post_doublet       = sum(int_df$cells_post_doublet,      na.rm = TRUE),
        doublets_found           = sum(int_df$doublets_found,          na.rm = TRUE),
        pct_doublet              = mean_pct_doublet,
        cells_removed_by_qc      = sum(int_df$cells_removed_by_qc,     na.rm = TRUE),
        cells_removed_by_doublet = sum(int_df$cells_removed_by_doublet,na.rm = TRUE),
        cells_removed_total      = sum(int_df$cells_removed_total,     na.rm = TRUE),
        pct_removed_by_qc        = NA_real_,
        pct_removed_by_doublet   = NA_real_,
        pct_removed_total        = NA_real_,
        pct_retained_final       = NA_real_,
        sample_short             = "TOTAL",
        stringsAsFactors = FALSE
    )
    total_row$pct_removed_by_qc      <- round(total_row$cells_removed_by_qc      / total_row$cells_raw * 100, 2)
    total_row$pct_removed_by_doublet <- round(total_row$cells_removed_by_doublet / total_row$cells_raw * 100, 2)
    total_row$pct_removed_total       <- round(total_row$cells_removed_total      / total_row$cells_raw * 100, 2)
    total_row$pct_retained_final      <- 100 - total_row$pct_removed_total

    int_df_full <- rbind(int_df, total_row)

    # ── 8.3 Save integrated CSV ───────────────────────────────────────────────
    INTEGRATED_SUMMARY_DIR <- file.path(BASE_OUT_DIR, "integrated_summary")
    dir.create(INTEGRATED_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

    # Requested integrated schema (study-level)
    int_export <- int_df_full[, c(
        "sample", "cells_raw", "cells_post_qc", "pct_removed_by_qc",
        "cells_post_doublet", "pct_doublet", "cells_post_doublet"
    )]
    colnames(int_export) <- c(
        "Study", "cells_raw", "cells_after_filtering", "pct_filtered",
        "cells_after_doublet_filtering", "pct_doublet", "cells_after_QC"
    )

    int_csv_path <- file.path(INTEGRATED_SUMMARY_DIR, "integrated_qc_doublet_summary.csv")
    write.csv(int_export, int_csv_path, row.names = FALSE)
    log_info(sprintf("Integrated summary CSV saved: %s", int_csv_path))

    # Print to log
    log_info("=== Integrated QC + Doublet Summary ===")
        for (ln in capture.output(print(int_export, row.names = FALSE))) log_info(ln)

    # ── 8.4 Build integrated dashboard plots ─────────────────────────────────
    # Work only with per-sample rows (exclude TOTAL) for per-sample plots
    pd <- int_df[complete.cases(int_df[, c("cells_raw", "cells_post_qc")]), ]

    if (nrow(pd) == 0) {
        log_warn("Section 8: No complete sample rows — skipping dashboard plots.")
    } else {

        # ── Panel 1: 3-stage waterfall bar chart ──────────────────────────────
        df_waterfall <- melt(
            pd[, c("sample_short", "cells_raw", "cells_post_qc", "cells_post_doublet")],
            id.vars = "sample_short"
        )
        df_waterfall$variable <- factor(df_waterfall$variable,
            levels = c("cells_raw", "cells_post_qc", "cells_post_doublet"),
            labels = c("Raw", "Post-QC", "Post-Doublet"))

        p8_1 <- ggplot(df_waterfall, aes(x = sample_short, y = value, fill = variable)) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.85, width = 0.7) +
            scale_fill_manual(values = c(
                "Raw"          = "#bdbdbd",
                "Post-QC"      = "#2c7fb8",
                "Post-Doublet" = "#1a9641")) +
            scale_y_continuous(labels = scales::comma) +
            labs(title    = "Cell Counts Across All Three Stages",
                 subtitle = "Raw → Post-QC → Post-Doublet removal",
                 x = "", y = "Number of Cells", fill = "Stage") +
            theme_minimal(base_size = 10) +
            theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
                  legend.position = "top",
                  plot.title   = element_text(face = "bold"))

        # ── Panel 2: Stacked % loss per stage ─────────────────────────────────
        df_pct_loss <- melt(
            pd[, c("sample_short", "pct_removed_by_qc", "pct_removed_by_doublet")],
            id.vars = "sample_short"
        )
        df_pct_loss$variable <- factor(df_pct_loss$variable,
            levels = c("pct_removed_by_qc", "pct_removed_by_doublet"),
            labels = c("QC Removal", "Doublet Removal"))

        p8_2 <- ggplot(df_pct_loss, aes(x = sample_short, y = value, fill = variable)) +
            geom_bar(stat = "identity", position = "stack", alpha = 0.85) +
            geom_hline(yintercept = median(pd$pct_removed_total, na.rm = TRUE),
                       linetype = "dashed", color = "black", linewidth = 0.5) +
            geom_hline(yintercept = 50, linetype = "dotted", color = "red", linewidth = 0.4) +
            annotate("text", x = 0.8, y = 51.5, label = "50% loss", size = 2.5, color = "red", hjust = 0) +
            scale_fill_manual(values = c("QC Removal" = "#fc8d59", "Doublet Removal" = "#d73027")) +
            labs(title    = "Stacked % Cell Loss by Stage",
                 subtitle = "Dashed = cohort median total loss",
                 x = "", y = "% Cells Removed", fill = "Stage") +
            theme_minimal(base_size = 10) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  legend.position = "top",
                  plot.title  = element_text(face = "bold"))

        # ── Panel 3: Cumulative funnel (connected lines) ───────────────────────
        df_funnel <- pd[, c("sample_short", "cells_raw", "cells_post_qc", "cells_post_doublet")]
        df_funnel_long <- melt(df_funnel, id.vars = "sample_short")
        df_funnel_long$stage <- factor(df_funnel_long$variable,
            levels = c("cells_raw", "cells_post_qc", "cells_post_doublet"),
            labels = c("Raw", "Post-QC", "Post-Doublet"))

        p8_3 <- ggplot(df_funnel_long, aes(x = stage, y = value,
                                            group = sample_short, color = sample_short)) +
            geom_line(linewidth = 0.8, alpha = 0.8) +
            geom_point(size = 2.5, alpha = 0.9) +
            scale_y_continuous(labels = scales::comma) +
            scale_color_viridis_d(option = "turbo", end = 0.92) +
            labs(title    = "Per-Sample Cell Count Funnel",
                 subtitle = "Each line = one sample tracked across stages",
                 x = "Stage", y = "Cells", color = "") +
            theme_minimal(base_size = 10) +
            theme(legend.position = "right",
                  legend.text     = element_text(size = 7),
                  plot.title      = element_text(face = "bold"))

        # ── Panel 4 & 5: log1p distributions — load RDS files for density ─────
        # We collect per-sample data from:
        #   raw RDS → rds_file_map (if Stage 1 ran)
        #   QC-filtered RDS → FILTERED_DIR
        #   Doublet-filtered (singlets) → DOUBLET_SINGLET_DIR
        density_list <- list()

        for (samp in pd$sample) {
            samp_short <- pd$sample_short[pd$sample == samp]

            # Raw
            raw_path <- if (exists("rds_file_map")) rds_file_map[[samp]] else NULL
            # QC-filtered
            qc_path  <- file.path(FILTERED_DIR, paste0(samp, "_filtered.rds"))
            # Doublet-filtered singlets (try two common naming patterns)
            dbl_path <- file.path(DOUBLET_SINGLET_DIR, paste0(samp, "_singlets.rds"))
            if (!file.exists(dbl_path))
                dbl_path <- file.path(DOUBLET_SINGLET_DIR, paste0(samp, "_harmonized_filtered_singlets.rds"))
            if (!file.exists(dbl_path))
                dbl_path <- file.path(DOUBLET_DIR, paste0(samp, "_singlets.rds"))
            if (!file.exists(dbl_path))
                dbl_path <- file.path(DOUBLET_DIR, paste0(samp, "_harmonized_filtered_singlets.rds"))

            for (stage_info in list(
                    list(path = raw_path,  label = "Raw"),
                    list(path = qc_path,   label = "Post-QC"),
                    list(path = dbl_path,  label = "Post-Doublet")
                )) {
                p <- stage_info$path; lbl <- stage_info$label
                if (is.null(p) || !file.exists(p)) next
                obj_tmp <- tryCatch(readRDS(p), error = function(e) NULL)
                if (is.null(obj_tmp)) next

                nc  <- as.numeric(obj_tmp$nCount_RNA)
                nf  <- as.numeric(obj_tmp$nFeature_RNA)
                pmt <- if ("percent.mito" %in% colnames(obj_tmp@meta.data))
                           as.numeric(obj_tmp$percent.mito) else rep(NA, ncol(obj_tmp))

                # Subsample large datasets to keep memory manageable (≤20 k cells)
                n_samp <- length(nc)
                idx    <- if (n_samp > 20000) sample(n_samp, 20000) else seq_len(n_samp)

                density_list[[paste(samp_short, lbl, sep = "||")]] <- data.frame(
                    sample_short  = samp_short,
                    stage         = lbl,
                    log1p_nCount  = log1p(nc[idx]),
                    log1p_nFeat   = log1p(nf[idx]),
                    percent_mito  = pmt[idx]
                )
                rm(obj_tmp); gc()
            }
        }

        has_density <- length(density_list) > 0
        if (has_density) {
            dens_df       <- do.call(rbind, density_list)
            dens_df$stage <- factor(dens_df$stage, levels = c("Raw", "Post-QC", "Post-Doublet"))
            stage_cols    <- c("Raw" = "#bdbdbd", "Post-QC" = "#2c7fb8", "Post-Doublet" = "#1a9641")

            # Panel 4: log1p(nCount_RNA) — global density overlay
            p8_4 <- ggplot(dens_df, aes(x = log1p_nCount, fill = stage, color = stage)) +
                geom_density(alpha = 0.30, linewidth = 0.6) +
                scale_fill_manual(values  = stage_cols) +
                scale_color_manual(values = stage_cols) +
                facet_wrap(~sample_short, scales = "free_y", ncol = 4) +
                labs(title    = "log1p(nCount_RNA) Distribution: Raw vs Post-QC vs Post-Doublet",
                     subtitle = "Overlapping density per dataset",
                     x = "log1p(nCount_RNA)", y = "Density",
                     fill = "Stage", color = "Stage") +
                theme_minimal(base_size = 9) +
                theme(legend.position = "top",
                      strip.text      = element_text(size = 7),
                      plot.title      = element_text(face = "bold"))

            # Panel 5: log1p(nFeature_RNA) — global density overlay
            p8_5 <- ggplot(dens_df, aes(x = log1p_nFeat, fill = stage, color = stage)) +
                geom_density(alpha = 0.30, linewidth = 0.6) +
                scale_fill_manual(values  = stage_cols) +
                scale_color_manual(values = stage_cols) +
                facet_wrap(~sample_short, scales = "free_y", ncol = 4) +
                labs(title    = "log1p(nFeature_RNA) Distribution: Raw vs Post-QC vs Post-Doublet",
                     subtitle = "Overlapping density per dataset",
                     x = "log1p(nFeature_RNA)", y = "Density",
                     fill = "Stage", color = "Stage") +
                theme_minimal(base_size = 9) +
                theme(legend.position = "top",
                      strip.text      = element_text(size = 7),
                      plot.title      = element_text(face = "bold"))

            # Panel 6: % mito overlay (Raw vs Post-QC only — doublet step doesn't change mito)
            dens_mito <- dens_df[dens_df$stage %in% c("Raw", "Post-QC") &
                                 !is.na(dens_df$percent_mito), ]
            p8_6 <- if (nrow(dens_mito) > 0) {
                ggplot(dens_mito, aes(x = percent_mito, fill = stage, color = stage)) +
                    geom_density(alpha = 0.35, linewidth = 0.6) +
                    scale_fill_manual(values  = stage_cols[c("Raw", "Post-QC")]) +
                    scale_color_manual(values = stage_cols[c("Raw", "Post-QC")]) +
                    facet_wrap(~sample_short, scales = "free_y", ncol = 4) +
                    labs(title    = "Mitochondrial % Distribution: Raw vs Post-QC",
                         subtitle = "Shift shows efficiency of mito-based filtering",
                         x = "% Mitochondrial Reads", y = "Density",
                         fill = "Stage", color = "Stage") +
                    theme_minimal(base_size = 9) +
                    theme(legend.position = "top",
                          strip.text      = element_text(size = 7),
                          plot.title      = element_text(face = "bold"))
            } else {
                ggplot() + annotate("text", x = 0.5, y = 0.5,
                    label = "% mito data unavailable", size = 5) +
                    theme_void() + ggtitle("Mitochondrial % — no data")
            }

        } else {
            # Placeholder panels when RDS not accessible
            make_na_panel <- function(title) {
                ggplot() +
                    annotate("text", x = 0.5, y = 0.5, label = paste(title, "\n(RDS files needed)"),
                             size = 4, color = "grey60") +
                    theme_void() + ggtitle(title)
            }
            p8_4 <- make_na_panel("log1p(nCount_RNA) density overlay")
            p8_5 <- make_na_panel("log1p(nFeature_RNA) density overlay")
            p8_6 <- make_na_panel("Mitochondrial % overlay")
        }

        # ── Panel 7: Doublet score ridge (if singlets RDS accessible) ─────────
        dbl_score_list <- list()
        for (samp in pd$sample) {
            samp_short <- pd$sample_short[pd$sample == samp]
            dbl_full   <- file.path(DOUBLET_WITH_DOUBLETS_DIR, paste0(samp, "_with_doublet_calls.rds"))
            if (!file.exists(dbl_full))
                dbl_full <- file.path(DOUBLET_DIR, paste0(samp, "_with_doublet_calls.rds"))
            if (!file.exists(dbl_full)) next
            obj_tmp <- tryCatch(readRDS(dbl_full), error = function(e) NULL)
            if (is.null(obj_tmp)) next
            if (!"scDblFinder.score" %in% colnames(obj_tmp@meta.data)) { rm(obj_tmp); next }
            n_tmp <- ncol(obj_tmp)
            idx   <- if (n_tmp > 10000) sample(n_tmp, 10000) else seq_len(n_tmp)
            dbl_score_list[[samp_short]] <- data.frame(
                sample_short         = samp_short,
                scDblFinder.score    = as.numeric(obj_tmp$scDblFinder.score)[idx],
                scDblFinder.class    = as.character(obj_tmp$scDblFinder.class)[idx]
            )
            rm(obj_tmp); gc()
        }

        p8_7 <- if (length(dbl_score_list) > 0) {
            dbl_score_df <- do.call(rbind, dbl_score_list)
            ggplot(dbl_score_df, aes(x = scDblFinder.score, y = sample_short,
                                     fill = sample_short)) +
                geom_density_ridges(alpha = 0.70, scale = 1.4,
                                    color = "white", linewidth = 0.35) +
                geom_vline(xintercept = 0.5, linetype = "dashed",
                           color = "firebrick3", linewidth = 0.6) +
                annotate("text", x = 0.52, y = 1,
                         label = "0.5 threshold", color = "firebrick3",
                         size = 2.8, hjust = 0) +
                scale_x_continuous(limits = c(0, 1)) +
                scale_fill_viridis_d(option = "mako", end = 0.88) +
                labs(title    = "scDblFinder Score Distribution per Sample",
                     subtitle = "Dashed red = 0.5 decision boundary",
                     x = "Doublet Score", y = "") +
                theme_minimal(base_size = 9) +
                theme(legend.position = "none",
                      axis.text.y     = element_text(size = 8),
                      plot.title      = element_text(face = "bold"))
        } else {
            ggplot() +
                annotate("text", x = 0.5, y = 0.5,
                         label = "Doublet scores unavailable\n(needs *_with_doublet_calls.rds)",
                         size = 4, color = "grey60") +
                theme_void() + ggtitle("scDblFinder Score Distribution")
        }

        # ── Panel 8: Doublet % vs QC removal % scatter ────────────────────────
        pd_scatter <- pd[!is.na(pd$pct_doublet) & !is.na(pd$pct_removed_by_qc), ]
        p8_8 <- if (nrow(pd_scatter) > 1) {
            ggplot(pd_scatter, aes(x = pct_removed_by_qc, y = pct_doublet,
                                   label = sample_short, color = pct_removed_total)) +
                geom_point(size = 4, alpha = 0.85) +
                ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15,
                                         show.legend = FALSE) +
                scale_color_viridis_c(option = "inferno", direction = -1,
                                      name = "Total\n% removed") +
                geom_smooth(method = "lm", se = TRUE, linewidth = 0.7,
                            color = "grey40", fill = "grey85") +
                labs(title    = "Doublet Rate vs QC Filtering Stringency",
                     subtitle = "Point colour = total % cells removed; line = OLS trend",
                     x = "% Removed by QC", y = "% Doublets Detected") +
                theme_minimal(base_size = 10) +
                theme(plot.title = element_text(face = "bold"))
        } else {
            ggplot(pd, aes(x = pct_removed_by_qc, y = pct_doublet)) +
                geom_point(size = 3, alpha = 0.7) +
                labs(title = "Doublet Rate vs QC Stringency",
                     x = "% Removed by QC", y = "% Doublets") +
                theme_minimal(base_size = 10)
        }

        # ── Panel 9: Donut / pie — total cells lost per stage across cohort ───
        total_qc_loss  <- sum(pd$cells_removed_by_qc,      na.rm = TRUE)
        total_dbl_loss <- sum(pd$cells_removed_by_doublet, na.rm = TRUE)
        total_retained <- sum(pd$cells_post_doublet,       na.rm = TRUE)
        donut_df <- data.frame(
            category = c("Retained", "Removed by QC", "Removed by Doublet"),
            count    = c(total_retained, total_qc_loss, total_dbl_loss)
        )
        donut_df$fraction <- donut_df$count / sum(donut_df$count)
        donut_df$label    <- sprintf("%s\n(%.1f%%)", scales::comma(donut_df$count),
                                     donut_df$fraction * 100)
        donut_df$category <- factor(donut_df$category,
            levels = c("Retained", "Removed by QC", "Removed by Doublet"))

        p8_9 <- ggplot(donut_df, aes(x = 2, y = fraction, fill = category)) +
            geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.6) +
            coord_polar(theta = "y", start = 0) +
            xlim(0.5, 2.5) +     # donut hole
            scale_fill_manual(values = c(
                "Retained"           = "#1a9641",
                "Removed by QC"      = "#fc8d59",
                "Removed by Doublet" = "#d73027")) +
            geom_text(aes(label = label),
                      position = position_stack(vjust = 0.5),
                      size = 3, color = "white", fontface = "bold") +
            labs(title    = "Cohort-Level Cell Budget",
                 subtitle = sprintf("Total input: %s cells", scales::comma(sum(donut_df$count))),
                 fill = "") +
            theme_void(base_size = 10) +
            theme(legend.position = "right",
                  plot.title      = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle   = element_text(hjust = 0.5))

        # ── 8.5 Assemble 9-panel dashboard ───────────────────────────────────
        # Layout:
        #   Row 1 (equal): P1 (waterfall) | P2 (stacked %) | P3 (funnel)
        #   Row 2 (equal): P4 (nCount density facets, spans 2) | P7 (ridge, 1)
        #   Row 3 (equal): P5 (nFeat density facets, spans 2)  | P8 (scatter, 1)
        #   Row 4 (equal): P6 (mito facets, spans 2)           | P9 (donut, 1)

        int_dashboard <- (p8_1 | p8_2 | p8_3) /
                         (p8_4 | p8_7) /
                         (p8_5 | p8_8) /
                         (p8_6 | p8_9) +
            plot_layout(heights = c(1, 2.2, 2.2, 2.2)) +
            plot_annotation(
                title    = "Integrated QC + Doublet Removal \u2014 End-to-End Summary",
                subtitle = sprintf(
                    "Pipeline run: %s  |  %d samples  |  Raw: %s  \u2192  Final: %s cells (%.1f%% retained)",
                    run_ts,
                    nrow(pd),
                    scales::comma(sum(pd$cells_raw,          na.rm = TRUE)),
                    scales::comma(sum(pd$cells_post_doublet, na.rm = TRUE)),
                    100 - mean(pd$pct_removed_total, na.rm = TRUE)
                ),
                theme = theme(
                    plot.title    = element_text(size = 16, face = "bold",  hjust = 0.5),
                    plot.subtitle = element_text(size = 10,                  hjust = 0.5)
                )
            )

        # ── 8.6 Save dashboard PNG ────────────────────────────────────────────
        int_png_path <- file.path(INTEGRATED_SUMMARY_DIR, "integrated_summary_dashboard.png")
        ggsave(int_png_path, plot = int_dashboard,
               width = 22, height = 40, dpi = 300, bg = "white")
        log_info(sprintf("Integrated dashboard saved: %s", int_png_path))

        # ── 8.7 Also embed dashboard in a PDF ─────────────────────────────────
        int_pdf_path <- file.path(INTEGRATED_SUMMARY_DIR, "integrated_summary_dashboard.pdf")
        pdf(int_pdf_path, width = 22, height = 40)
        print(int_dashboard)
        tryCatch(dev.off(), error = function(e) log_warn(paste("Integrated PDF close:", conditionMessage(e))))
        log_info(sprintf("Integrated dashboard PDF saved: %s", int_pdf_path))

    } # end nrow(pd) > 0 guard

    log_info("──────────────────────────────────────────────────")
    log_info("SECTION 8: INTEGRATED SUMMARY COMPLETE")
    log_info(sprintf("  CSV      : %s", int_csv_path))
    if (exists("int_png_path")) log_info(sprintf("  PNG      : %s", int_png_path))
    if (exists("int_pdf_path")) log_info(sprintf("  PDF      : %s", int_pdf_path))
    log_info("──────────────────────────────────────────────────")

} # end has_qc_summary || has_dbl_summary guard


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Final Pipeline Footer
# ══════════════════════════════════════════════════════════════════════════════

log_info("══════════════════════════════════════════════════")
log_info("ALL SECTIONS COMPLETE")
log_info(sprintf("  Stage run         : %s", RUN_STAGE))
log_info(sprintf("  Output base       : %s", BASE_OUT_DIR))
if (RUN_QC)      log_info(sprintf("  QC audit PDF      : %s", file.path(QC_DIR,              "QC_Before_After_Report.pdf")))
if (RUN_QC)      log_info(sprintf("  QC dashboard      : %s", file.path(QC_DIR,              "summary_qc_full_dashboard.png")))
if (RUN_QC)      log_info(sprintf("  QC summary CSV    : %s", file.path(QC_DIR,              "qc_summary_detailed.csv")))
if (RUN_DOUBLET) log_info(sprintf("  Doublet PDF       : %s", file.path(DOUBLET_SUMMARY_DIR, paste0("Doublet_Audit_Report_", run_ts, ".pdf"))))
if (RUN_DOUBLET) log_info(sprintf("  Doublet CSV       : %s", file.path(DOUBLET_SUMMARY_DIR, paste0("doublet_summary_",      run_ts, ".csv"))))
if (RUN_DOUBLET) log_info(sprintf("  Singlets dir      : %s", DOUBLET_SINGLET_DIR))
if (RUN_DOUBLET) log_info(sprintf("  With doublets dir : %s", DOUBLET_WITH_DOUBLETS_DIR))
log_info(sprintf("  Integrated CSV    : %s", file.path(BASE_OUT_DIR, "integrated_summary", "integrated_qc_doublet_summary.csv")))
log_info(sprintf("  Integrated PNG    : %s", file.path(BASE_OUT_DIR, "integrated_summary", "integrated_summary_dashboard.png")))
log_info(sprintf("  Integrated PDF    : %s", file.path(BASE_OUT_DIR, "integrated_summary", "integrated_summary_dashboard.pdf")))
log_info(sprintf("  Log file          : %s", log_file_path))
log_info("══════════════════════════════════════════════════")

} # end script execution guard