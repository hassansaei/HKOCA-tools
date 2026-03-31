# =============================================================================
# Single-Cell QC Pipeline — Full Integration
# Doublet Detection (scDblFinder)  →  QC Filtering
#
# Description:
#   Stage 1 (Doublet): Loads raw Seurat RDS files, runs scDblFinder (batch-aware,
#     with optional per-sample DBR via DCF dotted keys: sample.<name>.platform /
#     sample.<name>.dbr), saves singlet-only RDS + audit PDF + summary CSV.
#   Stage 2 (QC): Reads doublet-filtered RDS, applies config-driven thresholds
#     (nFeature, nCount, percent.mito), saves filtered RDS + audit PDF +
#     full 7-panel dashboard (bar, intensity, inflation, line plots, ridge plots).
#   Stage 3 (Summary): Integrated Doublet+QC end-to-end summary CSV + plots.
#   Both stages share a single config file, a single logger, and the same CLI.
#
# Config file (DCF format — qc_config.dcf):
#   rds_dir             = /path/to/raw_rds
#   qc_decisions_table  = CSV content embedded directly in this config (header + rows)
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
#   dbl_default_chemistry  = v3            (v2/v3 use same table; v4 has higher multiplet rates)
#   sample.<name>.platform = 10x | dropseq | indrops | ...
#   sample.<name>.chemistry = v3 | v4      (10x chemistry version — controls multiplet table)
#   sample.<name>.dbr      = 0.08   (omit for 10x auto-estimation)
#
# Usage:
#   Rscript QC_scdbl_Combined.R [--config PATH] [--stage all|qc|doublet] [OPTIONS]
#
# CLI options:
#   --config PATH           Config DCF file (default: qc_config.dcf next to script)
#   --stage STR             all | qc | doublet  (default: all)
#   --rds_dir PATH          Override raw RDS input directory
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
        "scales", "reshape2", "ggridges", "ggrepel", "SingleCellExperiment", 
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
        library(ggrepel)
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

.as_num <- function(x, default = NA_real_) {
    if (is.null(x) || length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") return(default)
    v <- suppressWarnings(as.numeric(x))
    if (is.na(v)) default else v
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

    # Tolerant DCF parser with multiline continuation support.
    # This accepts embedded CSV blocks (qc_decisions_table) where each row is
    # placed on its own subsequent line.
    lines <- readLines(config_path, warn = FALSE)
    if (length(lines) == 0)
        stop(sprintf("Config file is empty: %s", config_path))

    cfg <- list()
    current_key <- NULL
    current_val <- ""

    .flush_current <- function() {
        if (!is.null(current_key)) {
            cfg[[current_key]] <<- current_val
        }
    }

    for (ln in lines) {
        if (grepl("^\\s*#", ln) || !nzchar(trimws(ln))) next

        key_match <- regexec("^\\s*([^:]+):\\s*(.*)$", ln)
        key_parts <- regmatches(ln, key_match)[[1]]

        if (length(key_parts) == 3) {
            .flush_current()
            current_key <- trimws(key_parts[2])
            current_val <- key_parts[3]
        } else {
            if (is.null(current_key)) {
                stop(sprintf("Incorrect DCF format near line: %s", ln))
            }
            current_val <- paste(current_val, sub("\\s+$", "", ln), sep = "\n")
        }
    }

    .flush_current()

    if (length(cfg) == 0)
        stop(sprintf("Config file contains no valid key:value entries: %s", config_path))

    cfg
}

# ── 1.6 Null-coalesce operator ────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── 1.7 QC decisions loader (config-embedded table) ──────────────────────────
.load_qc_decisions <- function(cfg) {
    inline_table <- cfg$qc_decisions_table %||% NULL

    if (is.null(inline_table) || !nzchar(trimws(inline_table))) {
        stop("Config key 'qc_decisions_table' must be set in the config file.")
    }

    lines <- unlist(strsplit(gsub("\r\n?", "\n", inline_table), "\n", fixed = TRUE))
    lines <- lines[nzchar(trimws(lines))]

    if (length(lines) < 2) {
        stop("Config key 'qc_decisions_table' must contain a CSV header plus at least one data row.")
    }

    qc_decisions <- read.csv(
        text = paste(lines, collapse = "\n"),
        stringsAsFactors = FALSE,
        colClasses = "character",
        strip.white = TRUE
    )

    if (!("Dataset_Name" %in% colnames(qc_decisions))) {
        stop("Config key 'qc_decisions_table' must contain a 'Dataset_Name' column.")
    }

    required_cols <- c("Dataset_Name", "Lower_Feature_Method", "Upper_Feature_Method",
                       "Lower_Count_Method", "Upper_Count_Method", "Upper_Mito_Method")
    missing_cols  <- setdiff(required_cols, colnames(qc_decisions))
    if (length(missing_cols) > 0) {
        stop(sprintf("qc_decisions_table is missing required column(s): %s",
                     paste(missing_cols, collapse = ", ")))
    }

    list(data = qc_decisions, source = "config:qc_decisions_table")
}

# ── 1.8 Logger ────────────────────────────────────────────────────────────────
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

# ── 1.9 RDS file discovery ────────────────────────────────────────────────────
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

# Normalise sample names so config keys, QC outputs, and doublet inputs use the
# same lookup key regardless of suffixes such as .rds or _filtered.
.canon_sample_name <- function(x) {
    s <- tolower(trimws(as.character(x)))
    s <- gsub("\\.rds$", "", s, ignore.case = TRUE)
    s <- gsub("_harmonized_filtered$|_filtered$|_with_doublet_calls$|_singlets$", "", s, ignore.case = TRUE)
    s
}

# ── 1.10 Usage printer ────────────────────────────────────────────────────────
.print_usage <- function(default_config) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg_idx <- grep("^--file=", args)
    script_name <- if (length(file_arg_idx) == 1L) {
        basename(sub("^--file=", "", args[file_arg_idx]))
    } else {
        "QC_scdbl_Combined.R"
    }
    cat("Usage:\n")
    cat(sprintf("  Rscript %s [--config PATH] [--stage all|qc|doublet] [OPTIONS]\n\n", script_name))
    cat("Options:\n")
    cat(sprintf("  --config PATH            Config DCF file (default: %s)\n", default_config))
    cat("  --stage STR              all | qc | doublet  (default: all)\n")
    cat("  --rds_dir PATH           Override raw RDS input directory\n")
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

# ── 1.11 Parameter Metadata Checkpointing ────────────────────────────────────
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

.save_run_metadata <- function(target_rds, current_params, extra = list()) {
    meta_path <- .get_meta_path(target_rds)
    meta_data <- c(list(
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        param_hash = digest::digest(current_params),
        params = current_params
    ), extra)
    jsonlite::write_json(meta_data, meta_path, auto_unbox = TRUE, pretty = TRUE)
}

# ── 1.12 10x multiplet-table helpers (v2/v3 and v4) ─────────────────────────
# Table source: 10x Genomics user guides (loaded cells, recovered cells,
# and multiplet rate).
#
# IMPORTANT:
#   The published values are approximate ("~"). We therefore learn a smooth
#   pattern from the table and use it for prediction, including extrapolation.
#   In this pipeline ghost removal happens first, then the post-filter cell
#   count is used to estimate DBR, and scDblFinder runs on that same filtered
#   object so technical ghost cells do not distort the model.
#
# v2/v3 chemistry (default):
DBL_10X_V3_TABLE_LOADED    <- c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000)
DBL_10X_V3_TABLE_RECOVERED <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
DBL_10X_V3_TABLE_DBR       <- c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076)
#
# v4 chemistry (Single Cell 3' Reagent Kits v4, CG000731 Rev A — lower multiplet rate, ~0.58x v3):
DBL_10X_V4_TABLE_LOADED    <- c(725, 1450, 2900, 4350, 5800, 7250, 8700, 10150, 11600, 13050,
                                 14500, 15950, 17400, 18850, 20300, 21750, 23200, 24650, 26100, 27550, 29000)
DBL_10X_V4_TABLE_RECOVERED <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                                 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000)
DBL_10X_V4_TABLE_DBR       <- c(0.002, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.036,
                                 0.040, 0.044, 0.048, 0.052, 0.056, 0.060, 0.064, 0.068, 0.072, 0.076, 0.080)
#
# 5' chemistry (Single Cell 5' Reagent Kit — official 10x table, same curve as v4):
DBL_10X_5P_TABLE_LOADED    <- c(725, 1450, 2900, 4350, 5800, 7250, 8700, 10150, 11600, 13050,
                                 14500, 15950, 17400, 18850, 20300, 21750, 23200, 24650, 26100, 27550, 29000)
DBL_10X_5P_TABLE_RECOVERED <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                                 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000)
DBL_10X_5P_TABLE_DBR       <- c(0.002, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.036,
                                 0.040, 0.044, 0.048, 0.052, 0.056, 0.060, 0.064, 0.068, 0.072, 0.076, 0.080)

# Backward-compatible aliases so any existing code referencing old names still works:
DBL_10X_TABLE_LOADED    <- DBL_10X_V3_TABLE_LOADED
DBL_10X_TABLE_RECOVERED <- DBL_10X_V3_TABLE_RECOVERED
DBL_10X_TABLE_DBR       <- DBL_10X_V3_TABLE_DBR

# Generic piecewise-linear estimator over a table x-axis.
.estimate_dbr_from_table <- function(cell_count, tbl_counts, tbl_dbr) {
    x <- .as_num(cell_count, NA_real_)
    if (is.na(x) || x <= 0) return(NA_real_)

    # In-range: piecewise linear interpolation on the table.
    if (x >= min(tbl_counts) && x <= max(tbl_counts)) {
        return(as.numeric(stats::approx(
            x = tbl_counts,
            y = tbl_dbr,
            xout = x,
            method = "linear",
            rule = 2
        )$y))
    }

    # Out-of-range: extrapolate using edge segment slope (pattern continuation).
    if (x < min(tbl_counts)) {
        slope_lo <- (tbl_dbr[2] - tbl_dbr[1]) /
                    (tbl_counts[2] - tbl_counts[1])
        d <- tbl_dbr[1] + slope_lo * (x - tbl_counts[1])
    } else {
        n <- length(tbl_counts)
        slope_hi <- (tbl_dbr[n] - tbl_dbr[n - 1]) /
                    (tbl_counts[n] - tbl_counts[n - 1])
        d <- tbl_dbr[n] + slope_hi * (x - tbl_counts[n])
    }

    # Keep within valid scDblFinder boundary used elsewhere in this script.
    as.numeric(min(max(d, 0), 0.30))
}

# v2/v3 estimator (original behaviour)
.estimate_10x_dbr_from_count <- function(cell_count) {
    .estimate_dbr_from_table(cell_count, DBL_10X_V3_TABLE_RECOVERED, DBL_10X_V3_TABLE_DBR)
}

# v4 estimator
.estimate_10x_v4_dbr_from_count <- function(cell_count) {
    .estimate_dbr_from_table(cell_count, DBL_10X_V4_TABLE_RECOVERED, DBL_10X_V4_TABLE_DBR)
}

# 5' estimator
.estimate_10x_5p_dbr_from_count <- function(cell_count) {
    .estimate_dbr_from_table(cell_count, DBL_10X_5P_TABLE_RECOVERED, DBL_10X_5P_TABLE_DBR)
}

# Dispatcher: pick the right estimator based on chemistry version string.
.estimate_10x_dbr_by_chemistry <- function(cell_count, chemistry = "v3") {
    chem <- tolower(trimws(chemistry))
    if (chem == "v4") {
        .estimate_10x_v4_dbr_from_count(cell_count)
    } else if (chem %in% c("5p", "5prime", "5'")) {
        .estimate_10x_5p_dbr_from_count(cell_count)
    } else {
        # v2 and v3 share the same 10x table
        .estimate_10x_dbr_from_count(cell_count)
    }
}

# Parse an sc_protocol metadata string into a normalised platform + chemistry.
#
# Known protocol strings and their mappings:
#   "10x_3_v3.1" / "10x_3_v3"  → platform=10x, chemistry=v3
#   "10x_3_v2"   / "10x_v2"    → platform=10x, chemistry=v2
#   "10x_3_v1"                  → platform=10x, chemistry=v2  (v1 ≈ v2 multiplet rates)
#   "10x_v4"                    → platform=10x, chemistry=v4
#   "10X_5"                     → platform=10x, chemistry=5p
#   "Dropseq" / "Drop-seq"     → platform=dropseq
#   "PIPseq_V_T10_3"           → platform=pipseq
#
# Returns: list(platform = <string>, chemistry = <string or NA>)
.parse_sc_protocol <- function(protocol_string) {
    s <- tolower(trimws(protocol_string))

    # Drop-seq / Dropseq → dropseq
    if (grepl("drop.?seq", s)) {
        return(list(platform = "dropseq", chemistry = NA_character_))
    }

    # PIPseq → pipseq
    if (grepl("pipseq", s)) {
        return(list(platform = "pipseq", chemistry = NA_character_))
    }

    # 10x platforms
    if (grepl("^10x", s)) {
        # 10x 5' kit (e.g., "10x_5") → dedicated 5' multiplet rate table
        if (grepl("_5($|[_'])", s)) {
            return(list(platform = "10x", chemistry = "5p"))
        }

        # Extract version: v1, v2, v3, v3.1, v4
        ver_match <- regmatches(s, regexpr("v[0-9]+(\\.[0-9]+)?", s))
        if (length(ver_match) == 1 && nzchar(ver_match)) {
            chem <- ver_match
            if (chem == "v3.1") chem <- "v3"   # v3.1 shares v3 multiplet table
            if (chem == "v1")   chem <- "v2"   # v1 closest to v2 table
            return(list(platform = "10x", chemistry = chem))
        }

        # 10x but no recognisable version → default to v3
        return(list(platform = "10x", chemistry = "v3"))
    }

    # Unknown platform — return as-is
    list(platform = s, chemistry = NA_character_)
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Configuration & Path Setup
# ══════════════════════════════════════════════════════════════════════════════

# Run only if called as a script
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
BASE_OUT_DIR    <- .resolve_path(cfg$output_dir    %||% NULL,                  getwd())

summary_subdir  <- if (!is.null(cfg$summary_subdir)  && nzchar(cfg$summary_subdir))  cfg$summary_subdir  else "Summary"
filtered_subdir <- if (!is.null(cfg$filtered_subdir) && nzchar(cfg$filtered_subdir)) cfg$filtered_subdir else "qc_filtered_rds"
doublet_subdir  <- if (!is.null(cfg$doublet_subdir)  && nzchar(cfg$doublet_subdir))  cfg$doublet_subdir  else "doublet_filtered_rds"
integrated_summary_subdir <- if (!is.null(cfg$integrated_summary_subdir)) trimws(cfg$integrated_summary_subdir) else ""

QC_DIR              <- file.path(BASE_OUT_DIR, summary_subdir)
FILTERED_DIR        <- file.path(BASE_OUT_DIR, filtered_subdir)
DOUBLET_DIR         <- file.path(BASE_OUT_DIR, doublet_subdir)
DOUBLET_SUMMARY_DIR <- file.path(DOUBLET_DIR, "Summary")
INTEGRATED_SUMMARY_DIR <- if (nzchar(integrated_summary_subdir)) file.path(BASE_OUT_DIR, integrated_summary_subdir) else BASE_OUT_DIR

default_log   <- file.path(BASE_OUT_DIR, "logs", "qc_pipeline.log")
log_candidate <- if (!is.null(cfg$log_file) && nzchar(cfg$log_file)) cfg$log_file else default_log
log_file_path <- if (grepl("^(/|[A-Za-z]:)", log_candidate)) {
    normalizePath(log_candidate, mustWork = FALSE)
} else {
    normalizePath(file.path(BASE_OUT_DIR, log_candidate), mustWork = FALSE)
}

# ── scDblFinder tuning parameters ────────────────────────────────────────────
DBL_BATCH_COL      <- cfg$dbl_batch_col      %||% "sample_id"
DBL_MIN_COUNT      <- as.integer(cfg$dbl_min_count   %||% 100)
DBL_MIN_FEATURE    <- as.integer(cfg$dbl_min_feature %||% 50)
DBL_UMAP_DIMS      <- 1:as.integer(cfg$dbl_umap_dims %||% 20)
DBL_DEFAULT_PLATFORM  <- tolower(cfg$dbl_default_platform  %||% "10x")
DBL_DEFAULT_CHEMISTRY <- tolower(cfg$dbl_default_chemistry %||% "v3")

# ── Per-sample DBR table parsed from DCF keys ─────────────────────────────────
# Keys follow the convention:
#   sample.<sample_name>.platform  = 10x | dropseq | indrops | ...
#   sample.<sample_name>.chemistry = v3 | v4  (10x chemistry version; default v3)
#   sample.<sample_name>.dbr       = 0.08   (omit for 10x auto-estimation)
#
# Example in qc_config.dcf:
#   sample.my_sample_A.platform: 10x
#   sample.my_sample_A.chemistry: v4
#   sample.my_sample_A.dbr: 0.08
#   sample.my_sample_B.platform: dropseq
#   sample.my_sample_B.dbr: 0.05
#
# The parser scans all cfg keys that start with "sample." and builds a lookup
# list: DBL_SAMPLE_CFG[["my_sample_A"]] = list(platform="10x", chemistry="v4", dbr=0.08)
DBL_SAMPLE_CFG <- list()
sample_keys    <- grep("^sample\\.", names(cfg), value = TRUE)
for (sk in sample_keys) {
    key_match <- regexec("^sample\\.(.+)\\.(platform|chemistry|dbr)$", sk)
    key_parts <- regmatches(sk, key_match)[[1]]
    if (length(key_parts) != 3) next
    sname <- .canon_sample_name(key_parts[2])
    field <- key_parts[3]
    if (is.null(DBL_SAMPLE_CFG[[sname]])) DBL_SAMPLE_CFG[[sname]] <- list()
    DBL_SAMPLE_CFG[[sname]][[field]] <- cfg[[sk]]
}

# ── Platform-specific DBR defaults (non-10x) ─────────────────────────────────
# Keys in DCF:  platform_dbr.<platform>: <fraction>
# e.g.  platform_dbr.dropseq: 0.049
#       platform_dbr.pipseq:  0.07
DBL_PLATFORM_DBR  <- list()
pdbr_keys         <- grep("^platform_dbr\\.", names(cfg), value = TRUE)
for (pk in pdbr_keys) {
    pname <- tolower(sub("^platform_dbr\\.", "", pk))
    DBL_PLATFORM_DBR[[pname]] <- as.numeric(cfg[[pk]])
}

# ── Discovery settings ────────────────────────────────────────────────────────
rds_pattern         <- cfg$rds_pattern         %||% "\\.rds$"
recursive_discovery <- .as_bool(cfg$recursive_discovery %||% "FALSE", FALSE)

# ── Create output directories ─────────────────────────────────────────────────
for (d in c(QC_DIR, FILTERED_DIR, DOUBLET_DIR, DOUBLET_SUMMARY_DIR, INTEGRATED_SUMMARY_DIR, dirname(log_file_path)))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Initialise shared logger ──────────────────────────────────────────────────
.init_logger(log_file_path)
run_ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

log_info("══════════════════════════════════════════════════")
log_info(sprintf("QC PIPELINE START  [stage: %s]", RUN_STAGE))
log_info(sprintf("  Config       : %s", config_path))
log_info(sprintf("  Raw RDS dir  : %s", rds_dir))
log_info("  Decisions src: config:qc_decisions_table")
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
    ifelse(bound == "upper", Inf, -Inf)
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


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — STAGE 2: QC Filtering
# ══════════════════════════════════════════════════════════════════════════════

.run_stage_qc <- function(qc_input_dir = rds_dir) {

    log_info("──────────────────────────────────────────────────")
    log_info("STAGE 2: QC FILTERING")
    log_info("──────────────────────────────────────────────────")

    # ── 5.1 Discover input RDS files ──────────────────────────────────────────
    rds_file_map <- discover_rds_files(qc_input_dir, pattern = rds_pattern, recursive = recursive_discovery)
    sample_names <- sort(names(rds_file_map))
    log_info(sprintf("Discovered %d RDS files.", length(sample_names)))
    log_info(sprintf("Datasets: %s", paste(sample_names, collapse = ", ")))

    # ── 5.2 Load QC decisions (from config table) ──────────────────────────────
    decisions_payload <- .load_qc_decisions(cfg)
    qc_decisions <- decisions_payload$data
    log_info(sprintf("Loaded QC decisions from %s", decisions_payload$source))

    # ── 5.3 Dependency guard ──────────────────────────────────────────────────
    required_vars <- c("rds_dir", "QC_DIR", "FILTERED_DIR",
                       "s_median", "s_mean", "s_max")
    missing_vars  <- required_vars[!sapply(required_vars, function(v) exists(v, inherits = TRUE))]
    if (length(missing_vars) > 0)
        stop(sprintf("Missing required variables: %s", paste(missing_vars, collapse = ", ")))
    if (!is.list(rds_file_map) || length(rds_file_map) == 0)
        stop(sprintf("No RDS files discovered in QC input directory: %s", qc_input_dir))

    default_noisy_datasets <- c() # Example defaults; replace with actual names as needed

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

    qc_basenames <- gsub("\\.rds$", "", basename(qc_decisions$Dataset_Name), ignore.case = TRUE)
    configured_datasets <- intersect(qc_basenames, sample_names)
    # Also match via _singlets suffix (doublet-filtered outputs)
    if (length(configured_datasets) < length(qc_basenames)) {
        unmatched <- setdiff(qc_basenames, sample_names)
        singlets_matches <- intersect(paste0(unmatched, "_singlets"), sample_names)
        configured_datasets <- union(configured_datasets, singlets_matches)
    }

    if (length(configured_datasets) == 0)
        log_warn("No configured datasets overlap with discovered RDS files.")

    log_info(sprintf("Loaded %d QC decision rows. Running %d dataset(s).",
                     nrow(qc_decisions), length(configured_datasets)))

    # ── 5.5 Open audit PDF ────────────────────────────────────────────────────
    qc_pdf_path <- file.path(QC_DIR, "QC_Before_After_Report.pdf")
    pdf(qc_pdf_path, width = 14, height = 6)
    on.exit(tryCatch(dev.off(), error = function(e) NULL), add = TRUE)

    summary_list <- list()

    # ── 5.6 Per-dataset QC loop ───────────────────────────────────────────────
    for (i in seq_len(nrow(qc_decisions))) {

        nm       <- gsub("\\.rds$", "", basename(qc_decisions$Dataset_Name[i]), ignore.case = TRUE)
        rds_path <- rds_file_map[[nm]] # Retrieve full path from map
        if (is.null(rds_path) && identical(normalizePath(qc_input_dir, mustWork = FALSE),
                                           normalizePath(DOUBLET_DIR, mustWork = FALSE))) {
            nm_singlets <- sub("_filtered$", "_singlets", nm)
            rds_path <- rds_file_map[[nm_singlets]]
            # Also try appending _singlets (for harmonized datasets)
            if (is.null(rds_path)) {
                nm_singlets <- paste0(nm, "_singlets")
                rds_path <- rds_file_map[[nm_singlets]]
            }
            if (!is.null(rds_path)) {
                log_info(sprintf("[QC] Reverse-mode mapping: %s -> %s", nm, nm_singlets))
                nm <- nm_singlets
            }
        }
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
            # Try to read cell counts from cached metadata to avoid loading the raw RDS
            tryCatch({
                meta <- tryCatch(jsonlite::fromJSON(.get_meta_path(filtered_path)), error = function(e) list())
                fobj <- readRDS(filtered_path)
                cells_filtered <- ncol(fobj)

                if (!is.null(meta$cells_raw) && is.finite(as.numeric(meta$cells_raw))) {
                    cells_raw <- as.integer(meta$cells_raw)
                    log_info(sprintf("  Resume: read cells_raw=%d from metadata cache.", cells_raw))
                } else {
                    obj <- readRDS(rds_path)
                    cells_raw <- ncol(obj)
                    rm(obj); gc()
                    log_info(sprintf("  Resume: metadata cache miss — loaded raw RDS for cells_raw=%d.", cells_raw))
                }

                summary_list[[length(summary_list) + 1]] <- data.frame(
                    sample = nm,
                    algo_nCount_lower = row_cfg$Lower_Count_Method,
                    median_inflation_ratio = NA_real_,
                    cells_raw = cells_raw,
                    cells_filtered = cells_filtered,
                    cells_removed = cells_raw - cells_filtered,
                    pct_cells_removed = round((cells_raw - cells_filtered) / cells_raw * 100, 1),
                    genes_raw = NA_integer_,
                    genes_filtered = nrow(fobj),
                    median_nCount_raw = NA_real_,
                    median_nCount_filtered = median(fobj$nCount_RNA),
                    mean_nCount_raw = NA_real_,
                    mean_nCount_filtered = mean(fobj$nCount_RNA),
                    max_nCount_raw = NA_real_,
                    max_nCount_filtered = max(fobj$nCount_RNA),
                    cutoff_nCount_lower = NA, cutoff_nCount_upper = NA,
                    median_nFeature_raw = NA_real_,
                    median_nFeature_filtered = median(fobj$nFeature_RNA),
                    mean_nFeature_raw = NA_real_,
                    mean_nFeature_filtered = mean(fobj$nFeature_RNA),
                    max_nFeature_raw = NA_real_,
                    max_nFeature_filtered = max(fobj$nFeature_RNA),
                    cutoff_nFeature_lower = NA, cutoff_nFeature_upper = NA,
                    median_pct_mito_raw = NA_real_,
                    median_pct_mito_filtered = median(fobj$percent.mito),
                    mean_pct_mito_raw = NA_real_,
                    mean_pct_mito_filtered = mean(fobj$percent.mito),
                    max_pct_mito_raw = NA_real_,
                    max_pct_mito_filtered = max(fobj$percent.mito),
                    cutoff_mito_lower = NA, cutoff_mito_upper = NA,
                    median_pct_ribo_raw = NA_real_,
                    median_pct_ribo_filtered = median(fobj$percent.ribo),
                    mean_pct_ribo_raw = NA_real_,
                    mean_pct_ribo_filtered = mean(fobj$percent.ribo),
                    max_pct_ribo_raw = NA_real_,
                    max_pct_ribo_filtered = max(fobj$percent.ribo)
                )
                rm(fobj); gc()
            }, error = function(e) log_warn(sprintf("Could not restore summary stats for skipped %s: %s", nm, e$message)))
            next
        }

        dataset_ok <- tryCatch({

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
                as.numeric(obj@meta.data[["nFeature_RNA"]]) >= feat_lower  &
                as.numeric(obj@meta.data[["nFeature_RNA"]]) <= feat_upper  &
                as.numeric(obj@meta.data[["nCount_RNA"]])   >= count_lower &
                as.numeric(obj@meta.data[["nCount_RNA"]])   <= count_upper &
                as.numeric(obj@meta.data[["percent.mito"]]) <= mito_upper
            ]

            if (length(keep_cells) == 0) {
                log_error(sprintf("[ZERO CELLS] %s: all %d cells removed by QC thresholds", nm, cells_raw))
                log_error(sprintf("  Thresholds: nFeature [%g, %g] | nCount [%g, %g] | mito <= %g",
                                  feat_lower, feat_upper, count_lower, count_upper, mito_upper))
                log_error(sprintf("  Data ranges: nFeature [%g, %g] | nCount [%g, %g] | mito [%g, %g]",
                                  min(meta_raw$nFeature_RNA), max(meta_raw$nFeature_RNA),
                                  min(meta_raw$nCount_RNA),   max(meta_raw$nCount_RNA),
                                  min(meta_raw$percent.mito), max(meta_raw$percent.mito)))
                log_error("  → Review qc_decisions_table thresholds for this dataset.")
                rm(obj); gc()
                next
            }

            fobj <- subset(obj, cells = keep_cells)
            cells_filtered <- ncol(fobj)
            cells_removed  <- cells_raw - cells_filtered
            pct_removed    <- if (cells_raw > 0) round(100 * cells_removed / cells_raw, 1) else NA_real_

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
            .save_run_metadata(filtered_path, current_qc_params, extra = list(
                cells_raw = cells_raw, cells_filtered = cells_filtered
            ))

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

            # Summary row (all 40+ columns — identical to original Cell 4)
            summary_list[[length(summary_list) + 1]] <- data.frame(
                sample                    = nm,
                algo_nCount_lower         = row_cfg$Lower_Count_Method,
                algo_nFeature_lower       = row_cfg$Lower_Feature_Method,
                algo_mito_upper           = row_cfg$Upper_Mito_Method,
                median_inflation_ratio    = inflation,
                cells_raw                 = cells_raw,
                cells_filtered            = cells_filtered,
                cells_removed             = cells_removed,
                pct_cells_removed         = pct_removed,
                genes_raw                 = genes_raw,
                genes_filtered            = nrow(fobj),
                median_nCount_raw         = s_median(obj$nCount_RNA),
                median_nCount_filtered    = s_median(fobj$nCount_RNA),
                mean_nCount_raw           = s_mean(obj$nCount_RNA),
                mean_nCount_filtered      = s_mean(fobj$nCount_RNA),
                max_nCount_raw            = s_max(obj$nCount_RNA),
                max_nCount_filtered       = s_max(fobj$nCount_RNA),
                cutoff_nCount_lower       = round(count_lower, 2),
                cutoff_nCount_upper       = round(count_upper, 2),
                median_nFeature_raw       = s_median(obj$nFeature_RNA),
                median_nFeature_filtered  = s_median(fobj$nFeature_RNA),
                mean_nFeature_raw         = s_mean(obj$nFeature_RNA),
                mean_nFeature_filtered    = s_mean(fobj$nFeature_RNA),
                max_nFeature_raw          = s_max(obj$nFeature_RNA),
                max_nFeature_filtered     = s_max(fobj$nFeature_RNA),
                cutoff_nFeature_lower     = round(feat_lower, 2),
                cutoff_nFeature_upper     = round(feat_upper, 2),
                median_pct_mito_raw       = s_median(obj$percent.mito),
                median_pct_mito_filtered  = s_median(fobj$percent.mito),
                mean_pct_mito_raw         = s_mean(obj$percent.mito),
                mean_pct_mito_filtered    = s_mean(fobj$percent.mito),
                max_pct_mito_raw          = s_max(obj$percent.mito),
                max_pct_mito_filtered     = s_max(fobj$percent.mito),
                cutoff_mito_lower         = 0,
                cutoff_mito_upper         = round(mito_upper, 2),
                median_pct_ribo_raw       = s_median(obj$percent.ribo),
                median_pct_ribo_filtered  = s_median(fobj$percent.ribo),
                mean_pct_ribo_raw         = s_mean(obj$percent.ribo),
                mean_pct_ribo_filtered    = s_mean(fobj$percent.ribo),
                max_pct_ribo_raw          = s_max(obj$percent.ribo),
                max_pct_ribo_filtered     = s_max(fobj$percent.ribo),
                cutoff_ribo_lower         = 0,
                cutoff_ribo_upper         = Inf,
                stringsAsFactors = FALSE
            )

            rm(obj, fobj, meta_raw, meta_filt); gc()
            TRUE
        }, error = function(e) {
            log_error(sprintf("[FAILED] %s: %s", nm, conditionMessage(e)))
            log_error(sprintf("[CALLS] %s",
                paste(sapply(sys.calls(), function(x) deparse(x[[1]])[1]), collapse = " > ")))
            FALSE
        })
    }

    # TOTALS row
    summary_table <- if (length(summary_list) > 0) do.call(rbind, summary_list) else data.frame()
    if (nrow(summary_table) > 0) {
        totals_row                   <- summary_table[1, ]; totals_row[1, ] <- NA
        totals_row$sample            <- "TOTAL"
        totals_row$cells_raw         <- sum(summary_table$cells_raw,      na.rm = TRUE)
        totals_row$cells_filtered    <- sum(summary_table$cells_filtered, na.rm = TRUE)
        totals_row$cells_removed     <- sum(summary_table$cells_removed,  na.rm = TRUE)
        totals_row$pct_cells_removed <- round(100 * totals_row$cells_removed / totals_row$cells_raw, 1)
        summary_table                <- rbind(summary_table, totals_row)
    }

    tryCatch(dev.off(), error = function(e) log_warn(paste("PDF device close:", conditionMessage(e))))
    log_info(sprintf("QC audit PDF saved: %s", qc_pdf_path))

    # Print summary and save CSV
    log_info("=== QC Summary ===")
    if (nrow(summary_table) > 0) {
        print_cols <- intersect(
            c("sample", "algo_nCount_lower", "median_inflation_ratio",
              "cells_raw", "cells_filtered", "pct_cells_removed",
              "cutoff_nCount_lower", "cutoff_nCount_upper",
              "median_nCount_raw",   "median_nCount_filtered",
              "median_pct_mito_raw", "median_pct_mito_filtered"),
            colnames(summary_table)
        )
        for (ln in capture.output(print(summary_table[, print_cols], row.names = FALSE))) log_info(ln)
    } else {
        log_warn("No datasets processed.")
    }

    qc_summary_csv <- file.path(QC_DIR, "qc_summary_detailed.csv")
    write.csv(summary_table, qc_summary_csv, row.names = FALSE)
    log_info(sprintf("Full summary (%d rows, %d columns) saved: %s",
                     nrow(summary_table), ncol(summary_table), qc_summary_csv))

    # ── 5.7 Full QC Dashboard (Cell 5) ───────────────────────────────────────
    if (nrow(summary_table) == 0) {
        log_warn("Summary table empty — skipping dashboard.")
    } else {

        plot_data <- summary_table[summary_table$sample != "TOTAL", ]

        # Shorten labels for readability
        plot_data$sample_short <- plot_data$sample
        plot_data$sample_short <- gsub("d10_1016_j_",             "", plot_data$sample_short)
        plot_data$sample_short <- gsub("d10_1126_sciadv_",        "", plot_data$sample_short)
        plot_data$sample_short <- gsub("d10_1038_",               "", plot_data$sample_short)
        plot_data$sample_short <- gsub("dno_doi_kidney_organoid_", "kidney_", plot_data$sample_short)

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

        # P3: Median inflation ratio
        p3_inf <- ggplot(plot_data, aes(x = sample_short, y = median_inflation_ratio,
                                        fill = median_inflation_ratio > 1.5)) +
            geom_col(alpha = 0.85) +
            geom_hline(yintercept = 1.0, linetype = "solid",  color = "grey40", linewidth = 0.4) +
            geom_hline(yintercept = 1.4, linetype = "dashed", color = "orange", linewidth = 0.5) +
            geom_hline(yintercept = 1.5, linetype = "dashed", color = "red",    linewidth = 0.6) +
            scale_fill_manual(values = c("FALSE" = "#41b6c4", "TRUE" = "#e31a1c"),
                              labels = c("OK (\u22641.5x)", "Over-filtered (>1.5x)")) +
            labs(title = "Median nCount Inflation Ratio (post/pre)",
                 subtitle = "Target: 1.1\u20131.4\u00d7  |  Orange: 1.4\u00d7  |  Red: 1.5\u00d7 (over-filtering threshold)",
                 y = "Inflation Ratio", x = "", fill = "") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")

        # P4: nCount median shift (line plot)
        df_ncount      <- plot_data[, c("sample_short", "median_nCount_raw", "median_nCount_filtered")]
        df_long_ncount <- melt(df_ncount, id.vars = "sample_short")
        p4 <- ggplot(df_long_ncount, aes(x = sample_short, y = value, group = variable, color = variable)) +
            geom_line(linewidth = 0.8) + geom_point(size = 2) +
            scale_color_manual(values = c("median_nCount_raw" = "#fee0d2", "median_nCount_filtered" = "#de2d26"),
                               labels = c("Before", "After")) +
            labs(title = "Signal Recovery: Median nCount_RNA", y = "Counts", x = "", color = "Stage") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")

        # P5: nFeature median shift (line plot)
        df_nfeat      <- plot_data[, c("sample_short", "median_nFeature_raw", "median_nFeature_filtered")]
        df_long_nfeat <- melt(df_nfeat, id.vars = "sample_short")
        p5 <- ggplot(df_long_nfeat, aes(x = sample_short, y = value, group = variable, color = variable)) +
            geom_line(linewidth = 0.8) + geom_point(size = 2) +
            scale_color_manual(values = c("median_nFeature_raw" = "#e0f3db", "median_nFeature_filtered" = "#43a2ca"),
                               labels = c("Before", "After")) +
            labs(title = "Gene Coverage: Median nFeature_RNA", y = "Genes", x = "", color = "Stage") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")

        # P6: Mito/ribo contamination reduction (facet bar)
        df_qc_meta <- plot_data[, c("sample_short",
                                    "median_pct_mito_raw", "median_pct_mito_filtered",
                                    "median_pct_ribo_raw", "median_pct_ribo_filtered")]
        df_long_qc        <- melt(df_qc_meta, id.vars = "sample_short")
        df_long_qc$metric <- ifelse(grepl("mito", df_long_qc$variable), "Mitochondrial %", "Ribosomal %")
        df_long_qc$stage  <- ifelse(grepl("raw",  df_long_qc$variable), "Before", "After")
        df_long_qc$stage  <- factor(df_long_qc$stage, levels = c("Before", "After"))
        p6 <- ggplot(df_long_qc, aes(x = sample_short, y = value, fill = stage)) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
            facet_wrap(~metric, scales = "free_y", ncol = 1) +
            scale_fill_manual(values = c("Before" = "#bdbdbd", "After" = "#756bb1")) +
            labs(title = "Reduction in Contamination Metrics", y = "Percentage (%)", x = "", fill = "Stage") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "bottom")

        # Ridge plots — reload raw + filtered RDS for per-dataset distributions
        dist_df_list <- list()
        for (nm in configured_datasets) {
            raw_path  <- rds_file_map[[nm]]
            filt_path <- file.path(FILTERED_DIR, paste0(nm, "_filtered.rds"))

            if (is.null(raw_path) || !file.exists(raw_path) || !file.exists(filt_path)) next

            log_info(sprintf("Loading ridge data for: %s", nm))
            raw_obj  <- tryCatch(readRDS(raw_path),  error = function(e) {
                log_warn(sprintf("Could not load raw ridge data for %s: %s", nm, conditionMessage(e))); NULL })
            filt_obj <- tryCatch(readRDS(filt_path), error = function(e) {
                log_warn(sprintf("Could not load filtered ridge data for %s: %s", nm, conditionMessage(e))); NULL })

            if (is.null(raw_obj) || is.null(filt_obj)) next

            short_nm <- nm
            short_nm <- gsub("d10_1016_j_",             "", short_nm)
            short_nm <- gsub("d10_1126_sciadv_",        "", short_nm)
            short_nm <- gsub("d10_1038_",               "", short_nm)
            short_nm <- gsub("dno_doi_kidney_organoid_", "kidney_", short_nm)

            dist_df_list[[nm]] <- rbind(
                data.frame(sample = short_nm, status = "Before",
                           nCount_RNA = raw_obj$nCount_RNA, nFeature_RNA = raw_obj$nFeature_RNA),
                data.frame(sample = short_nm, status = "After",
                           nCount_RNA = filt_obj$nCount_RNA, nFeature_RNA = filt_obj$nFeature_RNA)
            )
            rm(raw_obj, filt_obj); gc()
        }

        if (length(dist_df_list) > 0) {
            dist_df              <- do.call(rbind, dist_df_list)
            dist_df$log1p_nCount <- log1p(dist_df$nCount_RNA)
            dist_df$log1p_nFeat  <- log1p(dist_df$nFeature_RNA)
            dist_df$status       <- factor(dist_df$status, levels = c("Before", "After"))

            cutoffs_df <- data.frame(
                sample       = plot_data$sample_short,
                nCount_lower = log1p(plot_data$cutoff_nCount_lower),
                nCount_upper = log1p(plot_data$cutoff_nCount_upper),
                nFeat_lower  = log1p(plot_data$cutoff_nFeature_lower),
                nFeat_upper  = log1p(plot_data$cutoff_nFeature_upper)
            )
            cutoffs_df$status <- factor("Before", levels = c("Before", "After"))

            # Per-dataset density: nCount
            p_counts_indiv <- ggplot(dist_df, aes(x = log1p_nCount, fill = status)) +
                geom_density(alpha = 0.5, linewidth = 0.2) +
                geom_vline(data = cutoffs_df, aes(xintercept = nCount_lower),
                           linetype = "dashed", color = "red", linewidth = 0.5) +
                geom_vline(data = cutoffs_df, aes(xintercept = nCount_upper),
                           linetype = "dashed", color = "red", linewidth = 0.5) +
                facet_wrap(~sample, scales = "free_y", ncol = 4) +
                scale_fill_manual(values = c("Before" = "#bdbdbd", "After" = "#e31a1c")) +
                theme_minimal() +
                labs(title = "log1p(nCount_RNA) Distribution per Dataset",
                     x = "log1p(nCount_RNA)", y = "Density") +
                theme(legend.position = "top", strip.text = element_text(size = 8))

            # Per-dataset density: nFeature
            p_genes_indiv <- ggplot(dist_df, aes(x = log1p_nFeat, fill = status)) +
                geom_density(alpha = 0.5, linewidth = 0.2) +
                geom_vline(data = cutoffs_df, aes(xintercept = nFeat_lower),
                           linetype = "dashed", color = "red", linewidth = 0.5) +
                geom_vline(data = cutoffs_df, aes(xintercept = nFeat_upper),
                           linetype = "dashed", color = "red", linewidth = 0.5) +
                facet_wrap(~sample, scales = "free_y", ncol = 4) +
                scale_fill_manual(values = c("Before" = "#bdbdbd", "After" = "#1f78b4")) +
                theme_minimal() +
                labs(title = "log1p(nFeature_RNA) Distribution per Dataset",
                     x = "log1p(nFeature_RNA)", y = "Density") +
                theme(legend.position = "top", strip.text = element_text(size = 8))

            # Global maxima for ridge x-axis alignment
            global_max_nCount <- max(dist_df$log1p_nCount, na.rm = TRUE)
            global_max_nFeat  <- max(dist_df$log1p_nFeat,  na.rm = TRUE)

            # Ridge: nCount across datasets
            p_counts_ridge <- ggplot(dist_df, aes(x = log1p_nCount, y = sample, fill = sample)) +
                geom_density_ridges(alpha = 0.7, scale = 1.5, color = "white", linewidth = 0.4) +
                geom_vline(xintercept = global_max_nCount, color = "black", linetype = "dashed", linewidth = 0.6) +
                scale_x_continuous(limits = c(0, global_max_nCount * 1.05)) +
                facet_wrap(~status, ncol = 2) +
                theme_minimal() +
                labs(title = sprintf("Vertical Across Datasets: log1p(nCount) [Global Max: %.2f]", global_max_nCount),
                     x = "log1p(nCount_RNA)", y = "") +
                theme(axis.text.y = element_text(size = 8), legend.position = "none")

            # Ridge: nFeature across datasets
            p_genes_ridge <- ggplot(dist_df, aes(x = log1p_nFeat, y = sample, fill = sample)) +
                geom_density_ridges(alpha = 0.7, scale = 1.5, color = "white", linewidth = 0.4) +
                geom_vline(xintercept = global_max_nFeat, color = "black", linetype = "dashed", linewidth = 0.6) +
                scale_x_continuous(limits = c(0, global_max_nFeat * 1.05)) +
                facet_wrap(~status, ncol = 2) +
                theme_minimal() +
                labs(title = sprintf("Vertical Across Datasets: log1p(nFeature) [Global Max: %.2f]", global_max_nFeat),
                     x = "log1p(nFeature_RNA)", y = "") +
                theme(axis.text.y = element_text(size = 8), legend.position = "none")

            # Full 7-panel composite layout
            final_layout <- (p1 | p2) / (p3_inf | p4) / (p5 | p6) /
                            p_counts_indiv / p_genes_indiv / p_counts_ridge / p_genes_ridge +
                plot_layout(heights = c(1, 1, 1.5, 2.5, 2.5, 2, 2)) +
                plot_annotation(
                    title    = "QC \u2014 Comprehensive Single-Cell QC Board",
                    subtitle = "CSV-driven thresholds | Violin + scatter PNGs per dataset | Inflation ratio diagnostic",
                    theme    = theme(
                        plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
                        plot.subtitle = element_text(size = 11, hjust = 0.5)
                    )
                )
        } else {
            # Fallback: 3-panel if ridge data unavailable
            final_layout <- (p1 | p2) / (p3_inf | p4) / (p5 | p6) +
                plot_layout(heights = c(1, 1, 1.5)) +
                plot_annotation(title = "QC \u2014 Comprehensive Single-Cell QC Board")
        }

        options(repr.plot.width = 18, repr.plot.height = 42)
        print(final_layout)

        ggsave(file.path(QC_DIR, "summary_qc_full_dashboard.png"),
               plot = final_layout, width = 18, height = 42, dpi = 300, bg = "white")
        log_info(sprintf("Full dashboard saved: %s", file.path(QC_DIR, "summary_qc_full_dashboard.png")))
    }

} # end .run_stage_qc


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — STAGE 1: Doublet Detection (scDblFinder)
# Input: Raw RDS files in rds_dir
# ══════════════════════════════════════════════════════════════════════════════

.run_stage_doublet <- function(doublet_input_dir = FILTERED_DIR) {

    log_info("──────────────────────────────────────────────────")
    log_info("STAGE 1: DOUBLET DETECTION")
    log_info(sprintf("  Input  : %s", doublet_input_dir))
    log_info(sprintf("  Output : %s", DOUBLET_DIR))
    log_info("──────────────────────────────────────────────────")

    # ── 6.1 Discover QC-filtered RDS files ───────────────────────────────────
    filtered_files <- list.files(doublet_input_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(filtered_files) == 0)
        stop(sprintf("No input RDS files found in: %s", doublet_input_dir))

    log_info(sprintf("Found %d input RDS files.", length(filtered_files)))

    # ── 6.2 Per-sample DBR config (read from DCF in Section 2) ──────────────
    # DBL_SAMPLE_CFG is already populated at startup from qc_config.dcf keys:
    #   sample.<sample_name>.platform  =  10x | dropseq | indrops | ...
    #   sample.<sample_name>.chemistry =  v3 | v4  (default: v3)
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
    on.exit(tryCatch(dev.off(), error = function(e) NULL), add = TRUE)

    dbl_summary_stats <- list()
    raw_cell_counts   <- list()  # cache raw ncol for Section 8

    # ── 6.4 Per-sample doublet detection loop ─────────────────────────────────
    for (file_path in filtered_files) {

        sample_nm <- gsub(
            "_harmonized_filtered\\.rds$|_filtered\\.rds$|\\.rds$",
            "",
            basename(file_path)
        )
        log_info(sprintf("[DOUBLET] Processing: %s", sample_nm))

        # Checkpoint / Resume logic (Parameter-Aware)
        singlets_path <- file.path(DOUBLET_DIR, paste0(sample_nm, "_singlets.rds"))
        dbl_run_params <- list(
            platform = DBL_DEFAULT_PLATFORM,
            chemistry = DBL_DEFAULT_CHEMISTRY,
            min_count = DBL_MIN_COUNT,
            min_feature = DBL_MIN_FEATURE,
            sample_cfg = DBL_SAMPLE_CFG[[.canon_sample_name(sample_nm)]]
        )
        # Parameter-aware skip: use .meta.json if available, else fall back to
        # simple file-existence check (preserves pre-upgrade outputs).
        skip_doublet <- .should_skip_run(singlets_path, dbl_run_params, cfg$force_overwrite)
        if (!skip_doublet && !.as_bool(cfg$force_overwrite, FALSE) &&
            file.exists(singlets_path) && file.info(singlets_path)$size > 0 &&
            !file.exists(.get_meta_path(singlets_path))) {
            skip_doublet <- TRUE
            log_info(sprintf("  [COMPAT] No .meta.json for %s — skipping via legacy file-existence check.", sample_nm))
        }
        if (skip_doublet) {
            log_info(sprintf("[DOUBLET SKIPPED] %s already processed with identical parameters.", sample_nm))
            # Recover cell counts from saved outputs if possible
            dbl_stats <- tryCatch({
                dbl_full_path <- file.path(DOUBLET_DIR, paste0(sample_nm, "_with_doublet_calls.rds"))
                if (file.exists(dbl_full_path)) {
                    obj_full <- readRDS(dbl_full_path)
                    n_total <- ncol(obj_full)
                    n_doublets <- sum(obj_full$scDblFinder.class == "doublet", na.rm = TRUE)
                    pct_dbl <- round((n_doublets / n_total) * 100, 2)
                    rm(obj_full); gc()
                    list(total = n_total, doublets = n_doublets, pct = pct_dbl)
                } else {
                    list(total = NA, doublets = NA, pct = NA)
                }
            }, error = function(e) list(total = NA, doublets = NA, pct = NA))
            dbl_summary_stats[[sample_nm]] <- data.frame(
                Sample           = sample_nm,
                Total_Cells      = dbl_stats$total,
                Doublets_Found   = dbl_stats$doublets,
                Percent_Doublet  = dbl_stats$pct,
                DBR_Used         = NA_real_,
                Cells_Loaded_Used = NA_real_,
                Cells_Recovered_Used = NA_real_,
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
        raw_cell_counts[[sample_nm]] <- n_before

        # S4SXP fix: force meta.data to base data.frame (Seurat v5 compatibility)
        obj@meta.data <- as.data.frame(obj@meta.data)

        batch_levels_meta <- character(0)
        if (DBL_BATCH_COL %in% colnames(obj@meta.data)) {
            batch_levels_meta <- unique(as.character(na.omit(obj@meta.data[[DBL_BATCH_COL]])))
            batch_levels_meta <- batch_levels_meta[nzchar(trimws(batch_levels_meta))]
        }
        has_batch <- length(batch_levels_meta) > 1

        # ── Resolve platform, chemistry, and DBR for this sample ──────────
        # Priority chain:
        #   1. Per-sample config override  (sample.<name>.platform / .chemistry / .dbr)
        #   2. sc_protocol metadata column (auto-parsed from the RDS object)
        #   3. Global defaults             (dbl_default_platform / dbl_default_chemistry)
        #
        # For 10x platforms the DBR is estimated from the chemistry-specific
        # multiplet rate table.  For non-10x platforms the DBR is looked up in
        # DBL_PLATFORM_DBR (config keys platform_dbr.<name>).
        # The pipeline never stops due to a missing value — it degrades gracefully.
        dbr               <- NULL
        dbr_source        <- "AUTO_10X"
        platform          <- DBL_DEFAULT_PLATFORM    # default: "10x"
        chemistry         <- DBL_DEFAULT_CHEMISTRY   # default: "v3"
        cells_loaded_used <- NA_real_
        cells_recovered_used <- NA_real_

        # (A) Auto-detect from sc_protocol metadata ────────────────────────
        proto_raw <- NULL
        if ("sc_protocol" %in% colnames(obj@meta.data)) {
            proto_vals <- unique(na.omit(obj@meta.data$sc_protocol))
            if (length(proto_vals) == 1L) {
                proto_raw <- as.character(proto_vals)
            } else if (length(proto_vals) > 1L) {
                proto_raw <- as.character(proto_vals[1])
                log_warn(sprintf(
                    "  Multiple sc_protocol values for '%s': %s — using first: '%s'",
                    sample_nm, paste(proto_vals, collapse = ", "), proto_raw))
            }
        }

        chemistry_display <- chemistry

        if (!is.null(proto_raw)) {
            parsed    <- .parse_sc_protocol(proto_raw)
            platform  <- parsed$platform
            if (!is.na(parsed$chemistry)) chemistry <- parsed$chemistry
            chemistry_display <- if (platform == "10x") chemistry else "N/A"
            dbr_source <- "SC_PROTOCOL_AUTO"
            log_info(sprintf("  sc_protocol='%s' -> platform=%s, chemistry=%s",
                             proto_raw, platform,
                             chemistry_display))
        } else {
            log_warn(sprintf(
                "  No sc_protocol metadata for '%s' — using global defaults (platform=%s, chemistry=%s)",
                sample_nm, platform, chemistry))
        }

        # (B) Per-sample config override (wins over sc_protocol) ──────────
        sample_cfg <- DBL_SAMPLE_CFG[[.canon_sample_name(sample_nm)]]
        if (!is.null(sample_cfg)) {
            if (!is.null(sample_cfg$platform) && nzchar(sample_cfg$platform)) {
                platform <- tolower(trimws(sample_cfg$platform))
                log_info(sprintf("  Config override: platform=%s", platform))
            }
            if (!is.null(sample_cfg$chemistry) && nzchar(sample_cfg$chemistry)) {
                chemistry <- tolower(trimws(sample_cfg$chemistry))
                log_info(sprintf("  Config override: chemistry=%s", chemistry))
            }
        }

        # Validate chemistry for 10x
        if (platform == "10x" && !chemistry %in% c("v2", "v3", "v4", "5p")) {
            log_warn(sprintf(
                "  Unknown 10x chemistry '%s' for '%s' — defaulting to v3",
                chemistry, sample_nm))
            chemistry <- "v3"
        }
        chemistry_display <- if (platform == "10x") chemistry else "N/A"

        # (C) Resolve DBR ─────────────────────────────────────────────────
        # Check for explicit user-supplied DBR first (config override)
        dbr_user <- NA_real_
        if (!is.null(sample_cfg) && !is.null(sample_cfg$dbr) && nzchar(sample_cfg$dbr)) {
            dbr_user <- suppressWarnings(as.numeric(sample_cfg$dbr))
        }

        if (platform == "10x") {
            if (!is.na(dbr_user)) {
                if (dbr_user > 0 && dbr_user <= 0.3) {
                    dbr <- dbr_user; dbr_source <- "USER_OVERRIDE_10X"
                } else {
                    log_warn(sprintf(
                        "  Invalid DBR %.4f for '%s' (must be in (0, 0.30]) — falling back to AUTO",
                        dbr_user, sample_nm))
                }
            } else if (has_batch) {
                dbr_source <- sprintf("TENX_%s_TABLE_PER_BATCH_FILTERED_CELLS", toupper(chemistry))
                cells_recovered_used <- ncol(obj)
            } else {
                # Automatic estimation from 10x multiplet rate table
                cells_count_for_dbr <- ncol(obj)
                est_dbr <- .estimate_10x_dbr_by_chemistry(cells_count_for_dbr, chemistry)
                if (!is.na(est_dbr) && est_dbr > 0 && est_dbr <= 0.3) {
                    dbr <- est_dbr
                    cells_recovered_used <- cells_count_for_dbr
                    dbr_source <- sprintf("TENX_%s_TABLE_FROM_FILTERED_CELLS", toupper(chemistry))
                } else {
                    dbr_source <- "AUTO_10X_FALLBACK"
                    log_warn(sprintf(
                        "  Could not estimate 10x %s table DBR for '%s' from %d cells — scDblFinder auto",
                        chemistry, sample_nm, cells_count_for_dbr))
                }
            }
        } else {
            # Non-10x platform — priority: user override → platform default → 10x fallback
            if (!is.na(dbr_user) && dbr_user > 0 && dbr_user <= 0.3) {
                dbr <- dbr_user; dbr_source <- "USER_OVERRIDE"
            } else if (!is.null(DBL_PLATFORM_DBR[[platform]])) {
                dbr <- DBL_PLATFORM_DBR[[platform]]
                dbr_source <- sprintf("PLATFORM_DEFAULT_%s", toupper(platform))
                log_info(sprintf("  Using platform default DBR=%.4f for %s", dbr, platform))
            } else {
                log_warn(sprintf(
                    "  Platform '%s' for '%s' has no default DBR — falling back to 10x AUTO",
                    platform, sample_nm))
                cells_count_for_dbr <- ncol(obj)
                est_dbr <- .estimate_10x_dbr_by_chemistry(cells_count_for_dbr, chemistry)
                if (!is.na(est_dbr) && est_dbr > 0 && est_dbr <= 0.3) {
                    dbr <- est_dbr
                    cells_recovered_used <- cells_count_for_dbr
                }
                dbr_source <- "AUTO_10X_FALLBACK"
            }
        }

        if (has_batch && platform == "10x" && is.na(dbr_user)) {
            log_info(sprintf("  Platform: %s | Chemistry: %s | DBR: per-batch auto [%s]",
                             platform, chemistry_display, dbr_source))
            log_info(sprintf("  10x %s merged input before per-batch split (filtered_cells_total): %.0f",
                             chemistry, ncol(obj)))
        } else {
            log_info(sprintf("  Platform: %s | Chemistry: %s | DBR: %s [%s]",
                             platform, chemistry_display,
                             ifelse(is.null(dbr), "AUTO", sprintf("%.4f", dbr)),
                             dbr_source))
        }

        if (!is.na(cells_recovered_used) && !(has_batch && platform == "10x" && is.na(dbr_user))) {
            tbl_max <- switch(chemistry,
                              v4 = max(DBL_10X_V4_TABLE_RECOVERED),
                              `5p` = max(DBL_10X_5P_TABLE_RECOVERED),
                              max(DBL_10X_V3_TABLE_RECOVERED))
            log_info(sprintf("  10x %s table input (filtered_cells_used): %.0f", chemistry, cells_recovered_used))
            if (cells_recovered_used > tbl_max) {
                log_info(sprintf("  10x pattern extrapolation: yes (above %.0f filtered cells)",
                                 tbl_max))
            }
        }

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
        class_results <- NULL
        score_results <- NULL

        if (has_batch) {
            batch_ids <- as.character(sce[[DBL_BATCH_COL]])
            batch_ids[is.na(batch_ids)] <- ""
            batch_levels <- unique(batch_ids[nzchar(batch_ids)])
            n_sub <- length(batch_levels)

            if (platform == "10x" && is.na(dbr_user)) {
                log_info(sprintf("  [BATCH MODE] %d internal samples — estimating DBR per batch.", n_sub))

                class_results <- setNames(rep(NA_character_, ncol(sce)), colnames(sce))
                score_results <- setNames(rep(NA_real_, ncol(sce)), colnames(sce))

                batch_dbrs <- setNames(rep(NA_real_, n_sub), batch_levels)
                batch_sizes <- setNames(rep(0L, n_sub), batch_levels)
                tbl_max <- switch(chemistry,
                                  v4 = max(DBL_10X_V4_TABLE_RECOVERED),
                                  `5p` = max(DBL_10X_5P_TABLE_RECOVERED),
                                  max(DBL_10X_V3_TABLE_RECOVERED))

                for (batch_nm in batch_levels) {
                    batch_idx <- which(batch_ids == batch_nm)
                    sub_sce <- sce[, batch_idx]
                    sub_n <- ncol(sub_sce)
                    batch_sizes[[batch_nm]] <- sub_n

                    sub_dbr <- .estimate_10x_dbr_by_chemistry(sub_n, chemistry)
                    if (!is.na(sub_dbr) && sub_dbr > 0 && sub_dbr <= 0.3) {
                        batch_dbrs[[batch_nm]] <- sub_dbr
                        log_info(sprintf("    [BATCH:%s] cells=%d | DBR=%.4f", batch_nm, sub_n, sub_dbr))
                    } else {
                        log_warn(sprintf("    [BATCH:%s] could not estimate 10x %s DBR from %d filtered cells — using scDblFinder auto",
                                         batch_nm, chemistry, sub_n))
                    }

                    if (sub_n > tbl_max) {
                        log_info(sprintf("    [BATCH:%s] extrapolation: yes (above %.0f filtered cells)",
                                         batch_nm, tbl_max))
                    }

                    sub_sce <- if (is.na(batch_dbrs[[batch_nm]])) scDblFinder(sub_sce)
                               else                               scDblFinder(sub_sce, dbr = batch_dbrs[[batch_nm]])

                    class_results[colnames(sub_sce)] <- as.character(sub_sce$scDblFinder.class)
                    score_results[colnames(sub_sce)] <- as.numeric(sub_sce$scDblFinder.score)
                }

                valid_batches <- is.finite(batch_dbrs)
                if (any(valid_batches)) {
                    dbr <- stats::weighted.mean(batch_dbrs[valid_batches], batch_sizes[valid_batches])
                    log_info(sprintf("  [BATCH MODE] Weighted DBR across %d internal samples: %.4f",
                                     sum(valid_batches), dbr))
                } else {
                    dbr <- NULL
                    dbr_source <- "AUTO_10X_FALLBACK"
                }
            } else {
                log_info(sprintf("  [BATCH MODE] %d internal samples — processing independently.", n_sub))
                sce <- if (is.null(dbr)) scDblFinder(sce, samples = DBL_BATCH_COL)
                       else              scDblFinder(sce, samples = DBL_BATCH_COL, dbr = dbr)
                class_results <- as.character(sce$scDblFinder.class)
                score_results <- as.numeric(sce$scDblFinder.score)
                names(class_results) <- names(score_results) <- colnames(sce)
            }
        } else {
            log_info("  [SINGLE MODE] No internal batching — running standard mode.")
            sce <- if (is.null(dbr)) scDblFinder(sce)
                   else              scDblFinder(sce, dbr = dbr)
            class_results <- as.character(sce$scDblFinder.class)
            score_results <- as.numeric(sce$scDblFinder.score)
            names(class_results) <- names(score_results) <- colnames(sce)
        }

        # Map results back to Seurat object
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
        saveRDS(obj,
                file.path(DOUBLET_DIR, paste0(sample_nm, "_with_doublet_calls.rds")))
        obj_singlets <- subset(obj, subset = scDblFinder.class == "singlet")
        saveRDS(obj_singlets,
                file.path(DOUBLET_DIR, paste0(sample_nm, "_singlets.rds")))

        log_info(sprintf("  Saved: %s_with_doublet_calls.rds | %s_singlets.rds", sample_nm, sample_nm))
        .save_run_metadata(singlets_path, dbl_run_params, extra = list(
            n_total = n_total, n_doublets = n_doublets, pct_doublet = pct_doublet
        ))

        dbl_summary_stats[[sample_nm]] <- data.frame(
            Sample           = sample_nm,
            Total_Cells      = n_total,
            Doublets_Found   = n_doublets,
            Percent_Doublet  = pct_doublet,
            DBR_Used         = ifelse(is.null(dbr), NA_real_, dbr),
            Cells_Loaded_Used = ifelse(is.na(cells_loaded_used), NA_real_, cells_loaded_used),
            Cells_Recovered_Used = ifelse(is.na(cells_recovered_used), NA_real_, cells_recovered_used),
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
                DBR_Used         = NA_real_,
                Cells_Loaded_Used = NA_real_,
                Cells_Recovered_Used = NA_real_,
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

    # Cache raw cell counts so Section 8 can skip re-reading RDS files
    if (length(raw_cell_counts) > 0) {
        rcc_df <- data.frame(
            sample = names(raw_cell_counts),
            cells_raw = as.integer(unlist(raw_cell_counts)),
            stringsAsFactors = FALSE
        )
        rcc_path <- file.path(DOUBLET_SUMMARY_DIR, "raw_cell_counts_cache.csv")
        write.csv(rcc_df, rcc_path, row.names = FALSE)
        log_info(sprintf("Raw cell count cache saved: %s", rcc_path))
    }

} # end .run_stage_doublet

# ── Stage execution dispatcher ───────────────────────────────────────────────
if (RUN_STAGE == "all") {
    .run_stage_doublet(rds_dir)
    .run_stage_qc(DOUBLET_DIR)
} else if (RUN_STAGE == "qc") {
    .run_stage_qc(DOUBLET_DIR)
} else if (RUN_STAGE == "doublet") {
    .run_stage_doublet(rds_dir)
}


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
# SECTION 8 — Integrated Doublet + QC End-to-End Summary
#
# Runs when both stages have been executed (stage = "all") or when outputs
# from both stages are already present on disk.
#
# Outputs:
#   integrated_qc_doublet_summary.csv  — one row per sample with:
#       cells_raw | cells_post_doublet | cells_post_qc |
#       pct_removed_by_doublet | pct_removed_by_qc | pct_removed_total
#   integrated_summary_dashboard.png   — 9-panel figure covering:
#       (1) 3-stage waterfall bar chart
#       (2) stacked % loss breakdown
#       (3) cumulative funnel (line)
#       (4) log1p(nCount_RNA) density overlay: raw → post-doublet → post-qc
#       (5) log1p(nFeature_RNA) density overlay
#       (6) % mito overlay (raw vs post-QC)
#       (7) doublet score distribution per sample (ridge plot)
#       (8) doublet % vs QC removal % scatter
#       (9) total cells lost per stage (pie / donut)
# ══════════════════════════════════════════════════════════════════════════════

log_info("──────────────────────────────────────────────────")
log_info("SECTION 8: INTEGRATED DOUBLET + QC SUMMARY")
log_info("──────────────────────────────────────────────────")

# ── 8.0 Guard: need at least one stage's outputs to proceed ──────────────────
qc_csv_path  <- file.path(QC_DIR, "qc_summary_detailed.csv")
dbl_csv_glob <- list.files(DOUBLET_SUMMARY_DIR, pattern = "^doublet_summary_.*\\.csv$", full.names = TRUE)
dbl_csv_path <- if (length(dbl_csv_glob) > 0) tail(sort(dbl_csv_glob), 1) else NULL

has_qc_summary  <- file.exists(qc_csv_path)
has_dbl_summary <- !is.null(dbl_csv_path) && file.exists(dbl_csv_path)

.safe_read_csv <- function(path, label) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    finfo <- file.info(path)
    if (is.na(finfo$size) || finfo$size <= 0) {
        log_warn(sprintf("Section 8: %s CSV is empty (0 bytes): %s", label, path))
        return(NULL)
    }
    lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
    non_blank <- lines[nzchar(trimws(lines))]
    if (length(non_blank) < 2) {
        log_warn(sprintf("Section 8: %s CSV has no data rows (only blank/header): %s", label, path))
        return(NULL)
    }
    tryCatch(
        read.csv(path, stringsAsFactors = FALSE),
        error = function(e) {
            log_warn(sprintf("Section 8: Failed to read %s CSV (%s): %s", label, path, conditionMessage(e)))
            NULL
        }
    )
}

if (!has_qc_summary && !has_dbl_summary) {
    log_warn("Section 8: Neither QC nor doublet summary CSVs found on disk — skipping integrated summary.")
} else {

    # ── 8.1 Load per-stage summaries ─────────────────────────────────────────

    # QC summary
    if (has_qc_summary) {
        qc_sum <- .safe_read_csv(qc_csv_path, "QC summary")
        if (!is.null(qc_sum) && "sample" %in% colnames(qc_sum)) {
            qc_sum <- qc_sum[qc_sum$sample != "TOTAL", ]   # drop totals row
        } else {
            qc_sum <- NULL
        }
    } else {
        qc_sum <- NULL
    }
    if (is.null(qc_sum)) {
        qc_sum <- data.frame(sample = character(0), cells_raw = integer(0), cells_filtered = integer(0))
        log_warn("Section 8: QC summary CSV missing/unreadable — cells_raw and cells_post_qc will be NA.")
    }

    # Doublet summary
    if (has_dbl_summary) {
        dbl_sum <- .safe_read_csv(dbl_csv_path, "Doublet summary")
        if (!is.null(dbl_sum) && "DBR_Source" %in% colnames(dbl_sum)) {
            dbl_sum <- dbl_sum[dbl_sum$DBR_Source != "ERROR", ]   # drop failed samples
        }
    } else {
        dbl_sum <- NULL
    }
    if (is.null(dbl_sum)) {
        dbl_sum <- data.frame(Sample = character(0), Total_Cells = integer(0), Doublets_Found = integer(0))
        log_warn("Section 8: Doublet summary CSV missing/unreadable — doublet columns will be NA.")
    }

    # ── 8.2 Build unified table ───────────────────────────────────────────────
    reverse_mode <- .as_bool(cfg$reverse_mode %||% "TRUE", TRUE)

    if (reverse_mode) {
        # Reverse mode sequence requested by user:
        # raw -> after ghost filter -> DBR used -> after doublet -> %doublet -> after QC -> %QC -> %overall

        # Raw counts: try cached CSV first, fall back to reading RDS files
        rcc_cache_path <- file.path(DOUBLET_SUMMARY_DIR, "raw_cell_counts_cache.csv")
        raw_tbl <- data.frame(sample = character(0), cells_raw = integer(0), stringsAsFactors = FALSE)
        if (file.exists(rcc_cache_path)) {
            rcc_cached <- tryCatch(read.csv(rcc_cache_path, stringsAsFactors = FALSE), error = function(e) NULL)
            if (!is.null(rcc_cached) && all(c("sample", "cells_raw") %in% colnames(rcc_cached))) {
                raw_tbl <- data.frame(
                    sample = .canon_sample_name(rcc_cached$sample),
                    cells_raw = as.integer(rcc_cached$cells_raw),
                    stringsAsFactors = FALSE
                )
                log_info(sprintf("Section 8: Read raw cell counts from cache (%d samples).", nrow(raw_tbl)))
            }
        }
        if (nrow(raw_tbl) == 0) {
            raw_map <- discover_rds_files(rds_dir, pattern = rds_pattern, recursive = recursive_discovery)
            if (length(raw_map) > 0) {
                raw_tbl <- do.call(rbind, lapply(names(raw_map), function(nm) {
                    p <- raw_map[[nm]]
                    n_raw <- tryCatch({
                        obj_tmp <- readRDS(p)
                        n <- as.integer(ncol(obj_tmp))
                        rm(obj_tmp); gc()
                        n
                    }, error = function(e) {
                        log_warn(sprintf("Section 8: Could not read raw cell count for %s (%s): %s", nm, p, conditionMessage(e)))
                        NA_integer_
                    })
                    data.frame(sample = .canon_sample_name(nm), cells_raw = n_raw, stringsAsFactors = FALSE)
                }))
                raw_tbl <- aggregate(cells_raw ~ sample, data = raw_tbl, FUN = function(v) max(v, na.rm = TRUE))
            }
        }

        # Doublet summary fields
        dbl_tbl <- data.frame(sample = character(0),
                              cells_after_ghost = numeric(0),
                              dbr_used = numeric(0),
                              doublets_found = numeric(0),
                              pct_doublet = numeric(0),
                              stringsAsFactors = FALSE)
        if (nrow(dbl_sum) > 0) {
            dbl_tbl <- data.frame(
                sample = .canon_sample_name(dbl_sum$Sample),
                cells_after_ghost = suppressWarnings(as.numeric(dbl_sum$Total_Cells)),
                dbr_used = if ("DBR_Used" %in% colnames(dbl_sum)) suppressWarnings(as.numeric(dbl_sum$DBR_Used)) else NA_real_,
                doublets_found = suppressWarnings(as.numeric(dbl_sum$Doublets_Found)),
                pct_doublet = suppressWarnings(as.numeric(dbl_sum$Percent_Doublet)),
                stringsAsFactors = FALSE
            )
        }

        # QC summary fields (QC is run after doublet in reverse mode)
        qc_tbl <- data.frame(sample = character(0), cells_after_qc = numeric(0), stringsAsFactors = FALSE)
        if (nrow(qc_sum) > 0 && all(c("sample", "cells_filtered") %in% colnames(qc_sum))) {
            qc_tbl <- data.frame(
                sample = .canon_sample_name(qc_sum$sample),
                cells_after_qc = suppressWarnings(as.numeric(qc_sum$cells_filtered)),
                stringsAsFactors = FALSE
            )
        }

        int_df <- merge(raw_tbl, dbl_tbl, by = "sample", all = TRUE)
        int_df <- merge(int_df, qc_tbl, by = "sample", all = TRUE)

        int_df$cells_post_qc <- int_df$cells_after_qc
        int_df$cells_post_doublet <- int_df$cells_after_ghost - int_df$doublets_found
        int_df$cells_after_doublet <- int_df$cells_post_doublet

        # If Percent_Doublet is missing, compute from ghost-filtered denominator.
        miss_dbl_pct <- is.na(int_df$pct_doublet) & !is.na(int_df$doublets_found) &
                        !is.na(int_df$cells_after_ghost) & int_df$cells_after_ghost > 0
        int_df$pct_doublet[miss_dbl_pct] <- round(int_df$doublets_found[miss_dbl_pct] /
                                                  int_df$cells_after_ghost[miss_dbl_pct] * 100, 2)

        int_df$cells_removed_by_doublet <- int_df$doublets_found
        int_df$pct_removed_by_doublet <- int_df$pct_doublet

        int_df$cells_removed_by_qc <- int_df$cells_after_doublet - int_df$cells_after_qc
        int_df$pct_removed_by_qc <- ifelse(!is.na(int_df$cells_after_doublet) & int_df$cells_after_doublet > 0,
                                           round(int_df$cells_removed_by_qc / int_df$cells_after_doublet * 100, 2),
                                           NA_real_)

        int_df$cells_removed_total <- int_df$cells_raw - int_df$cells_after_qc
        int_df$pct_removed_total <- ifelse(!is.na(int_df$cells_raw) & int_df$cells_raw > 0,
                                           round(int_df$cells_removed_total / int_df$cells_raw * 100, 2),
                                           NA_real_)
        int_df$pct_retained_final <- ifelse(!is.na(int_df$pct_removed_total),
                                            100 - int_df$pct_removed_total,
                                            NA_real_)

        log_info("Section 8: Reverse-mode integrated sequence applied: raw -> ghost filter -> doublet -> QC")
    } else {
        # Normal order: raw -> QC -> doublet
        if (has_qc_summary && nrow(qc_sum) > 0) {
            int_df <- qc_sum[, c("sample", "cells_raw", "cells_filtered")]
            colnames(int_df) <- c("sample", "cells_raw", "cells_post_qc")
        } else {
            int_df <- data.frame(
                sample       = dbl_sum$Sample,
                cells_raw    = NA_integer_,
                cells_post_qc = dbl_sum$Total_Cells
            )
        }

        if (has_dbl_summary && nrow(dbl_sum) > 0) {
            dbl_join <- dbl_sum[, c("Sample", "Total_Cells", "Doublets_Found", "Percent_Doublet")]
            colnames(dbl_join) <- c("sample", "cells_post_qc_dbl", "doublets_found", "pct_doublet")
            int_df <- merge(int_df, dbl_join, by = "sample", all = TRUE)
            int_df$cells_post_qc <- ifelse(!is.na(int_df$cells_post_qc), int_df$cells_post_qc, int_df$cells_post_qc_dbl)
            int_df$cells_post_qc_dbl <- NULL
            int_df$cells_post_doublet <- int_df$cells_post_qc - int_df$doublets_found
        } else {
            int_df$doublets_found     <- NA_integer_
            int_df$pct_doublet        <- NA_real_
            int_df$cells_post_doublet <- NA_integer_
        }

        fill_idx <- which(is.na(int_df$cells_raw) & !is.na(int_df$cells_post_qc))
        if (length(fill_idx) > 0) {
            int_df$cells_raw[fill_idx] <- int_df$cells_post_qc[fill_idx]
            log_warn(sprintf("Section 8: Filled missing cells_raw with cells_post_qc for %d sample(s)", length(fill_idx)))
        }
        int_df$cells_removed_by_qc      <- int_df$cells_raw      - int_df$cells_post_qc
        int_df$cells_removed_by_doublet <- int_df$cells_post_qc  - int_df$cells_post_doublet
        int_df$cells_removed_total       <- int_df$cells_raw      - int_df$cells_post_doublet
        int_df$pct_removed_by_qc      <- ifelse(!is.na(int_df$cells_raw) & int_df$cells_raw > 0,
                                                round(int_df$cells_removed_by_qc / int_df$cells_raw * 100, 2), NA_real_)
        int_df$pct_removed_by_doublet <- ifelse(!is.na(int_df$cells_raw) & int_df$cells_raw > 0,
                                                round(int_df$cells_removed_by_doublet / int_df$cells_raw * 100, 2), NA_real_)
        int_df$pct_removed_total      <- ifelse(!is.na(int_df$cells_raw) & int_df$cells_raw > 0,
                                                round(int_df$cells_removed_total / int_df$cells_raw * 100, 2), NA_real_)
        int_df$pct_retained_final     <- ifelse(!is.na(int_df$pct_removed_total), 100 - int_df$pct_removed_total, NA_real_)
    }

    # Keep explicit reverse-mode columns in output order requested by user
    if (!("cells_after_ghost" %in% colnames(int_df))) int_df$cells_after_ghost <- NA_real_
    if (!("dbr_used" %in% colnames(int_df)))          int_df$dbr_used <- NA_real_
    if (!("cells_after_doublet" %in% colnames(int_df))) int_df$cells_after_doublet <- int_df$cells_post_doublet
    if (!("cells_after_qc" %in% colnames(int_df)))      int_df$cells_after_qc <- int_df$cells_post_qc

    # Shorten sample names for plots
    int_df$sample_short <- int_df$sample
    int_df$sample_short <- gsub("d10_1016_j_",             "", int_df$sample_short)
    int_df$sample_short <- gsub("d10_1126_sciadv_",        "", int_df$sample_short)
    int_df$sample_short <- gsub("d10_1038_",               "", int_df$sample_short)
    int_df$sample_short <- gsub("dno_doi_kidney_organoid_", "kidney_", int_df$sample_short)

    # Add TOTAL row
    total_row <- data.frame(
        sample                   = "TOTAL",
        cells_raw                = sum(int_df$cells_raw,                na.rm = TRUE),
        cells_after_ghost        = sum(int_df$cells_after_ghost,        na.rm = TRUE),
        dbr_used                 = NA_real_,
        cells_after_doublet      = sum(int_df$cells_after_doublet,      na.rm = TRUE),
        cells_after_qc           = sum(int_df$cells_after_qc,           na.rm = TRUE),
        cells_post_qc            = sum(int_df$cells_post_qc,           na.rm = TRUE),
        cells_post_doublet       = sum(int_df$cells_post_doublet,      na.rm = TRUE),
        doublets_found           = sum(int_df$doublets_found,          na.rm = TRUE),
        pct_doublet              = NA_real_,
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
    total_row$pct_doublet <- ifelse(total_row$cells_after_ghost > 0,
                                    round(total_row$cells_removed_by_doublet / total_row$cells_after_ghost * 100, 2),
                                    NA_real_)
    total_row$pct_removed_by_qc <- ifelse(total_row$cells_after_doublet > 0,
                                          round(total_row$cells_removed_by_qc / total_row$cells_after_doublet * 100, 2),
                                          NA_real_)
    total_row$pct_removed_by_doublet <- total_row$pct_doublet
    total_row$pct_removed_total <- ifelse(total_row$cells_raw > 0,
                                          round(total_row$cells_removed_total / total_row$cells_raw * 100, 2),
                                          NA_real_)
    total_row$pct_retained_final <- ifelse(!is.na(total_row$pct_removed_total),
                                           100 - total_row$pct_removed_total,
                                           NA_real_)

    # Ensure row-bind is schema-safe even if upstream branches produced extra/missing columns
    miss_in_total <- setdiff(colnames(int_df), colnames(total_row))
    if (length(miss_in_total) > 0) {
        for (cn in miss_in_total) total_row[[cn]] <- NA
    }
    miss_in_int <- setdiff(colnames(total_row), colnames(int_df))
    if (length(miss_in_int) > 0) {
        for (cn in miss_in_int) int_df[[cn]] <- NA
    }
    total_row <- total_row[, colnames(int_df), drop = FALSE]

    int_df_full <- rbind(int_df, total_row)

    # ── 8.3 Save integrated CSV ───────────────────────────────────────────────
    dir.create(INTEGRATED_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

    int_csv_path <- file.path(INTEGRATED_SUMMARY_DIR, "integrated_qc_doublet_summary.csv")

    # Clean CSV view: one row per sample + TOTAL, with fixed output schema
    .first_non_na <- function(v) {
        vv <- v[!is.na(v)]
        if (length(vv) == 0) return(NA_real_)
        vv[[1]]
    }

    csv_rows <- lapply(unique(int_df_full$sample), function(s) {
        d <- int_df_full[int_df_full$sample == s, , drop = FALSE]
        get_num <- function(col, mode = "first") {
            if (!(col %in% colnames(d))) return(NA_real_)
            v <- suppressWarnings(as.numeric(d[[col]]))
            v <- v[!is.na(v)]
            if (length(v) == 0) return(NA_real_)
            if (mode == "max") return(max(v))
            if (mode == "sum") return(sum(v))
            v[[1]]
        }
        data.frame(
            sample = as.character(s),
            cells_raw = get_num("cells_raw", "max"),
            cells_after_ghost = get_num("cells_after_ghost", "max"),
            dbr_used = get_num("dbr_used", "first"),
            doublets_found = get_num("doublets_found", "max"),
            pct_doublets = get_num("pct_doublet", "first"),
            cells_after_doublet_removal = get_num("cells_after_doublet", "max"),
            cells_after_qc = get_num("cells_after_qc", "max"),
            pct_qc_removal = get_num("pct_removed_by_qc", "first"),
            total_cells_removal = get_num("cells_removed_total", "max"),
            pct_total_cells_removal = get_num("pct_removed_total", "first"),
            stringsAsFactors = FALSE
        )
    })

    int_csv_df <- do.call(rbind, csv_rows)
    int_csv_df <- int_csv_df[int_csv_df$sample != "TOTAL", , drop = FALSE]

    if (nrow(int_csv_df) > 0) {
        total_row_clean <- data.frame(
            sample = "TOTAL",
            cells_raw = sum(int_csv_df$cells_raw, na.rm = TRUE),
            cells_after_ghost = sum(int_csv_df$cells_after_ghost, na.rm = TRUE),
            dbr_used = NA_real_,
            doublets_found = sum(int_csv_df$doublets_found, na.rm = TRUE),
            pct_doublets = NA_real_,
            cells_after_doublet_removal = sum(int_csv_df$cells_after_doublet_removal, na.rm = TRUE),
            cells_after_qc = sum(int_csv_df$cells_after_qc, na.rm = TRUE),
            pct_qc_removal = NA_real_,
            total_cells_removal = sum(int_csv_df$total_cells_removal, na.rm = TRUE),
            pct_total_cells_removal = NA_real_,
            stringsAsFactors = FALSE
        )

        if (is.finite(total_row_clean$cells_after_ghost) && total_row_clean$cells_after_ghost > 0)
            total_row_clean$pct_doublets <- round(total_row_clean$doublets_found / total_row_clean$cells_after_ghost * 100, 2)
        if (is.finite(total_row_clean$cells_after_doublet_removal) && total_row_clean$cells_after_doublet_removal > 0)
            total_row_clean$pct_qc_removal <- round((total_row_clean$cells_after_doublet_removal - total_row_clean$cells_after_qc) / total_row_clean$cells_after_doublet_removal * 100, 2)
        if (is.finite(total_row_clean$cells_raw) && total_row_clean$cells_raw > 0)
            total_row_clean$pct_total_cells_removal <- round(total_row_clean$total_cells_removal / total_row_clean$cells_raw * 100, 2)

        int_csv_df <- rbind(int_csv_df, total_row_clean)
    }

    write.csv(int_csv_df, int_csv_path, row.names = FALSE)
    log_info(sprintf("Integrated summary CSV saved: %s", int_csv_path))

    # Print to log
    log_info("=== Integrated QC + Doublet Summary ===")
    print_int_cols <- c(
        "sample", "cells_raw", "cells_after_ghost", "dbr_used", "doublets_found",
        "pct_doublets", "cells_after_doublet_removal", "cells_after_qc",
        "pct_qc_removal", "total_cells_removal", "pct_total_cells_removal"
    )
    for (ln in capture.output(print(int_csv_df[, print_int_cols], row.names = FALSE))) log_info(ln)

    # ── 8.4 Build integrated dashboard plots ─────────────────────────────────
    # Work only with per-sample rows (exclude TOTAL) for per-sample plots
    pd <- int_df[complete.cases(int_df[, c("cells_raw", "cells_post_qc")]), ]

    if (nrow(pd) == 0) {
        log_warn("Section 8: No complete sample rows — skipping dashboard plots.")
    } else {

        # ── Panel 1: Waterfall bar chart ──────────────────────────────────────
        if (reverse_mode) {
            df_waterfall <- melt(
                pd[, c("sample_short", "cells_raw", "cells_after_ghost", "cells_after_doublet", "cells_after_qc")],
                id.vars = "sample_short"
            )
            df_waterfall$variable <- factor(df_waterfall$variable,
                levels = c("cells_raw", "cells_after_ghost", "cells_after_doublet", "cells_after_qc"),
                labels = c("Raw", "Post-Ghost", "Post-Doublet", "Post-QC"))
            wf_cols <- c("Raw" = "#bdbdbd", "Post-Ghost" = "#74add1", "Post-Doublet" = "#1a9641", "Post-QC" = "#2c7fb8")
            wf_subtitle <- "Raw -> ghost filter -> doublet removal -> QC"
        } else {
            df_waterfall <- melt(
                pd[, c("sample_short", "cells_raw", "cells_post_qc", "cells_post_doublet")],
                id.vars = "sample_short"
            )
            df_waterfall$variable <- factor(df_waterfall$variable,
                levels = c("cells_raw", "cells_post_qc", "cells_post_doublet"),
                labels = c("Raw", "Post-QC", "Post-Doublet"))
            wf_cols <- c("Raw" = "#bdbdbd", "Post-QC" = "#2c7fb8", "Post-Doublet" = "#1a9641")
            wf_subtitle <- "Raw -> Post-QC -> Post-Doublet removal"
        }

        p8_1 <- ggplot(df_waterfall, aes(x = sample_short, y = value, fill = variable)) +
            geom_bar(stat = "identity", position = "dodge", alpha = 0.85, width = 0.7) +
              scale_fill_manual(values = wf_cols) +
            scale_y_continuous(labels = scales::comma) +
            labs(title    = "Cell Counts Across All Three Stages",
                  subtitle = wf_subtitle,
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
            levels = c("pct_removed_by_doublet", "pct_removed_by_qc"),
            labels = c("Doublet Removal", "QC Removal"))

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
        if (reverse_mode) {
            df_funnel <- pd[, c("sample_short", "cells_raw", "cells_after_ghost", "cells_after_doublet", "cells_after_qc")]
            df_funnel_long <- melt(df_funnel, id.vars = "sample_short")
            df_funnel_long$stage <- factor(df_funnel_long$variable,
                levels = c("cells_raw", "cells_after_ghost", "cells_after_doublet", "cells_after_qc"),
                labels = c("Raw", "Post-Ghost", "Post-Doublet", "Post-QC"))
        } else {
            df_funnel <- pd[, c("sample_short", "cells_raw", "cells_post_qc", "cells_post_doublet")]
            df_funnel_long <- melt(df_funnel, id.vars = "sample_short")
            df_funnel_long$stage <- factor(df_funnel_long$variable,
                levels = c("cells_raw", "cells_post_qc", "cells_post_doublet"),
                labels = c("Raw", "Post-QC", "Post-Doublet"))
        }

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
        #   Doublet-filtered (singlets) → DOUBLET_DIR
        density_list <- list()

        for (samp in pd$sample) {
            samp_short <- pd$sample_short[pd$sample == samp]

            # Raw
            raw_path <- if (exists("rds_file_map")) rds_file_map[[samp]] else NULL
            # QC-filtered
            qc_path  <- file.path(FILTERED_DIR, paste0(samp, "_filtered.rds"))
            # Doublet-filtered singlets (try two common naming patterns)
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
            dbl_full   <- file.path(DOUBLET_DIR, paste0(samp, "_with_doublet_calls.rds"))
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
        total_retained <- if (reverse_mode) sum(pd$cells_after_qc, na.rm = TRUE) else sum(pd$cells_post_doublet, na.rm = TRUE)
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
                    ifelse(
                        sum(pd$cells_raw, na.rm = TRUE) > 0,
                        sum(pd$cells_post_doublet, na.rm = TRUE) / sum(pd$cells_raw, na.rm = TRUE) * 100,
                        NA_real_
                    )
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
log_info(sprintf("  Integrated CSV    : %s", file.path(INTEGRATED_SUMMARY_DIR, "integrated_qc_doublet_summary.csv")))
log_info(sprintf("  Integrated PNG    : %s", file.path(INTEGRATED_SUMMARY_DIR, "integrated_summary_dashboard.png")))
log_info(sprintf("  Integrated PDF    : %s", file.path(INTEGRATED_SUMMARY_DIR, "integrated_summary_dashboard.pdf")))
log_info(sprintf("  Log file          : %s", log_file_path))
log_info("══════════════════════════════════════════════════")

} # end sys.nframe() guard — only runs when script is executed directly
log_info("══════════════════════════════════════════════════")