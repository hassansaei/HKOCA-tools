# =============================================================================
# Single-Cell QC Pipeline — Full Integration
# Doublet Detection (scDblFinder) → QC Filtering
#
# 3 stages: (1) Doublet — runs scDblFinder, saves singlet-only RDS + audit PDF
# + summary CSV. (2) QC — applies threshold filters (nFeature/nCount/% mito),
# saves filtered RDS + audit PDF + dashboard. (3) Summary — combined CSV/plots
# across both stages.
#
# Config is a DCF file (qc_config.dcf) — see qc_config.dcf for all keys
# (paths, doublet tuning, per-sample platform/chemistry/dbr overrides, etc).
#
# Usage:
#   Rscript QC_scdbl_Combined.R [--config PATH] [--run all|qc|doublet] [OPTIONS]
#
# Run with --help for the full list of CLI flags.
# =============================================================================

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 0 — Package Loading
# ══════════════════════════════════════════════════════════════════════════════

load_packages_safely <- function() {
    required_packages <- c(
        "Seurat", "MASS", "ggplot2", "patchwork", "cowplot", 
        "scales", "reshape2", "ggridges", "ggrepel", "SingleCellExperiment", 
        "scDblFinder", "yaml", "dplyr", "jsonlite", "digest", "anndata", "reticulate"
    )

    missing <- character(0)
    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            missing <- c(missing, pkg)
        }
    }

    if (length(missing) > 0) {
        stop(sprintf("Missing required packages: %s\nPlease construct the environment using 'conda env create -f conda/environment_qc.yaml' and activate it.", paste(missing, collapse = ", ")))
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
        library(anndata)    
        library(reticulate) 
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
# Lenient DCF parser — handles multiline values and embedded CSV blocks.
.read_config_dcf <- function(config_path) {
    if (!file.exists(config_path))
        stop(sprintf("Config file not found: %s", config_path))

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

# ── 1.7 QC decisions loader ───────────────────────────────────────────────────
.load_qc_decisions <- function(cfg) {
    qc_decisions <- NULL
    source_label <- NULL

    # 1. Check if the user passed an external CSV via --filters
    if (!is.null(cfg$filters) && nzchar(trimws(cfg$filters))) {
        filter_file <- trimws(cfg$filters)
        if (!file.exists(filter_file)) {
            stop(sprintf("Filter CSV file not found: %s", filter_file))
        }
        qc_decisions <- read.csv(
            file = filter_file,
            stringsAsFactors = FALSE,
            colClasses = "character",
            strip.white = TRUE
        )
        source_label <- sprintf("file: %s", filter_file)
        
    # 2. Check for the inline DCF table
    } else if (!is.null(cfg$qc_decisions_table) && nzchar(trimws(cfg$qc_decisions_table))) {
        inline_table <- cfg$qc_decisions_table
        lines <- unlist(strsplit(gsub("\r\n?", "\n", inline_table), "\n", fixed = TRUE))
        lines <- lines[nzchar(trimws(lines))]

        if (length(lines) >= 2) {
            qc_decisions <- read.csv(
                text = paste(lines, collapse = "\n"),
                stringsAsFactors = FALSE,
                colClasses = "character",
                strip.white = TRUE
            )
            source_label <- "config: qc_decisions_table (inline)"
        }
    }

    # 3. FALLBACK: Dynamically generate a MAD3 table across discovered files
    if (is.null(qc_decisions)) {
        log_warn("No filter criteria provided. Auto-generating global 'mad3' default criteria for discovered datasets.")
        
        # Resolve inputs to discover files matching what the pipeline will process
        in_dir <- if (!is.null(cfg$doublet_subdir) && nzchar(cfg$doublet_subdir)) {
            # If running standard 'qc' stage standalone, it targets doublet dir
            file.path(cfg$output_dir %||% getwd(), cfg$doublet_subdir)
        } else {
            cfg$rds_dir %||% getwd()
        }
        
        # Safe detection of available files using existing logic pattern
        pat <- cfg$rds_pattern %||% "\\.rds$"
        rec <- .as_bool(cfg$recursive_discovery %||% "FALSE", FALSE)
        
        discovered_paths <- tryCatch({
            list.files(in_dir, pattern = pat, full.names = TRUE, recursive = rec, ignore.case = TRUE)
        }, error = function(e) character(0))
        
        if (length(discovered_paths) == 0) {
            stop("No filters provided and no input RDS files discovered to apply defaults to.")
        }
        
        basenames <- gsub(pat, "", basename(discovered_paths), ignore.case = TRUE)
        
        # Build a standard data frame matching the expected schema
        qc_decisions <- data.frame(
            Dataset_Name         = basenames,
            Lower_Feature_Method = rep("mad3", length(basenames)),
            Upper_Feature_Method = rep("mad3", length(basenames)),
            Lower_Count_Method   = rep("mad3", length(basenames)),
            Upper_Count_Method   = rep("mad3", length(basenames)),
            Upper_Mito_Method    = rep("mad3", length(basenames)),
            stringsAsFactors     = FALSE
        )
        source_label <- "pipeline_default: dynamic MAD3 auto-generation"
    }

    if (!("Dataset_Name" %in% colnames(qc_decisions))) {
        stop("QC decisions table must contain a 'Dataset_Name' column.")
    }

    required_cols <- c("Dataset_Name", "Lower_Feature_Method", "Upper_Feature_Method",
                       "Lower_Count_Method", "Upper_Count_Method", "Upper_Mito_Method")
    missing_cols  <- setdiff(required_cols, colnames(qc_decisions))
    if (length(missing_cols) > 0) {
        stop(sprintf("QC decisions table is missing required column(s): %s",
                     paste(missing_cols, collapse = ", ")))
    }

    list(data = qc_decisions, source = source_label)
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

# Harmonize writes <output_root>/<study>/rds/*_harmonized.rds — search recursively
# from output_root. Stage output dirs (doublet_filtered_rds, etc.) are flat.
.is_harmonized_rds_root <- function(input_dir) {
    identical(
        normalizePath(input_dir, mustWork = FALSE),
        normalizePath(rds_dir, mustWork = FALSE)
    )
}

.discover_input_rds <- function(input_dir) {
    if (.is_harmonized_rds_root(input_dir)) {
        log_info(sprintf(
            "  RDS discovery: recursive=%s, pattern=%s",
            recursive_discovery, rds_pattern
        ))
        discover_rds_files(input_dir, pattern = rds_pattern, recursive = recursive_discovery)
    } else {
        discover_rds_files(input_dir, pattern = "\\.rds$", recursive = FALSE)
    }
}

# Strip common harmonize / QC suffixes for consistent sample lookup.
.canon_sample_name <- function(x) {
    s <- tolower(trimws(as.character(x)))
    s <- gsub("\\.rds$", "", s, ignore.case = TRUE)
    s <- gsub(
        "_harmonized_filtered$|_harmonized$|_filtered$|_with_doublet_calls$|_singlets$",
        "",
        s,
        ignore.case = TRUE
    )
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
    cat(sprintf("  Rscript %s [--config PATH] [--run all|qc|doublet] [OPTIONS]\n\n", script_name))
    cat("Options:\n")
    cat(sprintf("  --config PATH            Config DCF file (default: %s)\n", default_config))
    cat("  --filters PATH           Path to an external CSV file for QC cutoffs\n")
    cat("  --run STR                all | qc | doublet  (default: all)\n")
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
# Tracks params in a hidden .meta.json so we skip re-running unless they change.
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

# ── 1.12 10x multiplet-rate tables ───────────────────────────────────────────
# Approximate published values; DBR is interpolated piecewise-linearly from
# post-ghost-filter cell count (with edge-slope extrapolation beyond range).
#
# v2/v3 chemistry:
DBL_10X_V3_TABLE_LOADED    <- c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000)
DBL_10X_V3_TABLE_RECOVERED <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
DBL_10X_V3_TABLE_DBR       <- c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076)
#
# v4 chemistry (CG000731 Rev A — ~0.58x v3 multiplet rate):
DBL_10X_V4_TABLE_LOADED    <- c(725, 1450, 2900, 4350, 5800, 7250, 8700, 10150, 11600, 13050,
                                 14500, 15950, 17400, 18850, 20300, 21750, 23200, 24650, 26100, 27550, 29000)
DBL_10X_V4_TABLE_RECOVERED <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                                 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000)
DBL_10X_V4_TABLE_DBR       <- c(0.002, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.036,
                                 0.040, 0.044, 0.048, 0.052, 0.056, 0.060, 0.064, 0.068, 0.072, 0.076, 0.080)
#
# 5' chemistry (same curve as v4):
DBL_10X_5P_TABLE_LOADED    <- c(725, 1450, 2900, 4350, 5800, 7250, 8700, 10150, 11600, 13050,
                                 14500, 15950, 17400, 18850, 20300, 21750, 23200, 24650, 26100, 27550, 29000)
DBL_10X_5P_TABLE_RECOVERED <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                                 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000)
DBL_10X_5P_TABLE_DBR       <- c(0.002, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.036,
                                 0.040, 0.044, 0.048, 0.052, 0.056, 0.060, 0.064, 0.068, 0.072, 0.076, 0.080)

# Backward-compatible aliases
DBL_10X_TABLE_LOADED    <- DBL_10X_V3_TABLE_LOADED
DBL_10X_TABLE_RECOVERED <- DBL_10X_V3_TABLE_RECOVERED
DBL_10X_TABLE_DBR       <- DBL_10X_V3_TABLE_DBR

# Piecewise-linear DBR estimator with edge-slope extrapolation outside range.
.estimate_dbr_from_table <- function(cell_count, tbl_counts, tbl_dbr) {
    x <- .as_num(cell_count, NA_real_)
    if (is.na(x) || x <= 0) return(NA_real_)

    if (x >= min(tbl_counts) && x <= max(tbl_counts)) {
        return(as.numeric(stats::approx(
            x = tbl_counts,
            y = tbl_dbr,
            xout = x,
            method = "linear",
            rule = 2
        )$y))
    }

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

    as.numeric(min(max(d, 0), 0.30))
}

.estimate_10x_dbr_from_count <- function(cell_count) {
    .estimate_dbr_from_table(cell_count, DBL_10X_V3_TABLE_RECOVERED, DBL_10X_V3_TABLE_DBR)
}

.estimate_10x_v4_dbr_from_count <- function(cell_count) {
    .estimate_dbr_from_table(cell_count, DBL_10X_V4_TABLE_RECOVERED, DBL_10X_V4_TABLE_DBR)
}

.estimate_10x_5p_dbr_from_count <- function(cell_count) {
    .estimate_dbr_from_table(cell_count, DBL_10X_5P_TABLE_RECOVERED, DBL_10X_5P_TABLE_DBR)
}

# Dispatch DBR estimator by chemistry string (v2/v3, v4, 5p).
.estimate_10x_dbr_by_chemistry <- function(cell_count, chemistry = "v3") {
    chem <- tolower(trimws(chemistry))
    if (chem == "v4") {
        .estimate_10x_v4_dbr_from_count(cell_count)
    } else if (chem %in% c("5p", "5prime", "5'")) {
        .estimate_10x_5p_dbr_from_count(cell_count)
    } else {
        .estimate_10x_dbr_from_count(cell_count)
    }
}

# Parse sc_protocol metadata string into normalised platform + chemistry.
# e.g. "10x_3_v3.1" → (10x, v3), "10X_5" → (10x, 5p), "Dropseq" → (dropseq, NA)
.parse_sc_protocol <- function(protocol_string) {
    s <- tolower(trimws(protocol_string))

    if (grepl("drop.?seq", s)) {
        return(list(platform = "dropseq", chemistry = NA_character_))
    }

    if (grepl("pipseq", s)) {
        return(list(platform = "pipseq", chemistry = NA_character_))
    }

    if (grepl("^10x", s)) {
        if (grepl("_5($|[_'])", s)) {
            return(list(platform = "10x", chemistry = "5p"))
        }

        ver_match <- regmatches(s, regexpr("v[0-9]+(\\.[0-9]+)?", s))
        if (length(ver_match) == 1 && nzchar(ver_match)) {
            chem <- ver_match
            if (chem == "v3.1") chem <- "v3"
            if (chem == "v1")   chem <- "v2"
            return(list(platform = "10x", chemistry = chem))
        }

        return(list(platform = "10x", chemistry = "v3"))
    }

    list(platform = s, chemistry = NA_character_)
}

# ── 1.13 Small-sample guards for doublet stage ──────────────────────────────
.make_all_singlet_calls <- function(cell_names, score_value = 0) {
    n_cells <- length(cell_names)
    list(
        class = setNames(rep("singlet", n_cells), cell_names),
        score = setNames(rep(score_value, n_cells), cell_names)
    )
}

.check_scdblfinder_input <- function(obj_or_sce, min_cells = NULL) {
    if (is.null(min_cells)) min_cells <- DBL_MIN_CELLS

    cell_count <- ncol(obj_or_sce)
    feature_count <- nrow(obj_or_sce)

    if (is.na(cell_count) || cell_count == 0)
        return(list(ok = FALSE, reason = "no cells remain after ghost-cell filtering"))
    if (cell_count < min_cells)
        return(list(ok = FALSE, reason = sprintf("too few cells for scDblFinder (%d < %d)", cell_count, min_cells)))
    if (is.na(feature_count) || feature_count < 2)
        return(list(ok = FALSE, reason = sprintf("too few features for scDblFinder (%d < 2)", feature_count)))

    list(ok = TRUE, reason = NULL)
}

.run_optional_doublet_reductions <- function(obj, requested_dims) {
    result <- list(obj = obj, pca_dims = integer(0), umap_ready = FALSE)

    cell_count <- ncol(obj)
    feature_count <- nrow(obj)
    if (cell_count < 3 || feature_count < 2) {
        log_warn(sprintf("  Skipping PCA/UMAP: too few cells/features (%d cells, %d features).",
                         cell_count, feature_count))
        return(result)
    }

    requested_dims <- sort(unique(as.integer(requested_dims[is.finite(requested_dims)])))
    requested_dims <- requested_dims[requested_dims >= 1L]
    if (length(requested_dims) == 0) {
        log_warn("  Skipping PCA/UMAP: no positive UMAP dimensions were requested.")
        return(result)
    }

    tryCatch({
        obj <- NormalizeData(obj, verbose = FALSE)
        obj <- FindVariableFeatures(obj, verbose = FALSE)

        var_features <- VariableFeatures(obj)
        if (length(var_features) < 2) {
            log_warn(sprintf("  Skipping PCA/UMAP: only %d variable feature(s) identified.",
                             length(var_features)))
            result$obj <- obj
            return(result)
        }

        obj <- ScaleData(obj, features = var_features, verbose = FALSE)

        max_pcs <- min(length(var_features), cell_count - 1L, feature_count - 1L, max(requested_dims))
        if (!is.finite(max_pcs) || max_pcs < 2) {
            log_warn(sprintf("  Skipping PCA/UMAP: only %d PC(s) are supportable for this sample.",
                             as.integer(max_pcs)))
            result$obj <- obj
            return(result)
        }

        obj <- RunPCA(obj, features = var_features, npcs = max_pcs, verbose = FALSE)

        available_pcs <- ncol(Embeddings(obj, "pca"))
        usable_dims <- requested_dims[requested_dims <= available_pcs]

        if (length(usable_dims) < length(requested_dims)) {
            log_warn(sprintf("  Clamping UMAP dims to available PCs: requested [%s], using [%s].",
                             paste(requested_dims, collapse = ","),
                             paste(usable_dims, collapse = ",")))
        }

        result$obj <- obj
        result$pca_dims <- usable_dims

        if (length(usable_dims) < 2) {
            log_warn(sprintf("  Skipping UMAP: only %d usable PC dimension(s) available.",
                             length(usable_dims)))
            return(result)
        }

        if (cell_count <= length(usable_dims)) {
            log_warn(sprintf("  Skipping UMAP: %d cells are insufficient for %d requested dimensions.",
                             cell_count, length(usable_dims)))
            return(result)
        }

        umap_neighbors <- min(30L, cell_count - 1L)
        if (umap_neighbors < 2L) {
            log_warn(sprintf("  Skipping UMAP: only %d neighbor(s) are supportable for %d cells.",
                             umap_neighbors, cell_count))
            return(result)
        }

        obj <- RunUMAP(obj, dims = usable_dims, n.neighbors = umap_neighbors, verbose = FALSE)
        result$obj <- obj
        result$umap_ready <- "umap" %in% names(obj@reductions)

        if (!result$umap_ready)
            log_warn("  UMAP reduction missing after RunUMAP; continuing without UMAP plot.")

        result
    }, error = function(e) {
        log_warn(sprintf("  Visualization reductions skipped: %s", conditionMessage(e)))
        result$obj <- obj
        result
    })
}

.run_scdblfinder_safe <- function(sce, label, dbr = NULL, samples = NULL) {
    viability <- .check_scdblfinder_input(sce)
    if (!viability$ok) {
        log_warn(sprintf("  %s: %s — marking all %d cells as singlets.",
                         label, viability$reason, ncol(sce)))
        singlet_calls <- .make_all_singlet_calls(colnames(sce))
        sce$scDblFinder.class <- singlet_calls$class
        sce$scDblFinder.score <- singlet_calls$score
        return(list(sce = sce, fallback = TRUE, reason = viability$reason))
    }

    tryCatch({
        sce_out <- if (!is.null(samples)) {
            if (is.null(dbr)) scDblFinder(sce, samples = samples)
            else              scDblFinder(sce, samples = samples, dbr = dbr)
        } else {
            if (is.null(dbr)) scDblFinder(sce)
            else              scDblFinder(sce, dbr = dbr)
        }
        list(sce = sce_out, fallback = FALSE, reason = NULL)
    }, error = function(e) {
        log_warn(sprintf("  %s: scDblFinder failed (%s) — marking all %d cells as singlets.",
                         label, conditionMessage(e), ncol(sce)))
        singlet_calls <- .make_all_singlet_calls(colnames(sce))
        sce$scDblFinder.class <- singlet_calls$class
        sce$scDblFinder.score <- singlet_calls$score
        list(sce = sce, fallback = TRUE, reason = conditionMessage(e))
    })
}

# ── 1.14 Data-driven ambiguous zone detection ─────────────────────────────────
# Finds where singlet/doublet KDE curves overlap (each normalised to its own
# peak). Returns list(lo, hi, method); (NA, NA, "NONE") if no overlap.
.compute_overlap_zone <- function(scores, classes,
                                  min_overlap_frac = 0.05,
                                  kde_n = 512L) {
    s_scores <- scores[!is.na(scores) & !is.na(classes) & classes == "singlet"]
    d_scores <- scores[!is.na(scores) & !is.na(classes) & classes == "doublet"]

    no_overlap <- list(lo = NA_real_, hi = NA_real_, method = "NONE")

    if (length(s_scores) < 5 || length(d_scores) < 5)
        return(no_overlap)

    grid <- seq(0, 1, length.out = kde_n)
    bw_s <- tryCatch(bw.nrd0(s_scores), error = function(e) 0.05)
    bw_d <- tryCatch(bw.nrd0(d_scores), error = function(e) 0.05)
    bw_s <- max(bw_s, 1e-4); bw_d <- max(bw_d, 1e-4)

    kde_s <- density(s_scores, bw = bw_s, from = 0, to = 1, n = kde_n)$y
    kde_d <- density(d_scores, bw = bw_d, from = 0, to = 1, n = kde_n)$y

    kde_s_norm <- kde_s / max(kde_s, na.rm = TRUE)
    kde_d_norm <- kde_d / max(kde_d, na.rm = TRUE)

    in_overlap <- kde_s_norm >= min_overlap_frac & kde_d_norm >= min_overlap_frac

    if (!any(in_overlap))
        return(no_overlap)

    overlap_idx <- which(in_overlap)
    lo <- grid[min(overlap_idx)]
    hi <- grid[max(overlap_idx)]

    if ((hi - lo) < 0.01)
        return(no_overlap)

    list(lo = lo, hi = hi,
         method = sprintf("KDE_OVERLAP(min_frac=%.3f,n=%d)", min_overlap_frac, kde_n))
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Configuration & Path Setup
# ══════════════════════════════════════════════════════════════════════════════

if (sys.nframe() == 0L || !exists(".read_config_dcf")) {
    script_dir     <- .get_script_dir()
    default_config <- file.path(script_dir, "qc_config.dcf")
    cli_args       <- .parse_cli_args(commandArgs(trailingOnly = TRUE))

    if (.as_bool(cli_args$help, FALSE)) { .print_usage(default_config); quit(save = "no", status = 0) }

config_path <- if (!is.null(cli_args$config)) .resolve_path(cli_args$config, getwd()) else default_config
cfg         <- .read_config_dcf(config_path)

# CLI overrides win over config file
for (nm in setdiff(names(cli_args), c("config", "help", "run", "stage"))) cfg[[nm]] <- cli_args[[nm]]

if (!is.null(cli_args$stage) && (is.null(cli_args$run) || !nzchar(cli_args$run)))
    message("Note: --stage is deprecated; use --run (all | qc | doublet)")

QC_RUN <- tolower(cli_args$run %||% cli_args$stage %||% "all")
if (!QC_RUN %in% c("all", "qc", "doublet"))
    stop(sprintf("Invalid --run '%s'. Choose: all | qc | doublet", QC_RUN))

RUN_QC      <- QC_RUN %in% c("all", "qc")
RUN_DOUBLET <- QC_RUN %in% c("all", "doublet")

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

RUN_H5AD        <- .as_bool(cfg$run_h5ad_conversion %||% "FALSE", FALSE)
h5ad_subdir     <- if (!is.null(cfg$h5ad_output_subdir) && nzchar(cfg$h5ad_output_subdir)) cfg$h5ad_output_subdir else "h5ad_converted"
H5AD_DIR        <- file.path(BASE_OUT_DIR, h5ad_subdir)

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
DBL_MIN_CELLS_RUN  <- as.integer(cfg$dbl_min_cells_run %||% 25)
DBL_MIN_CELLS      <- DBL_MIN_CELLS_RUN

# ── Ambiguous doublet zone ────────────────────────────────────────────────────
# Spans where singlet/doublet KDE curves overlap, per sample. Set
# dbl_ambiguous_min_overlap = 1.0 for pure binary mode (no ambiguous cells).
DBL_AMBIGUOUS_MIN_OVERLAP <- .as_num(cfg$dbl_ambiguous_min_overlap %||% "0.05", 0.05)
DBL_AMBIGUOUS_KDE_N       <- as.integer(cfg$dbl_ambiguous_kde_n     %||% "512")

# ── Per-sample DBR config ──────────────────────────────────────────────────────
# From DCF keys like sample.<name>.platform / .chemistry / .dbr
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
# DCF keys: platform_dbr.<platform>: <fraction>
DBL_PLATFORM_DBR  <- list()
pdbr_keys         <- grep("^platform_dbr\\.", names(cfg), value = TRUE)
for (pk in pdbr_keys) {
    pname <- tolower(sub("^platform_dbr\\.", "", pk))
    DBL_PLATFORM_DBR[[pname]] <- as.numeric(cfg[[pk]])
}

# ── Plotting parameters ──────────────────────────────────────────────────────
PLOT_DPI    <- .as_num(cfg$plot_dpi %||% "300", 300)
PLOT_BG     <- trimws(cfg$plot_bg %||% "white")
PLOT_WIDTH  <- .as_num(cfg$plot_base_width %||% "6", 6)
PLOT_HEIGHT <- .as_num(cfg$plot_base_height %||% "5", 5)

# ── Discovery settings ────────────────────────────────────────────────────────
rds_pattern         <- cfg$rds_pattern         %||% "\\.rds$"
recursive_discovery <- .as_bool(cfg$recursive_discovery %||% "FALSE", FALSE)

# ── Create output directories ─────────────────────────────────────────────────
for (d in c(QC_DIR, FILTERED_DIR, DOUBLET_DIR, DOUBLET_SUMMARY_DIR, INTEGRATED_SUMMARY_DIR, dirname(log_file_path)))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (RUN_H5AD) dir.create(H5AD_DIR, recursive = TRUE, showWarnings = FALSE)
# ── Initialise shared logger ──────────────────────────────────────────────────
.init_logger(log_file_path)
run_ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

log_info("══════════════════════════════════════════════════")
log_info(sprintf("QC PIPELINE START  [run: %s]", QC_RUN))
log_info(sprintf("  Config       : %s", config_path))
log_info(sprintf("  Raw RDS dir  : %s", rds_dir))
decisions_src_display <- "pipeline_default: dynamic MAD3 auto-generation"
if (!is.null(cfg$filters) && nzchar(trimws(cfg$filters))) {
    decisions_src_display <- sprintf("file: %s", trimws(cfg$filters))
} else if (!is.null(cfg$qc_decisions_table) && nzchar(trimws(cfg$qc_decisions_table))) {
    decisions_src_display <- "config: qc_decisions_table (inline)"
}
log_info(sprintf("  Decisions src: %s", decisions_src_display))
log_info(sprintf("  Output base  : %s", BASE_OUT_DIR))
log_info(sprintf("  Log file     : %s", log_file_path))
log_info("══════════════════════════════════════════════════")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Threshold Algorithm Implementations
# ══════════════════════════════════════════════════════════════════════════════

# ── Triangle threshold ────────────────────────────────────────────────────────
# Picks the bin farthest from the peak-to-endpoint line.
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

# ── Dispatcher: method string → numeric threshold ─────────────────────────────
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

# 2-D kernel density for scatter point colouring.
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

# Splits into one violin per group_col level (e.g. sample_id) when given.
make_violin <- function(obj, feat, title, max_y, cutoff_l = NULL, cutoff_h = NULL,
                        group_col = NULL) {
    has_group <- !is.null(group_col) && nzchar(group_col) &&
                 group_col %in% colnames(obj@meta.data) &&
                 length(unique(na.omit(obj@meta.data[[group_col]]))) > 1

    if (has_group) {
        grp <- as.character(obj@meta.data[[group_col]])
        grp[is.na(grp) | !nzchar(trimws(grp))] <- "unassigned"
        df  <- data.frame(x = grp, y = as.numeric(obj@meta.data[[feat]]),
                          stringsAsFactors = FALSE)
        df  <- df[is.finite(df$y), , drop = FALSE]
        df$x <- factor(df$x, levels = sort(unique(df$x)))

        n_grp <- nlevels(df$x)
        x_lbl <- group_col
        title <- sprintf("%s  (n_%s = %d)", title, group_col, n_grp)

        p <- ggplot(df, aes(x = x, y = y, fill = x)) +
            geom_violin(color = NA, alpha = 0.85, scale = "width", trim = FALSE) +
            geom_jitter(width = 0.2, size = 0.1, alpha = 0.10, color = "grey25") +
            scale_y_continuous(limits = c(0, max_y), oob = scales::squish) +
            scale_fill_viridis_d(option = "turbo", end = 0.92) +
            labs(x = x_lbl, y = feat, title = title) +
            theme_bw() +
            theme(
                plot.background  = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = "grey80"),
                axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
                axis.text.y      = element_text(size = 7),
                legend.position  = "none",
                plot.title       = element_text(size = 9, face = "bold")
            )
    } else {
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
    }

    if (!is.null(cutoff_l) && is.finite(cutoff_l))
        p <- p + geom_hline(yintercept = cutoff_l, linewidth = 0.5, color = "red")
    if (!is.null(cutoff_h) && is.finite(cutoff_h))
        p <- p + geom_hline(yintercept = cutoff_h, linewidth = 0.5, color = "red")
    p
}

# Density-coloured scatter split into one panel per group.
make_scatter_facet <- function(x, y, group, x_label, y_label, title,
                               max_x, max_y,
                               cutoff_l = NULL, cutoff_h = NULL,
                               show_cutoffs = TRUE) {
    grp <- as.character(group)
    grp[is.na(grp) | !nzchar(trimws(grp))] <- "unassigned"
    df_all <- data.frame(X1 = x, X2 = y, sample = grp,
                         stringsAsFactors = FALSE)
    df_all <- df_all[is.finite(df_all$X1) & is.finite(df_all$X2), , drop = FALSE]
    if (nrow(df_all) == 0)
        return(ggplot() + theme_void() + ggtitle(title))

    df_all <- do.call(rbind, lapply(split(df_all, df_all$sample), function(d) {
        d$dens <- point_density(d$X1, d$X2)
        d
    }))
    finite_dens <- df_all$dens[is.finite(df_all$dens)]
    dens_limit  <- if (length(finite_dens) == 0 || max(finite_dens) == 0) 1 else {
        dl <- quantile(finite_dens, 0.99)
        if (!is.finite(dl) || dl <= 0) max(finite_dens) else dl
    }

    n_grp <- length(unique(df_all$sample))
    ncol_facet <- min(4, max(1, n_grp))

    use_raster <- nrow(df_all) > 50000 && requireNamespace("ggrastr", quietly = TRUE)
    geom_pts   <- if (use_raster) ggrastr::geom_point_rast(size = 0.3, raster.dpi = 300)
                  else            geom_point(size = 0.3)

    p <- ggplot(df_all, aes(x = X1, y = X2, color = dens)) +
        geom_pts +
        facet_wrap(~sample, ncol = ncol_facet) +
        scale_x_continuous(limits = c(0, max_x), oob = scales::squish) +
        scale_y_continuous(limits = c(0, max_y), oob = scales::squish) +
                scale_color_viridis_c(option = "magma", na.value = "grey70",
                              limits = c(0, dens_limit), oob = scales::squish) +
        labs(x = x_label, y = y_label, color = "Density") +
        ggtitle(sprintf("%s  (per-sample)", title)) +

        theme_bw() +
        theme(
            plot.background  = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = "grey80"),
            panel.grid.major = element_line(color = "grey92"),
            panel.grid.minor = element_blank(),
            plot.title       = element_text(size = 9, face = "bold"),
            strip.text       = element_text(size = 7, face = "bold"),
            axis.title       = element_text(size = 8),
            axis.text        = element_text(size = 7),
            legend.position  = "bottom",
            legend.title     = element_text(size = 7, face = "bold"),
            legend.text      = element_text(size = 6)
        ) +
        guides(color = guide_colourbar(barwidth = 8, barheight = 0.5))

    if (show_cutoffs && !is.null(cutoff_l) && is.finite(cutoff_l))
        p <- p + geom_hline(yintercept = cutoff_l, linewidth = 0.5, color = "red", linetype = "dashed")
    if (show_cutoffs && !is.null(cutoff_h) && is.finite(cutoff_h))
        p <- p + geom_hline(yintercept = cutoff_h, linewidth = 0.5, color = "red", linetype = "dashed")
    p
}

# Returns the first metadata column from a candidate list that has >1 unique value.
.detect_qc_group_col <- function(meta) {
    candidates <- c("sample_id", "orig.ident", "donor", "donor_id",
                    "library_id", "Sample", "sample")
    for (col in candidates) {
        if (!(col %in% colnames(meta))) next
        vals <- as.character(meta[[col]])
        vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
        if (length(unique(vals)) > 1) return(col)
    }
    NULL
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

    rds_file_map <- .discover_input_rds(qc_input_dir)
    sample_names <- sort(names(rds_file_map))
    log_info(sprintf("Discovered %d RDS files.", length(sample_names)))
    log_info(sprintf("Datasets: %s", paste(sample_names, collapse = ", ")))

    decisions_payload <- .load_qc_decisions(cfg)
    qc_decisions <- decisions_payload$data
    log_info(sprintf("Loaded QC decisions from %s", decisions_payload$source))

    required_vars <- c("rds_dir", "QC_DIR", "FILTERED_DIR",
                       "s_median", "s_mean", "s_max")
    missing_vars  <- required_vars[!sapply(required_vars, function(v) exists(v, inherits = TRUE))]
    if (length(missing_vars) > 0)
        stop(sprintf("Missing required variables: %s", paste(missing_vars, collapse = ", ")))
    if (!is.list(rds_file_map) || length(rds_file_map) == 0)
        stop(sprintf("No RDS files discovered in QC input directory: %s", qc_input_dir))

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

    qc_basenames <- gsub("\\.rds$", "", basename(qc_decisions$Dataset_Name), ignore.case = TRUE)
    configured_datasets <- intersect(qc_basenames, sample_names)
    if (length(configured_datasets) < length(qc_basenames)) {
        unmatched <- setdiff(qc_basenames, sample_names)
        singlets_matches <- intersect(paste0(unmatched, "_singlets"), sample_names)
        configured_datasets <- union(configured_datasets, singlets_matches)
    }

    if (length(configured_datasets) == 0)
        log_warn("No configured datasets overlap with discovered RDS files.")

    log_info(sprintf("Loaded %d QC decision rows. Running %d dataset(s).",
                     nrow(qc_decisions), length(configured_datasets)))

    qc_pdf_path <- file.path(QC_DIR, "QC_Before_After_Report.pdf")
    pdf(qc_pdf_path, width = 14, height = 6)
    on.exit(tryCatch(dev.off(), error = function(e) NULL), add = TRUE)

    summary_list <- list()

    for (i in seq_len(nrow(qc_decisions))) {

        nm       <- gsub("\\.rds$", "", basename(qc_decisions$Dataset_Name[i]), ignore.case = TRUE)
        rds_path <- rds_file_map[[nm]]
        if (is.null(rds_path) && identical(normalizePath(qc_input_dir, mustWork = FALSE),
                                           normalizePath(DOUBLET_DIR, mustWork = FALSE))) {
            nm_singlets <- sub("_filtered$", "_singlets", nm)
            rds_path <- rds_file_map[[nm_singlets]]
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

        filtered_path <- file.path(FILTERED_DIR, paste0(nm, "_filtered.rds"))
        current_qc_params <- as.list(row_cfg)

        if (.should_skip_run(filtered_path, current_qc_params, cfg$force_overwrite)) {
            log_info(sprintf("[QC SKIPPED] %s already processed with identical thresholds.", nm))
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

            # Seurat v4/v5 count extraction: try @counts, then @layers, then GetAssayData().
            obj_raw <- readRDS(rds_path)
            log_info(sprintf("  Raw class: %s | Assay class: %s",
                             class(obj_raw), class(obj_raw[["RNA"]])))

            assay_obj <- obj_raw[["RNA"]]
            cm <- tryCatch(slot(assay_obj, "counts"), error = function(e) NULL)
            if (is.null(cm) || (is(cm, "Matrix") && length(cm@x) == 0) ||
                (is.matrix(cm) && length(cm) == 0)) {
                if (.hasSlot(assay_obj, "layers")) {
                    layers <- slot(assay_obj, "layers")
                    layer_nm <- if ("counts" %in% names(layers)) "counts"
                                else grep("^counts(\\.|$)", names(layers), value = TRUE)[1]
                    if (!is.null(layer_nm) && !is.na(layer_nm)) {
                        cm <- layers[[layer_nm]]
                        # v5 split layers can be missing dimnames; realign to full barcode set.
                        if (is.null(rownames(cm)) || is.null(colnames(cm))) {
                            feats <- slot(assay_obj, "features")
                            cells <- slot(assay_obj, "cells")
                            rownames(cm) <- rownames(feats)[feats[, layer_nm]]
                            colnames(cm) <- rownames(cells)[cells[, layer_nm]]
                        }
                    }
                }
            }
            if (is.null(cm)) {
                cm <- tryCatch(
                    GetAssayData(obj_raw, assay = "RNA", layer = "counts"),
                    error = function(e) GetAssayData(obj_raw, assay = "RNA", slot = "counts")
                )
            }

            meta_in  <- as.data.frame(obj_raw@meta.data)
            obj      <- CreateSeuratObject(counts = cm, meta.data = meta_in)
            rm(obj_raw); gc()

            if (is.null(obj$nCount_RNA)   || all(obj$nCount_RNA   == 0))
                obj$nCount_RNA   <- Matrix::colSums(cm)
            if (is.null(obj$nFeature_RNA) || all(obj$nFeature_RNA == 0))
                obj$nFeature_RNA <- Matrix::colSums(cm > 0)

            gene_names   <- rownames(obj)
            mito_genes   <- grep("^MT-", gene_names, value = TRUE, ignore.case = TRUE)
            ribo_genes   <- grep("^(RPL|RPS)", gene_names, value = TRUE, perl = TRUE, ignore.case = TRUE)

            total_counts <- as.numeric(obj$nCount_RNA)

            obj$percent.mito <- if (length(mito_genes) > 0) {
                as.numeric(Matrix::colSums(cm[mito_genes, , drop = FALSE]) / total_counts * 100)
            } else {
                rep(0, ncol(obj))
            }

            obj$percent.ribo <- if (length(ribo_genes) > 0) {
                as.numeric(Matrix::colSums(cm[ribo_genes, , drop = FALSE]) / total_counts * 100)
            } else {
                rep(0, ncol(obj))
            }

            nCount_vec   <- as.numeric(obj@meta.data[["nCount_RNA"]])
            nFeature_vec <- as.numeric(obj@meta.data[["nFeature_RNA"]])
            nonzero_keep <- colnames(obj)[nCount_vec > 0 & nFeature_vec > 0]
            if (length(nonzero_keep) < ncol(obj)) obj <- subset(obj, cells = nonzero_keep)

            meta_raw <- obj[[]]

            feat_lower  <- get_thresh(row_cfg$Lower_Feature_Method, meta_raw$nFeature_RNA, "lower")
            feat_upper  <- get_thresh(row_cfg$Upper_Feature_Method, meta_raw$nFeature_RNA, "upper")
            count_lower <- get_thresh(row_cfg$Lower_Count_Method,   meta_raw$nCount_RNA,   "lower")
            count_upper <- get_thresh(row_cfg$Upper_Count_Method,   meta_raw$nCount_RNA,   "upper")
            mito_upper  <- get_thresh(row_cfg$Upper_Mito_Method,    meta_raw$percent.mito, "upper")

            ribo_vals  <- meta_raw$percent.ribo
            ribo_lower <- if (median(ribo_vals, na.rm = TRUE) >= 10) quantile(ribo_vals, 0.05, na.rm = TRUE) else 0
            ribo_upper <- max(ribo_vals, na.rm = TRUE)

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

            # Raw axis maxima set the scale so pre/post plots are directly comparable.
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

            # Auto-detect a per-sample grouping column; emit additional per-group
            # violins and faceted scatters when >1 sample is present.
            qc_group_col <- .detect_qc_group_col(obj@meta.data)
            n_qc_groups  <- if (!is.null(qc_group_col))
                                length(unique(na.omit(obj@meta.data[[qc_group_col]])))
                            else 1L
            emit_per_sample <- !is.null(qc_group_col) && n_qc_groups > 1
            violin_width_grp <- min(max(7, 1.0 * n_qc_groups), 18)
            scatter_facet_w  <- min(4, n_qc_groups) * 3.2
            scatter_facet_h  <- ceiling(n_qc_groups / 4) * 3.0 + 1

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

                ggsave(file.path(QC_DIR, nm, paste0("Violin_qc.raw.",      feat, ".png")), p.v,   width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = PLOT_BG)
                ggsave(file.path(QC_DIR, nm, paste0("Plot_qc.raw.",        feat, ".png")), p.s,   width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = PLOT_BG)
                ggsave(file.path(QC_DIR, nm, paste0("Violin_qc.filtered.", feat, ".png")), p.v.f, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = PLOT_BG)
                ggsave(file.path(QC_DIR, nm, paste0("Plot_qc.filtered.",   feat, ".png")), p.s.f, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = PLOT_BG)

                if (emit_per_sample) {

                    grp_r  <- obj@meta.data[[qc_group_col]]
                    grp_f  <- fobj@meta.data[[qc_group_col]]

                    p.v.bs   <- make_violin(obj,  feat, nm, max.y, cutoff.l, cutoff.h,
                                            group_col = qc_group_col)
                    p.s.bs   <- make_scatter_facet(x_vals, y_vals, grp_r, x_lbl, feat, nm,
                                                   max.ct, max.y, cutoff.l, cutoff.h,
                                                   show_cutoffs = TRUE)
                    p.v.f.bs <- make_violin(fobj, feat, nm, max.y,
                                            group_col = qc_group_col)
                    p.s.f.bs <- make_scatter_facet(xf, yf, grp_f, x_lbl, feat, nm,
                                                   max.ct, max.y, NULL, NULL,
                                                   show_cutoffs = FALSE)

                    ggsave(file.path(QC_DIR, nm, paste0("Violin_qc.raw.by_sample.",      feat, ".png")),
                           p.v.bs,   width = violin_width_grp, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = PLOT_BG)
                    ggsave(file.path(QC_DIR, nm, paste0("Plot_qc.raw.by_sample.",        feat, ".png")),
                           p.s.bs,   width = scatter_facet_w,  height = scatter_facet_h, dpi = PLOT_DPI, bg = PLOT_BG)
                    ggsave(file.path(QC_DIR, nm, paste0("Violin_qc.filtered.by_sample.", feat, ".png")),
                           p.v.f.bs, width = violin_width_grp, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = PLOT_BG)
                    ggsave(file.path(QC_DIR, nm, paste0("Plot_qc.filtered.by_sample.",   feat, ".png")),
                           p.s.f.bs, width = scatter_facet_w,  height = scatter_facet_h, dpi = PLOT_DPI, bg = PLOT_BG)
                }           
            }
            log_info(sprintf("  PNGs: %d files → %s%s",
                             if (emit_per_sample) 32L else 16L,
                             file.path(QC_DIR, nm),
                             if (emit_per_sample)
                                 sprintf("  (per-sample plots split by '%s', %d samples)",
                                         qc_group_col, n_qc_groups)
                             else
                                 "  (single-sample dataset — no per-sample split)"))

            p1_b <- plot_dist(meta_raw, "nFeature_RNA", paste("Raw nFeature:", nm), feat_lower,  feat_upper,  TRUE,  "lightblue")
            p2_b <- plot_dist(meta_raw, "nCount_RNA",   paste("Raw nCount:",   nm), count_lower, count_upper, TRUE,  "lightgreen")
            p3_b <- plot_dist(meta_raw, "percent.mito", paste("Raw Mito %:",   nm), -Inf,        mito_upper,  FALSE, "lightcoral")
            print((p1_b | p2_b | p3_b) + plot_annotation(
                title    = paste("BEFORE QC:", nm),
                subtitle = sprintf("Lower Feat: %s | Upper Feat: %s | Upper Mito: %s",
                                   row_cfg$Lower_Feature_Method,
                                   row_cfg$Upper_Feature_Method,
                                   row_cfg$Upper_Mito_Method)))

            meta_filt <- fobj[[]]
            p1_a <- plot_dist(meta_filt, "nFeature_RNA", paste("Filtered nFeature:", nm), -Inf, Inf, TRUE,  "lightblue")
            p2_a <- plot_dist(meta_filt, "nCount_RNA",   paste("Filtered nCount:",   nm), -Inf, Inf, TRUE,  "lightgreen")
            p3_a <- plot_dist(meta_filt, "percent.mito", paste("Filtered Mito %:",   nm), -Inf, Inf, FALSE, "lightcoral")
            print((p1_a | p2_a | p3_a) + plot_annotation(
                title = paste("AFTER QC:", nm, sprintf("(%.1f%% cells removed)", pct_removed))))

            filtered_path <- file.path(FILTERED_DIR, paste0(nm, "_filtered.rds"))
            saveRDS(fobj, filtered_path)
            .save_run_metadata(filtered_path, current_qc_params, extra = list(
                cells_raw = cells_raw, cells_filtered = cells_filtered
            ))

            inflation      <- round(s_median(fobj$nCount_RNA) / s_median(obj$nCount_RNA), 2)
            inflation_flag <- if (inflation > 1.5)   " \u26a0 OVER-FILTERING"
                              else if (inflation < 1.05) " \u26a0 CHECK: very light filtering"
                              else ""
            log_info(sprintf("  Saved  : %s", filtered_path))
            log_info(sprintf("  Cells  : %d \u2192 %d  (%.1f%% removed)", cells_raw, cells_filtered, pct_removed))
            log_info(sprintf("  nCount : median %g \u2192 %g  [inflation %.2fx]%s",
                             s_median(obj$nCount_RNA), s_median(fobj$nCount_RNA),
                             inflation, inflation_flag))

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

    # ── Full QC Dashboard ─────────────────────────────────────────────────────
    if (nrow(summary_table) == 0) {
        log_warn("Summary table empty — skipping dashboard.")
    } else {

        plot_data <- summary_table[summary_table$sample != "TOTAL", ]

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

        # P4: nCount median shift
        df_ncount      <- plot_data[, c("sample_short", "median_nCount_raw", "median_nCount_filtered")]
        df_long_ncount <- melt(df_ncount, id.vars = "sample_short")
        p4 <- ggplot(df_long_ncount, aes(x = sample_short, y = value, group = variable, color = variable)) +
            geom_line(linewidth = 0.8) + geom_point(size = 2) +
            scale_color_manual(values = c("median_nCount_raw" = "#fee0d2", "median_nCount_filtered" = "#de2d26"),
                               labels = c("Before", "After")) +
            labs(title = "Signal Recovery: Median nCount_RNA", y = "Counts", x = "", color = "Stage") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")

        # P5: nFeature median shift
        df_nfeat      <- plot_data[, c("sample_short", "median_nFeature_raw", "median_nFeature_filtered")]
        df_long_nfeat <- melt(df_nfeat, id.vars = "sample_short")
        p5 <- ggplot(df_long_nfeat, aes(x = sample_short, y = value, group = variable, color = variable)) +
            geom_line(linewidth = 0.8) + geom_point(size = 2) +
            scale_color_manual(values = c("median_nFeature_raw" = "#e0f3db", "median_nFeature_filtered" = "#43a2ca"),
                               labels = c("Before", "After")) +
            labs(title = "Gene Coverage: Median nFeature_RNA", y = "Genes", x = "", color = "Stage") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "top")

        # P6: Mito/ribo contamination reduction
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

            global_max_nCount <- max(dist_df$log1p_nCount, na.rm = TRUE)
            global_max_nFeat  <- max(dist_df$log1p_nFeat,  na.rm = TRUE)

            p_counts_ridge <- ggplot(dist_df, aes(x = log1p_nCount, y = sample, fill = sample)) +
                geom_density_ridges(alpha = 0.7, scale = 1.5, color = "white", linewidth = 0.4) +
                geom_vline(xintercept = global_max_nCount, color = "black", linetype = "dashed", linewidth = 0.6) +
                scale_x_continuous(limits = c(0, global_max_nCount * 1.05)) +
                facet_wrap(~status, ncol = 2) +
                theme_minimal() +
                labs(title = sprintf("Vertical Across Datasets: log1p(nCount) [Global Max: %.2f]", global_max_nCount),
                     x = "log1p(nCount_RNA)", y = "") +
                theme(axis.text.y = element_text(size = 8), legend.position = "none")

            p_genes_ridge <- ggplot(dist_df, aes(x = log1p_nFeat, y = sample, fill = sample)) +
                geom_density_ridges(alpha = 0.7, scale = 1.5, color = "white", linewidth = 0.4) +
                geom_vline(xintercept = global_max_nFeat, color = "black", linetype = "dashed", linewidth = 0.6) +
                scale_x_continuous(limits = c(0, global_max_nFeat * 1.05)) +
                facet_wrap(~status, ncol = 2) +
                theme_minimal() +
                labs(title = sprintf("Vertical Across Datasets: log1p(nFeature) [Global Max: %.2f]", global_max_nFeat),
                     x = "log1p(nFeature_RNA)", y = "") +
                theme(axis.text.y = element_text(size = 8), legend.position = "none")

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
               plot = final_layout, width = 18, height = 42, dpi = PLOT_DPI, bg = PLOT_BG)
        log_info(sprintf("Full dashboard saved: %s", file.path(QC_DIR, "summary_qc_full_dashboard.png")))
    }

}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — STAGE 1: Doublet Detection (scDblFinder)
# ══════════════════════════════════════════════════════════════════════════════

.run_stage_doublet <- function(doublet_input_dir = FILTERED_DIR) {

    log_info("──────────────────────────────────────────────────")
    log_info("STAGE 1: DOUBLET DETECTION")
    log_info(sprintf("  Input  : %s", doublet_input_dir))
    log_info(sprintf("  Output : %s", DOUBLET_DIR))
    log_info("──────────────────────────────────────────────────")

    rds_map <- .discover_input_rds(doublet_input_dir)
    filtered_files <- unname(unlist(rds_map))
    if (length(filtered_files) == 0)
        stop(sprintf("No input RDS files found in: %s", doublet_input_dir))

    log_info(sprintf("Found %d input RDS files.", length(filtered_files)))

    if (length(DBL_SAMPLE_CFG) > 0) {
        log_info(sprintf("Per-sample DBR overrides loaded for %d sample(s): %s",
                         length(DBL_SAMPLE_CFG),
                         paste(names(DBL_SAMPLE_CFG), collapse = ", ")))
    } else {
        log_info(sprintf("No per-sample DBR overrides in config. Default platform: %s.",
                         DBL_DEFAULT_PLATFORM))
    }

    dbl_pdf <- file.path(DOUBLET_SUMMARY_DIR, paste0("Doublet_Audit_Report_", run_ts, ".pdf"))
    pdf(dbl_pdf, width = 16, height = 10)
    on.exit(tryCatch(dev.off(), error = function(e) NULL), add = TRUE)

    dbl_summary_stats <- list()
    raw_cell_counts   <- list()

    for (file_path in filtered_files) {

        sample_nm <- gsub(
            "_harmonized_filtered\\.rds$|_filtered\\.rds$|\\.rds$",
            "",
            basename(file_path)
        )
        log_info(sprintf("[DOUBLET] Processing: %s", sample_nm))

        singlets_path <- file.path(DOUBLET_DIR, paste0(sample_nm, "_singlets.rds"))

        skip_dbl_list <- .parse_csv_values(cfg$skip_scdblfinder_samples %||% "")
        is_explicitly_skipped <- .canon_sample_name(sample_nm) %in% .canon_sample_name(skip_dbl_list)

        dbl_run_params <- list(
            platform = DBL_DEFAULT_PLATFORM,
            chemistry = DBL_DEFAULT_CHEMISTRY,
            min_count = DBL_MIN_COUNT,
            min_feature = DBL_MIN_FEATURE,
            min_cells_run = DBL_MIN_CELLS,
            sample_cfg = DBL_SAMPLE_CFG[[.canon_sample_name(sample_nm)]],
            skipped_via_config = is_explicitly_skipped
        )

        # Skip if params unchanged; falls back to file-existence check pre-upgrade.
        skip_doublet <- .should_skip_run(singlets_path, dbl_run_params, cfg$force_overwrite)
        if (!skip_doublet && !.as_bool(cfg$force_overwrite, FALSE) &&
            file.exists(singlets_path) && file.info(singlets_path)$size > 0 &&
            !file.exists(.get_meta_path(singlets_path))) {
            skip_doublet <- TRUE
            log_info(sprintf("  [COMPAT] No .meta.json for %s — skipping via legacy file-existence check.", sample_nm))
        }
        if (skip_doublet) {
            log_info(sprintf("[DOUBLET SKIPPED] %s already processed with identical parameters.", sample_nm))
            dbl_stats <- tryCatch({
                dbl_full_path <- file.path(DOUBLET_DIR, paste0(sample_nm, "_with_doublet_calls.rds"))
                if (file.exists(dbl_full_path)) {
                    obj_full <- readRDS(dbl_full_path)
                    n_total <- ncol(obj_full)
                    if ("doublet_status" %in% colnames(obj_full@meta.data)) {
                        ds <- obj_full$doublet_status
                        n_dbl <- sum(ds == "doublet",   na.rm = TRUE)
                        n_amb <- sum(ds == "ambiguous",  na.rm = TRUE)
                        n_sng <- sum(ds == "singlet",    na.rm = TRUE)
                    } else {
                        n_dbl <- sum(obj_full$scDblFinder.class == "doublet", na.rm = TRUE)
                        n_amb <- NA_integer_
                        n_sng <- n_total - n_dbl
                    }
                    pct_dbl <- round((n_dbl / n_total) * 100, 2)
                    pct_amb <- if (!is.na(n_amb)) round((n_amb / n_total) * 100, 2) else NA_real_
                    rm(obj_full); gc()
                    list(total = n_total, doublets = n_dbl, ambiguous = n_amb,
                         singlets = n_sng, pct = pct_dbl, pct_amb = pct_amb)
                } else {
                    list(total = NA, doublets = NA, ambiguous = NA,
                         singlets = NA, pct = NA, pct_amb = NA)
                }
            }, error = function(e) list(total = NA, doublets = NA, ambiguous = NA,
                                        singlets = NA, pct = NA, pct_amb = NA))
            dbl_summary_stats[[sample_nm]] <- data.frame(
                Sample              = sample_nm,
                Total_Cells         = dbl_stats$total,
                Doublets_Found      = dbl_stats$doublets,
                Ambiguous_Found     = dbl_stats$ambiguous,
                Singlets_Retained   = dbl_stats$singlets,
                Cells_in_Singlets_RDS = (dbl_stats$singlets %||% 0) + (dbl_stats$ambiguous %||% 0),
                Percent_Doublet     = dbl_stats$pct,
                Percent_Ambiguous   = dbl_stats$pct_amb,
                Ambiguous_Score_Lo  = NA_real_,
                Ambiguous_Score_Hi  = NA_real_,
                DBR_Used            = NA_real_,
                Cells_Loaded_Used   = NA_real_,
                Cells_Recovered_Used = NA_real_,
                DBR_Source          = "SKIPPED_RESUME",
                stringsAsFactors    = FALSE
            )
            next 
        }

        dbl_ok <- tryCatch({
            obj <- readRDS(file_path)

        # Ghost-cell filter: very low-count cells (CellBender artifacts) break
        # scDblFinder size-factor estimation.
        n_before <- ncol(obj)
        obj      <- subset(obj, subset = nCount_RNA > DBL_MIN_COUNT & nFeature_RNA > DBL_MIN_FEATURE)
        n_after_ghost <- ncol(obj)
        log_info(sprintf("  Ghost-cell filter: %d \u2192 %d cells", n_before, n_after_ghost))
        raw_cell_counts[[sample_nm]] <- n_before

        if (n_after_ghost == 0) {
            log_warn(sprintf("  No cells remain after ghost-cell filtering for '%s' — skipping sample.", sample_nm))
            dbl_summary_stats[[sample_nm]] <- data.frame(
                Sample              = sample_nm,
                Total_Cells         = 0,
                Doublets_Found      = 0,
                Ambiguous_Found     = 0,
                Singlets_Retained   = 0,
                Cells_in_Singlets_RDS = 0,
                Percent_Doublet     = 0,
                Percent_Ambiguous   = 0,
                Ambiguous_Score_Lo  = NA_real_,
                Ambiguous_Score_Hi  = NA_real_,
                DBR_Used            = NA_real_,
                Cells_Loaded_Used   = n_before,
                Cells_Recovered_Used = 0,
                DBR_Source          = "SKIPPED_EMPTY_AFTER_GHOST",
                stringsAsFactors    = FALSE
            )
            rm(obj); gc()
            TRUE
        } else {

        # Force meta.data to base data.frame (Seurat v5 S4 compatibility).
        obj@meta.data <- as.data.frame(obj@meta.data)

        viability <- .check_scdblfinder_input(obj)

        if (is_explicitly_skipped) {
            log_info(sprintf("  [SKIP scDblFinder] %s explicitly configured to bypass doublet calculation.", sample_nm))
            viability <- list(ok = FALSE, reason = "explicitly skipped via config 'skip_scdblfinder_samples'")
        }

        run_doublet_detection <- viability$ok
        used_singlet_fallback <- FALSE
        fallback_reason <- NULL

        if (!run_doublet_detection) {
            log_warn(sprintf("  %s — scDblFinder will be skipped and all remaining cells will be marked as singlets.",
                             viability$reason))
        }

        batch_levels_meta <- character(0)
        if (DBL_BATCH_COL %in% colnames(obj@meta.data)) {
            batch_levels_meta <- unique(as.character(na.omit(obj@meta.data[[DBL_BATCH_COL]])))
            batch_levels_meta <- batch_levels_meta[nzchar(trimws(batch_levels_meta))]
        }
        has_batch <- length(batch_levels_meta) > 1

        # ── Resolve platform, chemistry, and DBR ──────────────────────────────
        # Priority: per-sample config > sc_protocol metadata > global defaults.
        # 10x DBR comes from the chemistry multiplet table; non-10x uses
        # DBL_PLATFORM_DBR. Missing values never stop the run.
        dbr               <- NULL
        dbr_source        <- if (run_doublet_detection) "AUTO_10X" else "SKIPPED_TOO_FEW_CELLS"
        platform          <- DBL_DEFAULT_PLATFORM
        chemistry         <- DBL_DEFAULT_CHEMISTRY
        cells_loaded_used <- NA_real_
        cells_recovered_used <- NA_real_

        # (A) Auto-detect from sc_protocol metadata
        proto_raw <- NULL
        cell_protocols <- NULL
        if ("sc_protocol" %in% colnames(obj@meta.data)) {
            cell_protocols <- as.character(obj@meta.data$sc_protocol)
            proto_vals <- unique(na.omit(cell_protocols))

            if (length(proto_vals) == 1L) {
                proto_raw <- as.character(proto_vals)
            } else if (length(proto_vals) > 1L) {
                proto_raw <- "MULTIPLE"
                log_info(sprintf(
                    "  Multiple sc_protocol values for '%s': %s — will resolve dynamically per batch.",
                    sample_nm, paste(proto_vals, collapse = ", ")))

                # If no batching is currently configured, force it to split by sc_protocol
                # so scDblFinder evaluates v2 and v3 cell pools independently.
                if (!has_batch) {
                    log_info("  Forcing DBL_BATCH_COL = 'sc_protocol' to accurately process mixed chemistries.")
                    DBL_BATCH_COL <- "sc_protocol"
                    batch_levels_meta <- proto_vals[nzchar(trimws(proto_vals))]
                    has_batch <- length(batch_levels_meta) > 1
                }
            }
        }

        chemistry_display <- chemistry

        if (run_doublet_detection && !is.null(proto_raw) && proto_raw != "MULTIPLE") {
            parsed    <- .parse_sc_protocol(proto_raw)
            platform  <- parsed$platform
            if (!is.na(parsed$chemistry)) chemistry <- parsed$chemistry
            chemistry_display <- if (platform == "10x") chemistry else "N/A"
            dbr_source <- "SC_PROTOCOL_AUTO"
            log_info(sprintf("  sc_protocol='%s' -> platform=%s, chemistry=%s",
                             proto_raw, platform, chemistry_display))
        } else if (run_doublet_detection && !is.null(proto_raw) && proto_raw == "MULTIPLE") {
            chemistry_display <- "MIXED"
            dbr_source <- "SC_PROTOCOL_AUTO_PER_BATCH"
        } else if (run_doublet_detection) {
            log_warn(sprintf(
                "  No sc_protocol metadata for '%s' — using global defaults (platform=%s, chemistry=%s)",
                sample_nm, platform, chemistry))
        }

        # (B) Per-sample config override (wins over sc_protocol)
        sample_cfg <- DBL_SAMPLE_CFG[[.canon_sample_name(sample_nm)]]
        if (run_doublet_detection && !is.null(sample_cfg)) {
            if (!is.null(sample_cfg$platform) && nzchar(sample_cfg$platform)) {
                platform <- tolower(trimws(sample_cfg$platform))
                log_info(sprintf("  Config override: platform=%s", platform))
            }
            if (!is.null(sample_cfg$chemistry) && nzchar(sample_cfg$chemistry)) {
                chemistry <- tolower(trimws(sample_cfg$chemistry))
                log_info(sprintf("  Config override: chemistry=%s", chemistry))
            }
        }

        if (run_doublet_detection && platform == "10x" && !chemistry %in% c("v2", "v3", "v4", "5p")) {
            log_warn(sprintf(
                "  Unknown 10x chemistry '%s' for '%s' — defaulting to v3",
                chemistry, sample_nm))
            chemistry <- "v3"
        }
        chemistry_display <- if (platform == "10x") chemistry else "N/A"

        # (C) Resolve DBR
        dbr_user <- NA_real_
        if (run_doublet_detection && !is.null(sample_cfg) && !is.null(sample_cfg$dbr) && nzchar(sample_cfg$dbr)) {
            dbr_user <- suppressWarnings(as.numeric(sample_cfg$dbr))
        }

        if (run_doublet_detection && platform == "10x") {
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
        } else if (run_doublet_detection) {
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

        if (run_doublet_detection && has_batch && platform == "10x" && is.na(dbr_user)) {
            log_info(sprintf("  Platform: %s | Chemistry: %s | DBR: per-batch auto [%s]",
                             platform, chemistry_display, dbr_source))
            log_info(sprintf("  10x %s merged input before per-batch split (filtered_cells_total): %.0f",
                             chemistry, ncol(obj)))
        } else if (run_doublet_detection) {
            log_info(sprintf("  Platform: %s | Chemistry: %s | DBR: %s [%s]",
                             platform, chemistry_display,
                             ifelse(is.null(dbr), "AUTO", sprintf("%.4f", dbr)),
                             dbr_source))
        }

        if (run_doublet_detection && !is.na(cells_recovered_used) && !(has_batch && platform == "10x" && is.na(dbr_user))) {
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

        class_results <- NULL
        score_results <- NULL

        if (!run_doublet_detection) {
            singlet_calls <- .make_all_singlet_calls(colnames(obj))
            class_results <- singlet_calls$class
            score_results <- singlet_calls$score
            used_singlet_fallback <- TRUE
            fallback_reason <- viability$reason
        } else {
            log_info("  Normalizing | PCA | UMAP ...")
            reduction_state <- .run_optional_doublet_reductions(obj, DBL_UMAP_DIMS)
            obj <- reduction_state$obj

            sce <- as.SingleCellExperiment(obj)

            # Batch-aware detection prevents the ~30% doublet artifact that
            # appears when multi-sample merged objects are run as one.
            if (has_batch) {
                batch_ids <- as.character(sce[[DBL_BATCH_COL]])
                batch_ids[is.na(batch_ids)] <- ""
                batch_levels <- unique(batch_ids[nzchar(batch_ids)])
                n_sub <- length(batch_levels)
                batch_counts <- table(batch_ids[nzchar(batch_ids)])
                has_small_batches <- any(batch_counts < DBL_MIN_CELLS)

                if (platform == "10x" && is.na(dbr_user)) {
                    log_info(sprintf("  [BATCH MODE] %d internal samples — estimating DBR per batch.", n_sub))

                    class_results <- setNames(rep(NA_character_, ncol(sce)), colnames(sce))
                    score_results <- setNames(rep(NA_real_, ncol(sce)), colnames(sce))

                    batch_dbrs <- setNames(rep(NA_real_, n_sub), batch_levels)
                    batch_sizes <- setNames(rep(0L, n_sub), batch_levels)

                    for (batch_nm in batch_levels) {
                        batch_idx <- which(batch_ids == batch_nm)
                        sub_sce <- sce[, batch_idx]
                        sub_n <- ncol(sub_sce)
                        batch_sizes[[batch_nm]] <- sub_n

                        # Dynamically determine the chemistry for THIS specific batch
                        batch_chem <- chemistry # fallback to global
                        if (!is.null(cell_protocols)) {
                            batch_protos <- na.omit(cell_protocols[batch_idx])
                            if (length(batch_protos) > 0) {
                                top_proto <- names(sort(table(batch_protos), decreasing = TRUE))[1]
                                parsed_proto <- .parse_sc_protocol(top_proto)
                                if (parsed_proto$platform == "10x" && !is.na(parsed_proto$chemistry)) {
                                    batch_chem <- parsed_proto$chemistry
                                }
                            }
                        }

                        sub_dbr <- .estimate_10x_dbr_by_chemistry(sub_n, batch_chem)
                        if (!is.na(sub_dbr) && sub_dbr > 0 && sub_dbr <= 0.3) {
                            batch_dbrs[[batch_nm]] <- sub_dbr
                            log_info(sprintf("    [BATCH:%s] chemistry=%s | cells=%d | DBR=%.4f", 
                                             batch_nm, batch_chem, sub_n, sub_dbr))
                        } else {
                            log_warn(sprintf("    [BATCH:%s] could not estimate 10x %s DBR from %d filtered cells — using scDblFinder auto",
                                             batch_nm, batch_chem, sub_n))
                        }

                        tbl_max <- switch(batch_chem,
                                          v4 = max(DBL_10X_V4_TABLE_RECOVERED),
                                          `5p` = max(DBL_10X_5P_TABLE_RECOVERED),
                                          max(DBL_10X_V3_TABLE_RECOVERED))

                        if (sub_n > tbl_max) {
                            log_info(sprintf("    [BATCH:%s] extrapolation: yes (above %.0f filtered cells)",
                                             batch_nm, tbl_max))
                        }

                        sub_result <- .run_scdblfinder_safe(
                            sub_sce,
                            label = sprintf("[BATCH:%s]", batch_nm),
                            dbr = if (is.na(batch_dbrs[[batch_nm]])) NULL else batch_dbrs[[batch_nm]]
                        )
                        sub_sce <- sub_result$sce
                        used_singlet_fallback <- used_singlet_fallback || isTRUE(sub_result$fallback)
                        if (isTRUE(sub_result$fallback) && is.null(fallback_reason))
                            fallback_reason <- sprintf("batch %s: %s", batch_nm, sub_result$reason)

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
                } else if (has_small_batches) {
                    log_warn(sprintf("  [BATCH MODE] At least one internal sample has < %d cells — processing batches independently with fallback singlet calls for undersized batches.",
                                     DBL_MIN_CELLS))

                    class_results <- setNames(rep(NA_character_, ncol(sce)), colnames(sce))
                    score_results <- setNames(rep(NA_real_, ncol(sce)), colnames(sce))

                    for (batch_nm in batch_levels) {
                        batch_idx <- which(batch_ids == batch_nm)
                        sub_sce <- sce[, batch_idx]
                        sub_result <- .run_scdblfinder_safe(
                            sub_sce,
                            label = sprintf("[BATCH:%s]", batch_nm),
                            dbr = dbr
                        )
                        sub_sce <- sub_result$sce
                        used_singlet_fallback <- used_singlet_fallback || isTRUE(sub_result$fallback)
                        if (isTRUE(sub_result$fallback) && is.null(fallback_reason))
                            fallback_reason <- sprintf("batch %s: %s", batch_nm, sub_result$reason)

                        class_results[colnames(sub_sce)] <- as.character(sub_sce$scDblFinder.class)
                        score_results[colnames(sub_sce)] <- as.numeric(sub_sce$scDblFinder.score)
                    }
                } else {
                    log_info(sprintf("  [BATCH MODE] %d internal samples — processing independently.", n_sub))
                    batch_result <- .run_scdblfinder_safe(sce, label = "[BATCH MODE]", dbr = dbr, samples = DBL_BATCH_COL)
                    sce <- batch_result$sce
                    used_singlet_fallback <- used_singlet_fallback || isTRUE(batch_result$fallback)
                    if (isTRUE(batch_result$fallback) && is.null(fallback_reason))
                        fallback_reason <- batch_result$reason

                    class_results <- as.character(sce$scDblFinder.class)
                    score_results <- as.numeric(sce$scDblFinder.score)
                    names(class_results) <- names(score_results) <- colnames(sce)
                }
            } else {
                log_info("  [SINGLE MODE] No internal batching — running standard mode.")
                single_result <- .run_scdblfinder_safe(sce, label = "[SINGLE MODE]", dbr = dbr)
                sce <- single_result$sce
                used_singlet_fallback <- used_singlet_fallback || isTRUE(single_result$fallback)
                if (isTRUE(single_result$fallback) && is.null(fallback_reason))
                    fallback_reason <- single_result$reason

                class_results <- as.character(sce$scDblFinder.class)
                score_results <- as.numeric(sce$scDblFinder.score)
                names(class_results) <- names(score_results) <- colnames(sce)
            }
        }

        if (used_singlet_fallback && !dbr_source %in% c("SKIPPED_TOO_FEW_CELLS", "SKIPPED_EMPTY_AFTER_GHOST"))
            dbr_source <- "FALLBACK_SINGLET_CALLS"

        obj$scDblFinder.class <- class_results[colnames(obj)]
        obj$scDblFinder.score <- score_results[colnames(obj)]

        # ── 3-way doublet status: singlet / ambiguous / doublet ───────────────
        # Ambiguous zone derived per-sample from KDE overlap of score distributions.
        scores  <- obj$scDblFinder.score
        classes <- obj$scDblFinder.class

        overlap_zone <- .compute_overlap_zone(
            scores  = scores,
            classes = classes,
            min_overlap_frac = DBL_AMBIGUOUS_MIN_OVERLAP,
            kde_n            = DBL_AMBIGUOUS_KDE_N
        )
        amb_lo <- overlap_zone$lo
        amb_hi <- overlap_zone$hi

        if (!is.na(amb_lo)) {
            log_info(sprintf(
                "  Ambiguous zone (data-driven): score in [%.3f, %.3f]  [%s]",
                amb_lo, amb_hi, overlap_zone$method))
            doublet_status <- ifelse(
                !is.na(scores) & scores >= amb_lo & scores <= amb_hi,
                "ambiguous",
                ifelse(
                    !is.na(classes) & classes == "doublet" &
                    !is.na(scores)  & scores  >  amb_hi,
                    "doublet",
                    "singlet"
                )
            )
        } else {
            log_info(sprintf(
                "  Ambiguous zone: none detected (singlet/doublet distributions %s) — binary mode.",
                if (overlap_zone$method == "NONE") "do not overlap above threshold"
                else "could not be estimated"))
            doublet_status <- ifelse(
                !is.na(classes) & classes == "doublet", "doublet", "singlet")
            amb_lo <- NA_real_; amb_hi <- NA_real_
        }
        obj$doublet_status <- doublet_status

        n_total     <- ncol(obj)
        n_doublets  <- sum(doublet_status == "doublet",  na.rm = TRUE)
        n_ambiguous <- sum(doublet_status == "ambiguous", na.rm = TRUE)
        n_singlets  <- sum(doublet_status == "singlet",   na.rm = TRUE)
        pct_doublet   <- round((n_doublets  / n_total) * 100, 2)
        pct_ambiguous <- round((n_ambiguous / n_total) * 100, 2)
        log_info(sprintf(
            "  Total: %d | Confident doublets: %d (%.2f%%) | Ambiguous: %d (%.2f%%) | Singlets: %d",
            n_total, n_doublets, pct_doublet, n_ambiguous, pct_ambiguous, n_singlets))
        if (used_singlet_fallback && !is.null(fallback_reason))
            log_warn(sprintf("  Doublet fallback reason: %s", fallback_reason))

        # Audit plots (4-panel): UMAP, feature violin, complexity scatter, score density
        plot_df    <- as.data.frame(obj@meta.data)
        dbl_colors <- c("singlet"   = "grey75",
                        "ambiguous" = "#f0a500",
                        "doublet"   = "firebrick3")

        plot_df$doublet_status <- factor(plot_df$doublet_status,
                                         levels = c("singlet", "ambiguous", "doublet"))

        if ("umap" %in% names(obj@reductions)) {
            umap_coords    <- as.data.frame(Embeddings(obj, "umap"))
            plot_df$UMAP_1 <- umap_coords[, 1]
            plot_df$UMAP_2 <- umap_coords[, 2]

            p_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = doublet_status)) +
                geom_point(size = 0.5, alpha = 0.8) +
                scale_color_manual(values = dbl_colors) +
                labs(title    = paste("Doublet Classification \u2014", sample_nm),
                     subtitle = sprintf(
                         "Confident doublets: %d (%.2f%%)  |  Ambiguous: %d (%.2f%%)",
                         n_doublets, pct_doublet, n_ambiguous, pct_ambiguous),
                     color = "Status") +
                theme_minimal()
        } else {
            umap_coords <- NULL
            p_umap <- ggplot() +
                annotate("text", x = 0.5, y = 0.5,
                         label = "UMAP unavailable for this sample",
                         size = 5, color = "grey50") +
                labs(title = paste("Doublet Classification \u2014", sample_nm),
                     subtitle = sprintf(
                         "Confident doublets: %d (%.2f%%)  |  Ambiguous: %d (%.2f%%)",
                         n_doublets, pct_doublet, n_ambiguous, pct_ambiguous)) +
                theme_void()
        }

        p_vln <- ggplot(plot_df, aes(doublet_status, nFeature_RNA, fill = doublet_status)) +
            geom_violin(alpha = 0.7, trim = FALSE) +
            scale_fill_manual(values = dbl_colors) +
            labs(title = "Feature Count Distribution", x = NULL, y = "nFeature_RNA") +
            theme_minimal() + theme(legend.position = "none")

        p_sct <- ggplot(plot_df, aes(nCount_RNA, nFeature_RNA, color = doublet_status)) +
            geom_point(size = 0.5, alpha = 0.5) +
            scale_color_manual(values = dbl_colors) +
            labs(title = "Complexity Scatter", color = "Status") +
            theme_minimal()

        p_dns <- ggplot(plot_df, aes(scDblFinder.score, fill = scDblFinder.class)) +
            geom_density(alpha = 0.5) +
            scale_fill_manual(values = c("singlet" = "grey75", "doublet" = "firebrick3")) +
            labs(title = "Score Distribution", fill = NULL) +
            theme_minimal()

        print((p_umap + p_sct) / (p_vln + p_dns))

        # Full object retains all cells with scDblFinder.class (binary) and
        # doublet_status (singlet/ambiguous/doublet).
        saveRDS(obj,
                file.path(DOUBLET_DIR, paste0(sample_nm, "_with_doublet_calls.rds")))

        # Singlets RDS retains singlets + ambiguous; removes only confident doublets.
        obj_singlets <- subset(obj, subset = doublet_status %in% c("singlet", "ambiguous"))
        saveRDS(obj_singlets,
                file.path(DOUBLET_DIR, paste0(sample_nm, "_singlets.rds")))

        log_info(sprintf(
            "  Saved: %s_with_doublet_calls.rds (all %d cells) | %s_singlets.rds (%d cells: %d singlets + %d ambiguous)",
            sample_nm, n_total,
            sample_nm, ncol(obj_singlets), n_singlets, n_ambiguous))
        .save_run_metadata(singlets_path, dbl_run_params, extra = list(
            n_total = n_total, n_doublets = n_doublets, pct_doublet = pct_doublet
        ))

        dbl_summary_stats[[sample_nm]] <- data.frame(
            Sample              = sample_nm,
            Total_Cells         = n_total,
            Doublets_Found      = n_doublets,
            Ambiguous_Found     = n_ambiguous,
            Singlets_Retained   = n_singlets,
            Cells_in_Singlets_RDS = n_singlets + n_ambiguous,
            Percent_Doublet     = pct_doublet,
            Percent_Ambiguous   = pct_ambiguous,
            Ambiguous_Score_Lo  = amb_lo,
            Ambiguous_Score_Hi  = amb_hi,
            DBR_Used            = ifelse(is.null(dbr), NA_real_, dbr),
            Cells_Loaded_Used   = ifelse(is.na(cells_loaded_used), NA_real_, cells_loaded_used),
            Cells_Recovered_Used = ifelse(is.na(cells_recovered_used), n_total, cells_recovered_used),
            DBR_Source          = dbr_source,
            stringsAsFactors    = FALSE
        )

        rm(obj, obj_singlets, sce, plot_df, umap_coords, p_umap, p_vln, p_sct, p_dns)
        gc()
        TRUE
        }

        }, error = function(e) {
            log_error(sprintf("[FAILED DOUBLET] %s: %s", sample_nm, conditionMessage(e)))

            dbl_summary_stats[[sample_nm]] <<- data.frame(
                Sample              = sample_nm,
                Total_Cells         = NA,
                Doublets_Found      = NA,
                Ambiguous_Found     = NA,
                Singlets_Retained   = NA,
                Cells_in_Singlets_RDS = NA,
                Percent_Doublet     = NA,
                Percent_Ambiguous   = NA,
                Ambiguous_Score_Lo  = NA_real_,
                Ambiguous_Score_Hi  = NA_real_,
                DBR_Used            = NA_real_,
                Cells_Loaded_Used   = NA_real_,
                Cells_Recovered_Used = NA_real_,
                DBR_Source          = "ERROR",
                stringsAsFactors    = FALSE
            )
            FALSE
        })
    }

    tryCatch(dev.off(), error = function(e) log_warn(paste("PDF device close:", conditionMessage(e))))
    log_info(sprintf("Doublet audit PDF saved: %s", dbl_pdf))

    if (length(dbl_summary_stats) > 0) {
        dbl_summary_df  <- do.call(rbind, dbl_summary_stats)
        dbl_summary_csv <- file.path(DOUBLET_SUMMARY_DIR, paste0("doublet_summary_", run_ts, ".csv"))
        write.csv(dbl_summary_df, dbl_summary_csv, row.names = FALSE)
        log_info(sprintf("Doublet summary CSV saved: %s", dbl_summary_csv))
        print(dbl_summary_df)
    }

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

}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6.5 — H5AD Conversion
# ══════════════════════════════════════════════════════════════════════════════
.run_stage_h5ad <- function(input_dir, output_dir) {
    log_info("──────────────────────────────────────────────────")
    log_info("STAGE 3: H5AD CONVERSION")
    log_info(sprintf("  Input  : %s", input_dir))
    log_info(sprintf("  Output : %s", output_dir))
    log_info("──────────────────────────────────────────────────")

    rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(rds_files) == 0) {
        log_warn(sprintf("No RDS files found in %s for H5AD conversion.", input_dir))
        return()
    }

    for (rds_path in rds_files) {
        base_name <- sub("\\.rds$", "", basename(rds_path), ignore.case = TRUE)
        h5ad_path <- file.path(output_dir, paste0(base_name, ".h5ad"))

        if (!.as_bool(cfg$force_overwrite, FALSE) && file.exists(h5ad_path) && file.info(h5ad_path)$size > 0) {
            log_info(sprintf("[H5AD SKIPPED] %s.h5ad already exists.", base_name))
            next
        }

        log_info(sprintf("[H5AD] Converting: %s", base_name))
        success <- tryCatch({
            seurat_obj <- readRDS(rds_path)

            # Handle Seurat v5 Assays (convert to strictly safe Assay object to flatten layers)
            if (inherits(seurat_obj[["RNA"]], "Assay5")) {
                seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
            }

            # Extract count matrix and transpose (AnnData expects Cells x Genes)
            counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
            counts_t <- t(counts)

            # Extract metadata
            meta_data <- seurat_obj@meta.data

            # Create AnnData object
            ad <- anndata::AnnData(
                X = counts_t,
                obs = meta_data
            )

            # Write to h5ad
            ad$write_h5ad(h5ad_path)
            log_info(sprintf("  Saved: %s", h5ad_path))
            rm(seurat_obj, counts, counts_t, meta_data, ad); gc()
            TRUE
        }, error = function(e) {
            log_error(sprintf("[FAILED H5AD] %s: %s", base_name, conditionMessage(e)))
            FALSE
        })
    }
    log_info("H5AD Conversion Complete.")
}

# ── Stage execution dispatcher ─────────────────────────────────────────────────
reverse_mode_run <- .as_bool(cfg$reverse_mode %||% "TRUE", TRUE)

if (QC_RUN == "all") {
    if (reverse_mode_run) {
        .run_stage_doublet(doublet_input_dir = rds_dir)
        .run_stage_qc(qc_input_dir = DOUBLET_DIR)
    } else {
        .run_stage_qc(qc_input_dir = rds_dir)
        .run_stage_doublet(doublet_input_dir = FILTERED_DIR)
    }
} else if (QC_RUN == "qc") {
    .run_stage_qc(qc_input_dir = if (reverse_mode_run) DOUBLET_DIR else rds_dir)
} else if (QC_RUN == "doublet") {
    .run_stage_doublet(doublet_input_dir = if (reverse_mode_run) rds_dir else FILTERED_DIR)
}

# ── Optional Stage: H5AD conversion ───────────────────────────────────────────
if (RUN_H5AD) {
    h5ad_input_dir <- if (QC_RUN == "qc") {
        FILTERED_DIR
    } else if (QC_RUN == "doublet") {
        DOUBLET_DIR
    } else if (reverse_mode_run) {
        FILTERED_DIR
    } else {
        DOUBLET_DIR
    }
    .run_stage_h5ad(h5ad_input_dir, H5AD_DIR)
}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — Pipeline Complete
# ══════════════════════════════════════════════════════════════════════════════

log_info("══════════════════════════════════════════════════")
log_info("PIPELINE COMPLETE")
log_info(sprintf("  QC run     : %s", QC_RUN))
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
# Writes a per-sample CSV (cells_raw/post_doublet/post_qc + % removed per
# stage), per-plot PNGs, and a multi-page summary PDF.
# ══════════════════════════════════════════════════════════════════════════════

log_info("──────────────────────────────────────────────────")
log_info("SECTION 8: INTEGRATED DOUBLET + QC SUMMARY")
log_info("──────────────────────────────────────────────────")

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

    if (has_qc_summary) {
        qc_sum <- .safe_read_csv(qc_csv_path, "QC summary")
        if (!is.null(qc_sum) && "sample" %in% colnames(qc_sum)) {
            qc_sum <- qc_sum[qc_sum$sample != "TOTAL", ]
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

    if (has_dbl_summary) {
        dbl_sum <- .safe_read_csv(dbl_csv_path, "Doublet summary")
        if (!is.null(dbl_sum) && "DBR_Source" %in% colnames(dbl_sum)) {
            dbl_sum <- dbl_sum[dbl_sum$DBR_Source != "ERROR", ]
        }
    } else {
        dbl_sum <- NULL
    }
    if (is.null(dbl_sum)) {
        dbl_sum <- data.frame(Sample = character(0), Total_Cells = integer(0), Doublets_Found = integer(0))
        log_warn("Section 8: Doublet summary CSV missing/unreadable — doublet columns will be NA.")
    }

    reverse_mode <- .as_bool(cfg$reverse_mode %||% "TRUE", TRUE)

    if (reverse_mode) {
        # Sequence: raw → ghost filter → DBR → doublet removal → QC

        # Raw counts: try cache CSV first, fall back to reading RDS files.
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

        dbl_tbl <- data.frame(sample = character(0),
                              cells_after_ghost = numeric(0),
                              dbr_used = numeric(0),
                              doublets_found = numeric(0),
                              ambiguous_found = numeric(0),
                              pct_doublet = numeric(0),
                              pct_ambiguous = numeric(0),
                              stringsAsFactors = FALSE)
        if (nrow(dbl_sum) > 0) {
            dbl_tbl <- data.frame(
                sample = .canon_sample_name(dbl_sum$Sample),
                cells_after_ghost = suppressWarnings(as.numeric(dbl_sum$Total_Cells)),
                dbr_used = if ("DBR_Used" %in% colnames(dbl_sum)) suppressWarnings(as.numeric(dbl_sum$DBR_Used)) else NA_real_,
                doublets_found = suppressWarnings(as.numeric(dbl_sum$Doublets_Found)),
                ambiguous_found = if ("Ambiguous_Found" %in% colnames(dbl_sum)) suppressWarnings(as.numeric(dbl_sum$Ambiguous_Found)) else NA_real_,
                pct_doublet = suppressWarnings(as.numeric(dbl_sum$Percent_Doublet)),
                pct_ambiguous = if ("Percent_Ambiguous" %in% colnames(dbl_sum)) suppressWarnings(as.numeric(dbl_sum$Percent_Ambiguous)) else NA_real_,
                stringsAsFactors = FALSE
            )
        }

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
        # Normal order: raw → QC → doublet
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
            dbl_join_cols <- intersect(c("Sample", "Total_Cells", "Doublets_Found",
                                         "Ambiguous_Found", "Percent_Doublet", "Percent_Ambiguous"),
                                       colnames(dbl_sum))
            dbl_join <- dbl_sum[, dbl_join_cols, drop = FALSE]
            if (!"Ambiguous_Found"   %in% colnames(dbl_join)) dbl_join$Ambiguous_Found   <- NA_real_
            if (!"Percent_Ambiguous" %in% colnames(dbl_join)) dbl_join$Percent_Ambiguous <- NA_real_
            colnames(dbl_join)[colnames(dbl_join) == "Sample"]            <- "sample"
            colnames(dbl_join)[colnames(dbl_join) == "Total_Cells"]       <- "cells_post_qc_dbl"
            colnames(dbl_join)[colnames(dbl_join) == "Doublets_Found"]    <- "doublets_found"
            colnames(dbl_join)[colnames(dbl_join) == "Ambiguous_Found"]   <- "ambiguous_found"
            colnames(dbl_join)[colnames(dbl_join) == "Percent_Doublet"]   <- "pct_doublet"
            colnames(dbl_join)[colnames(dbl_join) == "Percent_Ambiguous"] <- "pct_ambiguous"
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

    if (!("cells_after_ghost" %in% colnames(int_df))) int_df$cells_after_ghost <- NA_real_
    if (!("dbr_used" %in% colnames(int_df)))          int_df$dbr_used <- NA_real_
    if (!("cells_after_doublet" %in% colnames(int_df))) int_df$cells_after_doublet <- int_df$cells_post_doublet
    if (!("cells_after_qc" %in% colnames(int_df)))      int_df$cells_after_qc <- int_df$cells_post_qc

    int_df$sample_short <- int_df$sample
    int_df$sample_short <- gsub("d10_1016_j_",             "", int_df$sample_short)
    int_df$sample_short <- gsub("d10_1126_sciadv_",        "", int_df$sample_short)
    int_df$sample_short <- gsub("d10_1038_",               "", int_df$sample_short)
    int_df$sample_short <- gsub("dno_doi_kidney_organoid_", "kidney_", int_df$sample_short)

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

    # Schema-safe rbind: add any missing columns to either frame before binding.
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

    dir.create(INTEGRATED_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

    int_csv_path <- file.path(INTEGRATED_SUMMARY_DIR, "integrated_qc_doublet_summary.csv")

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
            ambiguous_found = get_num("ambiguous_found", "max"),
            pct_doublets = get_num("pct_doublet", "first"),
            pct_ambiguous = get_num("pct_ambiguous", "first"),
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
            ambiguous_found = sum(int_csv_df$ambiguous_found, na.rm = TRUE),
            pct_doublets = NA_real_,
            pct_ambiguous = NA_real_,
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

    log_info("=== Integrated QC + Doublet Summary ===")
    print_int_cols <- c(
        "sample", "cells_raw", "cells_after_ghost", "dbr_used", "doublets_found",
        "ambiguous_found", "pct_doublets", "pct_ambiguous",
        "cells_after_doublet_removal", "cells_after_qc",
        "pct_qc_removal", "total_cells_removal", "pct_total_cells_removal"
    )
    for (ln in capture.output(print(int_csv_df[, print_int_cols], row.names = FALSE))) log_info(ln)

    # ── Dashboard plots ───────────────────────────────────────────────────────
    pd <- int_df[complete.cases(int_df[, c("cells_raw", "cells_post_qc")]), ]

    if (nrow(pd) == 0) {
        log_warn("Section 8: No complete sample rows — skipping dashboard plots.")
    } else {

        # Plot 1: Stacked % loss per stage
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

        # Plot 2: Cumulative funnel (connected lines)
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

        # Plot 3: Doublet % vs QC removal % scatter
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

        # Plot 4: Donut — cohort-level cell budget
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
            xlim(0.5, 2.5) +
            scale_fill_manual(values = c(
                "Retained"           = "#1a9641",
                "Removed by QC"      = "#fc8d59",
                "Removed by Doublet" = "#d73027")) +
            geom_text(aes(label = label),
                      position = position_stack(vjust = 0.5),
                      size = 3, color = "white", fontface = "bold") +
            labs(title    = "Cohort-Level Cell Budget",
                 subtitle = sprintf("Total input is: %s cells", scales::comma(sum(donut_df$count))),
                 fill = "") +
            theme_void(base_size = 10) +
            theme(legend.position = "right",
                  plot.title      = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle   = element_text(hjust = 0.5))

        int_png_dir <- file.path(INTEGRATED_SUMMARY_DIR, "plots")
        dir.create(int_png_dir, recursive = TRUE, showWarnings = FALSE)

        ggsave(file.path(int_png_dir, "01_stacked_loss.png"),           plot = p8_2, width = 12, height = 7, dpi = PLOT_DPI, bg = PLOT_BG)
        ggsave(file.path(int_png_dir, "02_funnel.png"),                 plot = p8_3, width = 12, height = 7, dpi = PLOT_DPI, bg = PLOT_BG)
        ggsave(file.path(int_png_dir, "03_scatter_doublet_vs_qc.png"),  plot = p8_8, width = 9,  height = 7, dpi = PLOT_DPI, bg = PLOT_BG)
        ggsave(file.path(int_png_dir, "04_donut_budget.png"),           plot = p8_9, width = 9,  height = 7, dpi = PLOT_DPI, bg = PLOT_BG)

        log_info(sprintf("Integrated dashboard plots saved as individual PNGs in: %s", int_png_dir))

        int_pdf_path <- file.path(INTEGRATED_SUMMARY_DIR, "integrated_summary_plots.pdf")
        pdf(int_pdf_path, width = 14, height = 9, onefile = TRUE)

        print(p8_2)
        print(p8_3)
        print(p8_8)
        print(p8_9)

        tryCatch(dev.off(), error = function(e) log_warn(paste("Integrated PDF close:", conditionMessage(e))))
        log_info(sprintf("Integrated dashboard PDF saved: %s", int_pdf_path))

    }

    log_info("──────────────────────────────────────────────────")
    log_info("SECTION 8: INTEGRATED SUMMARY COMPLETE")
    log_info(sprintf("  CSV       : %s", int_csv_path))
    if (exists("int_png_dir"))  log_info(sprintf("  PNG DIR   : %s", int_png_dir))
    if (exists("int_pdf_path")) log_info(sprintf("  PDF       : %s", int_pdf_path))
    log_info("──────────────────────────────────────────────────")

}

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Final Pipeline Footer
# ══════════════════════════════════════════════════════════════════════════════

log_info("══════════════════════════════════════════════════")
log_info("ALL SECTIONS COMPLETE")
log_info(sprintf("  QC run            : %s", QC_RUN))
log_info(sprintf("  Output base       : %s", BASE_OUT_DIR))
if (RUN_QC)      log_info(sprintf("  QC audit PDF      : %s", file.path(QC_DIR,              "QC_Before_After_Report.pdf")))
if (RUN_QC)      log_info(sprintf("  QC dashboard      : %s", file.path(QC_DIR,              "summary_qc_full_dashboard.png")))
if (RUN_QC)      log_info(sprintf("  QC summary CSV    : %s", file.path(QC_DIR,              "qc_summary_detailed.csv")))
if (RUN_DOUBLET) log_info(sprintf("  Doublet PDF       : %s", file.path(DOUBLET_SUMMARY_DIR, paste0("Doublet_Audit_Report_", run_ts, ".pdf"))))
if (RUN_DOUBLET) log_info(sprintf("  Doublet CSV       : %s", file.path(DOUBLET_SUMMARY_DIR, paste0("doublet_summary_",      run_ts, ".csv"))))
if (RUN_H5AD)    log_info(sprintf("  H5AD outputs      : %s", H5AD_DIR))
log_info(sprintf("  Integrated CSV    : %s", file.path(INTEGRATED_SUMMARY_DIR, "integrated_qc_doublet_summary.csv")))
log_info(sprintf("  Integrated PNG dir: %s", file.path(INTEGRATED_SUMMARY_DIR, "plots")))
log_info(sprintf("  Integrated PDF    : %s", file.path(INTEGRATED_SUMMARY_DIR, "integrated_summary_plots.pdf")))
log_info(sprintf("  Log file          : %s", log_file_path))
log_info("══════════════════════════════════════════════════")

}
log_info("══════════════════════════════════════════════════")