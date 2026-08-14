# Export PCA / Harmony / RPCA / CCA embeddings + cell metadata for scIB.

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
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

.log_info <- function(...) cat(sprintf("[%s] INFO  %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))
.log_warn <- function(...) cat(sprintf("[%s] WARN  %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

.reduction_for <- function(method) {
    switch(method,
        unintegrated = "pca",
        harmony = "harmony",
        rpca = "integrated.rpca",
        cca = "integrated.cca",
        NA_character_
    )
}

.write_embeddings <- function(obj, reduction, out_csv, npcs) {
    if (!reduction %in% Reductions(obj)) {
        stop(sprintf("Reduction '%s' not found in object.", reduction))
    }
    emb <- as.matrix(Embeddings(obj, reduction))
    if (ncol(emb) > npcs) emb <- emb[, seq_len(npcs), drop = FALSE]
    write.csv(emb, out_csv, row.names = TRUE)
    .log_info("Wrote %s (%d cells x %d dims).", out_csv, nrow(emb), ncol(emb))
}

.write_metadata <- function(obj, out_csv) {
    keep <- intersect(
        c("sample_id", "orig.ident", "Level_1", "Level_2", "Level_3",
          "Level_3_latest", "celltype_final"),
        colnames(obj@meta.data)
    )
    meta <- obj@meta.data[, keep, drop = FALSE]
    write.csv(meta, out_csv, row.names = TRUE)
    .log_info("Wrote cell metadata: %s (columns: %s).", out_csv, paste(keep, collapse = ", "))
}

.metadata_has_labels <- function(meta_csv) {
    if (!file.exists(meta_csv) || file.info(meta_csv)$size <= 0) return(FALSE)
    meta <- tryCatch(
        read.csv(meta_csv, row.names = 1, check.names = FALSE, nrows = 1),
        error = function(e) NULL
    )
    if (is.null(meta)) return(FALSE)
    label_cols <- c("Level_3", "Level_3_latest", "celltype_final", "Level_2", "Level_1")
    any(label_cols %in% colnames(meta))
}

tryCatch({
    suppressPackageStartupMessages(library(Seurat))
    cli <- .parse_cli_args(commandArgs(trailingOnly = TRUE))
    output_dir <- cli$output_dir %||% stop("Provide --output_dir")
    prepared_rds <- cli$prepared_rds %||% NA_character_
    methods_raw <- tolower(trimws(cli$methods %||% "harmony,rpca,cca"))
    methods <- trimws(unlist(strsplit(methods_raw, ",", fixed = TRUE)))
    methods <- methods[nzchar(methods)]
    force_flag <- .as_bool(cli$force_overwrite, FALSE)
    npcs <- as.integer(cli$npcs %||% 30L)

    bench_dir <- file.path(output_dir, "benchmark")
    obj_dir <- file.path(output_dir, "objects")
    dir.create(bench_dir, recursive = TRUE, showWarnings = FALSE)

    meta_csv <- file.path(bench_dir, "cell_metadata.csv")
    unint_csv <- file.path(bench_dir, "unintegrated_embeddings.csv")

    if (!is.na(prepared_rds) && nzchar(prepared_rds) && file.exists(prepared_rds)) {
        need_meta <- force_flag || !file.exists(meta_csv) || file.info(meta_csv)$size <= 0 ||
            !.metadata_has_labels(meta_csv)
        need_pca <- force_flag || !file.exists(unint_csv) || file.info(unint_csv)$size <= 0
        if (need_meta || need_pca) {
            .log_info("Loading prepared RDS for unintegrated PCA / metadata.")
            obj_prep <- readRDS(prepared_rds)
            if (need_meta) .write_metadata(obj_prep, meta_csv)
            if (need_pca) .write_embeddings(obj_prep, "pca", unint_csv, npcs)
            rm(obj_prep); gc(verbose = FALSE)
        } else {
            .log_info("Skipping unintegrated export; files exist.")
        }
    } else {
        .log_warn("prepared_rds missing; unintegrated PCA will not be exported.")
    }

    for (method in methods) {
        out_csv <- file.path(bench_dir, sprintf("%s_embeddings.csv", method))
        if (!force_flag && file.exists(out_csv) && file.info(out_csv)$size > 0) {
            .log_info("[%s] Skipping embedding export; exists: %s", method, out_csv)
            next
        }
        rds <- file.path(obj_dir, sprintf("integrated_%s.rds", method))
        if (!file.exists(rds) || file.info(rds)$size <= 0) {
            .log_warn("[%s] Integrated RDS missing: %s", method, rds)
            next
        }
        red <- .reduction_for(method)
        .log_info("[%s] Loading %s", method, rds)
        obj <- readRDS(rds)
        .write_embeddings(obj, red, out_csv, npcs)
        if (!file.exists(meta_csv) || file.info(meta_csv)$size <= 0) {
            .write_metadata(obj, meta_csv)
        }
        rm(obj); gc(verbose = FALSE)
    }

    .log_info("Benchmark embedding export complete: %s", bench_dir)
    quit(save = "no", status = 0)
}, error = function(e) {
    cat(sprintf("[%s] ERROR %s\n", format(Sys.time(), "%H:%M:%S"), conditionMessage(e)), file = stderr())
    quit(save = "no", status = 1)
})
