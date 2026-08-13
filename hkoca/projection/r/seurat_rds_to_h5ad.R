# Export RNA counts from a Seurat RDS (e.g. integration prep sct_prepared.rds)
# to a sparse MTX + obs/var tables that Python can load as AnnData.

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

tryCatch({
    args <- .parse_cli_args(commandArgs(trailingOnly = TRUE))
    input_rds <- args$input_rds %||% args$input
    out_dir <- args$output_dir %||% args$output
    if (is.null(input_rds) || is.null(out_dir))
        stop("Usage: seurat_rds_to_h5ad.R --input_rds PATH --output_dir DIR")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    suppressPackageStartupMessages({
        library(Seurat)
        library(Matrix)
    })
    obj <- readRDS(input_rds)
    if (!inherits(obj, "Seurat")) stop("input_rds must be a Seurat object.")

    assay <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
    DefaultAssay(obj) <- assay
    assay_obj <- obj[[assay]]
    if (inherits(assay_obj, "Assay5") && exists("JoinLayers", mode = "function")) {
        obj[[assay]] <- JoinLayers(assay_obj)
    }
    cm <- tryCatch(
        GetAssayData(obj, assay = assay, layer = "counts"),
        error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
    )
    if (is.null(cm) || !nrow(cm) || !ncol(cm))
        stop("Could not extract RNA counts from the Seurat object.")

    Matrix::writeMM(cm, file.path(out_dir, "counts.mtx"))
    write.table(
        data.frame(gene = rownames(cm), stringsAsFactors = FALSE),
        file.path(out_dir, "var.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )
    meta <- obj@meta.data
    meta$barcode <- colnames(obj)
    write.csv(meta, file.path(out_dir, "obs.csv"), row.names = FALSE)
    cat(sprintf("Wrote counts from assay '%s': %d genes x %d cells\n", assay, nrow(cm), ncol(cm)))
}, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)), file = stderr())
    quit(save = "no", status = 1)
})
