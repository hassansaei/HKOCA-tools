# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

# Define paths
INPUT_DIR  <- "/data-master/pure-workspace/labss/hmami/new_data/HKOCA-tools/scripts/qc/test/qc_results_after_harmonization_final/qc_filtered_rds"
OUTPUT_DIR <- "/data-master/pure-workspace/labss/hmami/new_data/HKOCA-tools/scripts/qc/test/qc_results_after_harmonization_final/qc_filtered_h5ad_converted"

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  message("Created output directory: ", OUTPUT_DIR)
}

# Get list of RDS files
rds_files <- list.files(INPUT_DIR, pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) {
  stop("No .rds files found in the specified INPUT_DIR.")
}

message(sprintf("Found %d RDS files to process.", length(rds_files)))
message("--------------------------------------------------")

for (rds_path in rds_files) {
  
  file_name <- basename(rds_path)
  base_name <- tools::file_path_sans_ext(file_name)
  
  h5seurat_path <- file.path(OUTPUT_DIR, paste0(base_name, ".h5Seurat"))
  h5ad_path     <- file.path(OUTPUT_DIR, paste0(base_name, ".h5ad"))
  
  message("Processing: ", file_name)
  message("File size: ", round(file.info(rds_path)$size / (1024^3), 2), " GB")
  
  # ==========================================
  # 1. PRE-CONVERSION SANITY CHECKS
  # ==========================================
  message("  -> Loading RDS object...")
  obj <- readRDS(rds_path)
  
  # Ensure it is a Seurat object
  if (!inherits(obj, "Seurat")) {
    warning("  -> Object is not a Seurat object. Skipping ", file_name)
    rm(obj); gc()
    next
  }
  
  # If using Seurat v5, SeuratDisk requires standard v3 assays to convert properly
  # This cast operation is completely lossless for the data you want to export
  if (inherits(obj[["RNA"]], "Assay5")) {
    message("  -> Seurat v5 Assay detected. Casting to v3 Assay for compatibility with SeuratDisk...")
    obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay")
  }
  
  # Record metrics for sanity check
  pre_cells <- ncol(obj)
  pre_features <- nrow(obj)
  
  # CAPTURE COLUMN NAMES BEFORE DELETION
  pre_meta_colnames <- colnames(obj@meta.data) 
  pre_meta_cols <- length(pre_meta_colnames)
  
  message(sprintf("  -> Pre-conversion metrics: %d features, %d cells, %d metadata columns", 
                  pre_features, pre_cells, pre_meta_cols))
  
  # ==========================================
  # 2. CONVERSION
  # ==========================================
  # Step A: Convert RDS to h5Seurat (intermediate format)
  message("  -> Saving to intermediate h5Seurat format...")
  SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE)
  
  # Clear the heavy Seurat object from RAM immediately before the final h5ad conversion
  rm(obj)
  gc(reset = TRUE) # Double garbage collection to force RAM release to OS
  gc()
  
  # Step B: Convert h5Seurat to h5ad
  message("  -> Converting h5Seurat to h5ad...")
  Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
  
  # Clean up the intermediate file to save disk space
  if (file.exists(h5seurat_path)) {
    file.remove(h5seurat_path)
  }
  
  # ==========================================
  # 3. POST-CONVERSION SANITY CHECKS
  # ==========================================
  if (file.exists(h5ad_path)) {
    h5ad_size <- round(file.info(h5ad_path)$size / (1024^3), 2)
    
    if (h5ad_size > 0) {
      message("  -> [SUCCESS] h5ad file created successfully.")
      message("  -> Final file size: ", h5ad_size, " GB")
      
      # ---------------------------------------------------------
      # Robust Metadata Comparison using hdf5r
      # ---------------------------------------------------------
      message("  -> Verifying metadata integrity in h5ad...")
      
      tryCatch({
        # Open the h5ad file in read-only mode (extremely memory efficient)
        h5_file <- hdf5r::H5File$new(h5ad_path, mode = "r")
        
        # Cell metadata in h5ad is stored in the 'obs' group
        if ("obs" %in% names(h5_file)) {
          obs_obj <- h5_file[["obs"]]
          
          # Extract column names (h5ad stores dataframe columns as datasets within the 'obs' group)
          h5ad_meta_colnames <- names(obs_obj)
          
          # Clean up AnnData-specific internal keys that shouldn't be compared
          h5ad_meta_colnames <- h5ad_meta_colnames[!h5ad_meta_colnames %in% c("_index", "__categories")]
          
          # Compare Original RDS columns vs H5AD columns
          missing_cols <- setdiff(pre_meta_colnames, h5ad_meta_colnames)
          extra_cols   <- setdiff(h5ad_meta_colnames, pre_meta_colnames)
          
          if (length(missing_cols) == 0) {
            message("  -> [SUCCESS] All ", length(pre_meta_colnames), " metadata columns successfully transferred.")
          } else {
            warning("  -> [WARNING] Missing metadata columns in h5ad: ", paste(missing_cols, collapse = ", "))
          }
          
          # Optional: Flag if SeuratDisk generated any unexpected extra columns
          if (length(extra_cols) > 0) {
            message("  -> [NOTE] Extra columns found in h5ad (often index or converted factors): ", 
                    paste(extra_cols, collapse = ", "))
          }
          
        } else {
          warning("  -> [WARNING] 'obs' (metadata) group not found in the h5ad file!")
        }
        
        # Always close the HDF5 connection
        h5_file$close_all()
        
      }, error = function(e) {
        warning("  -> [WARNING] Could not read h5ad file to verify metadata. Error: ", e$message)
        # Ensure connection is closed even if it fails
        if (exists("h5_file") && inherits(h5_file, "H5File")) {
          h5_file$close_all()
        }
      })
      
    } else {
      warning("  -> [WARNING] h5ad file was created but size is 0 bytes!")
    }
  } else {
    warning("  -> [ERROR] h5ad file was not created!")
  }
  
  message("--------------------------------------------------")
}

message("All files processed successfully.")