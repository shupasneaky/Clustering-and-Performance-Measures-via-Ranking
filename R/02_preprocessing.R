#===================================================================================================
# 02_preprocess_data.R is meant to preprocess each single cell dataset:
# - Loads raw RData files containing matrices named "NAME_data"
# - Converts matrices to numeric and replaces NAs with 0
# - Filters out genes (rows) and cells (columns) with no counts
# - Adds a pseudo count and log-normalizes the data
# - Saves the processed SingleCellExperiment as "NAME_data_preprocessed.RData"
# - Cleans up the environment to free memory
#===================================================================================================

# Packages
library(SingleCellExperiment)
library(scuttle)

# Preprocessing function: filter zero counts, add pseudo count, and log-normalize
preprocess <- function(d) {
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(d)))
  sce <- scuttle::addPerFeatureQCMetrics(sce)
  sce <- scuttle::addPerCellQCMetrics(sce)
  sce <- sce[rowSums(counts(sce)) > 0, ]
  sce <- sce[, colSums(counts(sce)) > 0]
  counts(sce) <- counts(sce) + 0.0001
  sce <- scuttle::logNormCounts(sce)
  return(sce)
}

# List of files to process
files <- c(
  "Beutner_mice_age_sc_data/BEUTNER_data.RData",
  "Yan_sc_data/YAN_data.RData",
  "Pollen_sc_data/POLLEN_data.RData",
  "Goolam_sc_data/GOOLAM_data.RData",
  "Deng_sc_data/DENG_data.RData",
  "Baise_sc_data/BAISE_data.RData"
)

for(file_path in files) {
  load(file_path)  # Loads an object named like "NAME_data"
  
  # Extract base name (e.g., "BEUTNER_data") from file path
  base_name <- tools::file_path_sans_ext(basename(file_path))
  message("Processing ", base_name)
  
  # Retrieve, convert, and clean the matrix
  mat_obj <- get(base_name)
  mat_obj <- as.matrix(mat_obj)
  storage.mode(mat_obj) <- "numeric"
  mat_obj[is.na(mat_obj)] <- 0
  
  # Preprocess and save the SingleCellExperiment object
  sce <- preprocess(mat_obj)
  new_file_name <- file.path(dirname(file_path), paste0(base_name, "_preprocessed.RData"))
  save(sce, file = new_file_name)
  
  # Remove objects to free up memory
  rm(list = c(base_name, "mat_obj", "sce"))
  gc()
}
