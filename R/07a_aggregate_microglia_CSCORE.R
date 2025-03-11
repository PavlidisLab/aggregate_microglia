## Aggregate human and mouse microglial CSCORE output. This script requires that
## all microglial dataset processing scripts in R/CSCORE/ have been ran.
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

dat_dir <- file.path(data_out_dir, "CSCORE")
mcg_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)
pc_hg <- read.delim(ref_hg_path)
pc_mm <- read.delim(ref_mm_path)

# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# Dataset IDs and ensure CSCORE has been ran for each 
ids_hg <- filter(mcg_meta, Species == "Human")$ID
ids_mm <- filter(mcg_meta, Species == "Mouse")$ID
dat_paths <- list.files(dat_dir)
stopifnot(all(paste0(c(ids_hg, ids_mm), "_CSCORE.RDS") %in% dat_paths))

# Genes to keep after filtering
keep_hg <- count_summ$Human$Filter_genes
keep_mm <- count_summ$Mouse$Filter_genes 


# Here aggregating using all rank approach
aggregate_cscore <- function(ids, dat_dir, gene_df) {
  
  # Tracking matrix
  amat <- init_agg_mat(gene_df)
  
  for (id in ids) {
    
    message(paste(id, Sys.time()))
    
    # Read in CSCORE result list and get estimated correlation matrix
    path <- file.path(dat_dir, paste0(id, "_CSCORE.RDS"))
    res <- readRDS(path)
    mat <- res$est
    
    # Matrix of input genes, inserting current values where possible
    mat_keep <- matrix(0, nrow  = nrow(gene_df), ncol = nrow(gene_df))
    rownames(mat_keep) <- colnames(mat_keep) <- gene_df$Symbol
    common <- intersect(rownames(mat), gene_df$Symbol)
    mat_keep[common, common] <- mat[common, common]
    
    # Prepare mat for ranking
    mat_keep <- mat_keep %>% 
      diag_to_one() %>%
      uppertri_to_na()
    
    # Rank + standardize  
    rmat <- allrank_mat(-mat_keep, ties_arg = "min")
    rmat <- rmat / max(rmat, na.rm = TRUE)
    
    # Add current matrix to tracking matrix
    amat <- amat + rmat
    rm(mat_keep, mat, res, rmat)
    gc(verbose = FALSE)
  }
  
  # Finalize agg mat
  amat <- allrank_mat(-amat, ties_arg = "min")
  amat <- amat / max(amat, na.rm = TRUE)
  amat <- diag_to_one(amat)  # Make aggregate symmetric and self cor = 1
  amat <- lowertri_to_symm(amat)
  
  return(amat)
}



save_function_results(
  path = cscore_hg_path,
  fun = aggregate_cscore,
  args = list(
    ids = ids_hg,
    dat_dir = dat_dir,
    gene_df = keep_hg
  )
)



save_function_results(
  path = cscore_mm_path,
  fun = aggregate_cscore,
  args = list(
    ids = ids_mm,
    dat_dir = dat_dir,
    gene_df = keep_mm
  )
)
