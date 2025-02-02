## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
# library(aggtools)
source("/home/amorin/Projects/aggtools/R/aggregate_functions.R")
source("R/00_config.R")
source("R/utils/functions.R")

dat_dir <- file.path(data_out_dir, "CSCORE")
mcg_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)
pc_hg <- read.delim(ref_hg_path)
pc_mm <- read.delim(ref_mm_path)

# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# Completed IDs
ids_hg <- filter(mcg_meta, Species == "Human")$ID
ids_mm <- filter(mcg_meta, Species == "Mouse")$ID
dat_paths <- list.files(dat_dir)
comp_ids_hg <- intersect(ids_hg, str_replace_all(dat_paths, "_CSCORE.RDS", ""))
comp_ids_mm <- intersect(ids_mm, str_replace_all(dat_paths, "_CSCORE.RDS", ""))

# Outpaths
cscore_hg_path <- file.path(data_out_dir, "CSCORE_aggregate_microglia_hg.RDS")
cscore_mm_path <- file.path(data_out_dir, "CSCORE_aggregate_microglia_mm.RDS")


# Min genes to keep
keep_hg <- count_summ$Human$Summ_df %>% 
  filter(N_msr >= floor(max(N_msr) * 1/3)) %>%
  select(Symbol)

keep_mm <- count_summ$Mouse$Summ_df %>% 
  filter(N_msr >= floor(max(N_msr) * 1/3)) %>%
  select(Symbol)


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
    ids = comp_ids_hg,
    dat_dir = dat_dir,
    gene_df = keep_hg
  )
)



save_function_results(
  path = cscore_mm_path,
  fun = aggregate_cscore,
  args = list(
    ids = comp_ids_mm,
    dat_dir = dat_dir,
    gene_df = keep_mm
  )
)
