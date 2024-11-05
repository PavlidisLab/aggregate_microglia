## TODO: save .RDS or just .tsv?
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_df <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
meta <- read.delim(mcg_meta_path) %>% filter(Species == "Mouse")

# Paths for the Fisher's Z, Rank Sum Rank, and NA tracking matrices
fz_path <- file.path(cmat_dir_mm, "aggregate_cormat_FZ_mm.RDS")
rsr_path <- file.path(cmat_dir_mm, "aggregate_cormat_RSR_mm.RDS")

# TODO: keep?
fz_tsv_path <- file.path(cmat_dir_mm, "aggregate_cormat_FZ_mm.tsv")
rsr_tsv_path <- file.path(cmat_dir_mm, "aggregate_cormat_RSR_mm.tsv")
na_tsv_path <- file.path(cmat_dir_mm, "aggregate_cormat_NA_mm.tsv")




# Generate and save
# ------------------------------------------------------------------------------



save_function_results(
  path = fz_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df,
    cor_method = "pearson",
    agg_method = "FZ"
  )
)



save_function_results(
  path = rsr_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)



# FZ matrix as .tsv files
if (!file.exists(fz_tsv_path)) {
  
  agg_l <- readRDS(fz_path)
  fwrite_mat(agg_l$Agg_mat, fz_tsv_path)
  
  if (!file.exists(na_tsv_path)) {
    fwrite_mat(agg_l$NA_mat, na_tsv_path)
  }

}



# RSR matrix as .tsv files
if (!file.exists(rsr_tsv_path)) {
  
  agg_l <- readRDS(fz_path)
  fwrite_mat(agg_l$Agg_mat, rsr_tsv_path)
  
  if (!file.exists(na_tsv_path)) {
    fwrite_mat(agg_l$NA_mat, na_tsv_path)
  }
  
}
