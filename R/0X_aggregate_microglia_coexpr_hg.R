## TODO: save .RDS or just .tsv?
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_df <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
meta <- read.delim(mcg_meta_path) %>% filter(Species == "Human")

# Paths for the Fisher's Z, Rank Sum Rank, and NA tracking matrices
fz_path <- file.path(cmat_dir_hg, "aggregate_cormat_FZ_hg.RDS")
rsr_path <- file.path(cmat_dir_hg, "aggregate_cormat_RSR_hg.RDS")

# TODO: keep?
fz_tsv_path <- file.path(cmat_dir_hg, "aggregate_cormat_FZ_hg.tsv")
rsr_tsv_path <- file.path(cmat_dir_hg, "aggregate_cormat_RSR_hg.tsv")
na_tsv_path <- file.path(cmat_dir_hg, "aggregate_cormat_NA_hg.tsv")



prepare_celltype_mat <- function(id, mat, meta, pc_df, cell_type, min_cell = 20) {
  
  stopifnot(inherits(mat, "dgCMatrix"),
            c("ID", "Cell_type") %in% colnames(meta),
            cell_type %in% meta[["Cell_type"]],
            identical(rownames(mat), pc_df$Symbol))
  
  # ids <- meta[meta$Cell_type %in% cell_type, "ID"][["ID"]]
  ids <- filter(meta, Cell_type %in% cell_type) %>% pull(ID)
  
  if (!all(ids %in% colnames(mat))) {
    message(paste(id, "did not have all meta IDs in mat column names"))
    ids <- intersect(ids, colnames(mat))
  }
  
  ct_mat <- t(mat[pc_df$Symbol, ids])
  ct_mat <- zero_sparse_cols(ct_mat, min_cell)
  
  return(ct_mat)
}




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
    agg_method = "RSR"
  )
)



# FZ matrix as .tsv files
if (!file.exists(fz_tsv_path)) {
  
  agg_l <- readRDS(fz_path)
  fwrite_mat(aggr_l$Agg_mat, fz_tsv_path)
  
  if (!file.exists(na_tsv_path)) {
    fwrite_mat(aggr_l$NA_mat, na_tsv_path)
  }

}



# RSR matrix as .tsv files
if (!file.exists(fz_tsv_path)) {
  
  agg_l <- readRDS(fz_path)
  fwrite_mat(aggr_l$Agg_mat, fz_tsv_path)
  
  if (!file.exists(na_tsv_path)) {
    fwrite_mat(aggr_l$NA_mat, na_tsv_path)
  }
  
}
