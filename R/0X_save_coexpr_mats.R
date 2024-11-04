library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table
pc_df_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_df_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)





# Microglia meta with file paths
mcg_meta <- read.delim(mcg_meta_path) %>% distinct(ID, .keep_all = TRUE)
meta_hg <- filter(mcg_meta, Species == "Human")
meta_mm <- filter(mcg_meta, Species == "Mouse")


##


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



save_coexpr_multi_dataset <- function(input_df,
                                      pc_df,
                                      cor_method,
                                      min_cell = 20,
                                      verbose = TRUE,
                                      out_dir) {
  
  stopifnot(all(c("ID", "Cell_type", "Path") %in% colnames(input_df)),
            cor_method %in% c("pearson", "spearman"))
  
  data_ids <- unique(input_df[["ID"]])
  
  for (id in data_ids) {
    
    file <- file.path(out_dir, paste0(id, "_cormat.tsv"))
    
    if (verbose) message(paste(id, Sys.time()))
    if (file.exists(file)) next()
    
    ct <- input_df[input_df$ID == id, "Cell_type"]
    dat_path <- unique(input_df[input_df$ID == id, "Path"])
    
    # Load dataset and get count matrix for current cell type
    
    dat <- load_scdat(dat_path)
    
    ct_mat <- prepare_celltype_mat(id = id,
                                   mat = dat[["Mat"]],
                                   meta = dat[["Meta"]],
                                   pc_df = pc_df,
                                   cell_type = ct,
                                   min_cell = min_cell)
    
    # Check if filtering removed all genes
    
    no_msr <- all(ct_mat == 0)
    
    if (no_msr) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    # Get cell-type cor matrix and increment count of NAs before imputing to 0
    
    cmat <- calc_sparse_correlation(ct_mat, cor_method)
    
    # Transform raw correlation matrix, add to aggregate and clean up
    
    fwrite_mat(cmat, file)
    rm(cmat, ct_mat)
    gc(verbose = FALSE)
    
  }
  
  return(NULL)
}


##





save_coexpr_multi_dataset(
  input_df = meta_hg,
  pc_df = pc_df_hg,
  cor_method = "pearson",
  out_dir = cmat_dir_hg
)


save_coexpr_multi_dataset(
  input_df = meta_mm,
  pc_df = pc_df_mm,
  cor_method = "pearson",
  out_dir = cmat_dir_mm
)

