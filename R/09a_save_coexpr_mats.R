## To perform differential coexpression, need to have the individual coexpr
## mats. This performs the same filtering as the aggregation workflow, but 
## saves the coexpr matrices for microglia and macrophages instead of aggregating.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene tables
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# Meta with file paths
mcg_meta <- read.delim(mcg_meta_dedup_path)
mcg_meta_hg <- filter(mcg_meta, Species == "Human")
mcg_meta_mm <- filter(mcg_meta, Species == "Mouse")

macro_meta <- read.delim(macro_meta_dedup_path)
macro_meta_hg <- filter(macro_meta, Species == "Human")
macro_meta_mm <- filter(macro_meta, Species == "Mouse")


# Iterate through each dataset in the input/meta df, loading the count matrix,
# performing coexpression, and saving out the result
# ------------------------------------------------------------------------------

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
    ct <- unlist(str_split(ct, ";"))
    
    dat_path <- unique(input_df[input_df$ID == id, "Path"])
    
    # Load dataset and get count matrix for current cell type
    
    dat <- load_scdat(dat_path)
    
    ct_mat <- prepare_celltype_mat(mat = dat[["Mat"]],
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





# Run/save
# ------------------------------------------------------------------------------


# Microglia human
save_coexpr_multi_dataset(
  input_df = mcg_meta_hg,
  pc_df = pc_hg,
  cor_method = "pearson",
  out_dir = cmat_dir_hg
)


# Microglia mouse
save_coexpr_multi_dataset(
  input_df = mcg_meta_mm,
  pc_df = pc_mm,
  cor_method = "pearson",
  out_dir = cmat_dir_mm
)



# Macrophage human
save_coexpr_multi_dataset(
  input_df = macro_meta_hg,
  pc_df = pc_hg,
  cor_method = "pearson",
  out_dir = "/space/scratch/amorin/aggregate_microglia/Cormats/Macrophage_hg"
)


# Macrophage mouse
save_coexpr_multi_dataset(
  input_df = macro_meta_mm,
  pc_df = pc_mm,
  cor_method = "pearson",
  out_dir = "/space/scratch/amorin/aggregate_microglia/Cormats/Macrophage_mm"
)
