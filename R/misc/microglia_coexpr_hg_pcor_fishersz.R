library(WGCNA)
library(tidyverse)
library(data.table)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/cell_type_aggr_coexpr.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# List of cell types
ct_l <- readRDS(celltype_list_path)

pc_df <- read.delim(ref_hg_path, stringsAsFactors = FALSE)

# Save out
out_dir <- "/space/scratch/amorin/TR_singlecell/Microglia/Hg_pcor"
dir.create(out_dir, showWarnings = FALSE)
aggr_path <- file.path(out_dir, "aggregate_matrix_fishersz.tsv")
na_path <- file.path(out_dir, "NA_matrix.tsv")



# Extract IDs of specific cell type of interest
# ------------------------------------------------------------------------------


ct <- "microglia|mcg|mgl|microg"  


ct_df <- lapply(ct_l, function(x) {
  ct_vec <- str_to_lower(x$Ct_count$Cell_type)
  ct_which <- str_detect(ct_vec, ct)
  
  if (sum(ct_which) == 0) {
    return(NA)
  }
  
  x$Ct_count[ct_which, ]
  
})


ct_df <- ct_df[!is.na(ct_df)]



ct_df <- data.frame(
  do.call(rbind, ct_df)
) %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(ID = str_replace(ID, "\\.[:digit:]+$", "")) %>% 
  left_join(., sc_meta[, c("ID", "Species")], by = "ID")


ct_df <- filter(ct_df, Species == "Human")



# Generate and save
# ------------------------------------------------------------------------------



# Check if there are IDs without cormat, prompting regen of aggregate

all_ids <- unique(ct_df[["ID"]])
cmat_paths <- list.files(path = out_dir, pattern = "cormat")
cmat_ids <- str_replace(cmat_paths, "_cormat.tsv", "")
all_ids_have_cor <- all(all_ids %in% cmat_ids)



if (!file.exists(aggr_path) || !all_ids_have_cor) {
  
  
  aggr_l <- aggr_celltype_coexpr_fishersZ(ct_df = ct_df,
                                          pc_df = pc_df,
                                          cor_method = "pearson",
                                          ncores = ncore,
                                          in_dir = amat_dir,
                                          out_dir = out_dir)
  
  fwrite(
    data.frame(aggr_l$Avg_nonNA, check.names = FALSE),
    sep = "\t",
    row.names = TRUE,
    quote = FALSE,
    verbose = FALSE,
    showProgress = FALSE,
    file = aggr_path
  )
  
  
  # fwrite(
  #   data.frame(aggr_l$NA_mat, check.names = FALSE),
  #   sep = "\t",
  #   row.names = TRUE,
  #   quote = FALSE,
  #   verbose = FALSE,
  #   showProgress = FALSE,
  #   file = na_path
  # )
}
