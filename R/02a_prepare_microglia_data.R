## Organizing metadata and the underlying count matrices for microglia
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# List of cell types
ct_list <- readRDS(celltype_list_path)

# Based on manual inspection of the cell type list
ct_str <- ".*microglia.*|^microg.*|^mgs$|^mg$|^micro$"

# Create a metadata table of the relevant datasets
ct_df <- create_celltype_df(pattern = ct_str, ct_list, sc_meta)
ct_df$Path <- paste0(amat_dir, ct_df$ID, "/",  ct_df$ID, "_clean_mat_and_meta_CPM.RDS")




write.table(ct_df, sep = "\t", quote = FALSE, row.names = FALSE, 
            file = mcg_meta_path)