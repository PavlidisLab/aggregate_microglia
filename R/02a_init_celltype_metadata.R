## Organizing tables corresponding to microglia cell types
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")


# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# List of cell types
celltype_l <- readRDS(celltype_list_path)


# Extract datasets and fuzzy cell type names corresponding to microglia
# ------------------------------------------------------------------------------


create_celltype_df <- function(celltype_l, sc_meta, str_vec) {
  
  # Return list of relevant cell types matching string search
  ct_l <- lapply(celltype_l, function(x) {
    ct_vec <- str_to_lower(x$Ct_count$Cell_type)
    ct_which <- str_detect(ct_vec, str_vec)
    
    if (sum(ct_which) == 0) {
      return(NA)
    }
    
    x$Ct_count[ct_which, ]
    
  })
  
  ct_l <- ct_l[!is.na(ct_l)]
  
  # Bind into a df and join with relevant details from meta
  ct_df <- do.call(rbind, ct_l) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "ID") %>% 
    mutate(ID = str_replace(ID, "\\.[:digit:]+$", "")) %>% 
    left_join(., 
              sc_meta[, c("ID", "Species", "Platform", "GEO_link", "Data_link")], 
              by = "ID")
  
  return(ct_df)
}




# This was based off of manual inspection of the cell type list
mcg_str <- ".*microglia.*|^microg.*|^mgs$|^mg$|^micro$"


mcg_df <- create_celltype_df(celltype_l, sc_meta, mcg_str)


mcg_df$Path <- paste0(amat_dir, mcg_df$ID, "/",  mcg_df$ID, "_clean_mat_and_meta_CPM.RDS")


write.table(mcg_df, sep = "\t", quote = FALSE, row.names = FALSE, 
            file = mcg_meta_path)
