## TODO
## https://journals.asm.org/doi/10.1128/microbiolspec.mchd-0015-2015
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")


# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# List of cell types
celltype_l <- readRDS(celltype_list_path)



# 
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



all_ct_df <- do.call(rbind, lapply(names(celltype_l), function(x) {
  data.frame(celltype_l[[x]]$Ct_count, ID = x)
}))


tally_ct <- count(all_ct_df, Cell_type)



check <- c(
  "Macro",  # keep: GSE134003, GSE217368, and GSE183310
  "MAC",  # keep: GSE222510
  "Mac",  # keep: GSE198582
  "Mac-I",  # keep: GSE198582 (Mac-I to V and Mac-a to e)
  "Mac_",  # keep: GSE198582
  "Mon/Mac",  # rm: unclear if distinct from monocytes
  "Monocyte/Macrophage",  # rm: unclear if distinct from monocytes
  "endo-Mac",  # keep: GSE198582
  "epi-Mac",  # keep: GSE198582
  "macrophage dendritic cell progenitor",  # rm: unclear if distinct from monocytes
  "Macrophages_progenitors",  # rm: unclear if distinct from monocytes
  "MonocyteMacrophages",  # rm: unclear if distinct from monocytes
  "prol.Mac",  # keep: GSE198582
  "recMac",  # keep: GSE225664 (recruited macrophage)
  "IMs",  # keep: GSE225667, GSE225664, GSE225662, GSE225666, GSE225663, GSE225665 (interstitial macrophages)
  "Macrophage/cDC"  # rm: unclear if distinct from conventional dendritic cells
)


check_df <- filter(all_ct_df, Cell_type %in% check)


# This was based off of manual inspection of the cell type list
ct_str <- ".*macrophage.*|.*kupffer.*|^macro$|^mac$|mac-.*|mac_.*|.*-mac$|^recmac$|^ims$"


rm_str <- c("macrophage dendritic cell progenitor",
            "MonocyteMacrophages",
            "Mon/Mac",
            "Monocyte/Macrophage",
            "Macrophages_progenitors")


ct_df <- create_celltype_df(celltype_l, sc_meta, ct_str)

ct_df <- filter(ct_df, Cell_type %!in% rm_str)

ct_df$Path <- paste0(amat_dir, ct_df$ID, "/",  ct_df$ID, "_clean_mat_and_meta_CPM.RDS")

write.table(ct_df, sep = "\t", quote = FALSE, row.names = FALSE, file = macro_meta_path)
