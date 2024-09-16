library(tidyverse)
library(data.table)
source("R/00_config.R")
source("R/utils/dev_functions.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# List of cell types
ct_l <- readRDS(celltype_list_path)

pc_df <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# Save out
out_dir <- "/space/scratch/amorin/TR_singlecell/Microglia/Mm_pcor"
dir.create(out_dir, showWarnings = FALSE)
aggr_path <- file.path(out_dir, "test_aggregate_matrix_allrank.tsv")
na_path <- file.path(out_dir, "test_NA_matrix_allrank.tsv")



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


ct_df <- filter(ct_df, Species == "Mouse")

# Add # GSE118020: Micro (micro grep grabs non-microglia) 
ct_df <- rbind(ct_df, c("GSE118020", "Micro", 3230, "Mouse"))


ct_df$Path <- paste0(amat_dir, ct_df$ID, "/",  ct_df$ID, "_clean_mat_and_meta_CPM.RDS")




# Generate and save
# ------------------------------------------------------------------------------



aggr_l <- aggr_coexpr_across_datasets(
  input_df = ct_df,
  pc_df = pc_df,
  cor_method = "pearson",
  agg_method = "allrank"
)



fwrite(
  data.frame(aggr_l$Agg_mat, check.names = FALSE),
  sep = "\t",
  row.names = TRUE,
  quote = FALSE,
  verbose = FALSE,
  showProgress = FALSE,
  file = aggr_path
)



fwrite(
  data.frame(aggr_l$NA_mat, check.names = FALSE),
  sep = "\t",
  row.names = TRUE,
  quote = FALSE,
  verbose = FALSE,
  showProgress = FALSE,
  file = na_path
)
