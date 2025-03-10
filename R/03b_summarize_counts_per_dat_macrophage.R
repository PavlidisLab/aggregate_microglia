## Summarize the CPM-normalized count data across macrophage datasets, saving
## this info out as a list.
## -----------------------------------------------------------------------------

library(tidyverse)
library(preprocessCore)
library(matrixStats)
source("R/00_config.R")
source("R/utils/functions.R")

# Dataset meta and human/mouse data IDs
macro_meta <- read.delim(macro_meta_dedup_path)
ids_hg <- unique(filter(macro_meta, Species == "Human")$ID)
ids_mm <- unique(filter(macro_meta, Species == "Mouse")$ID)


if (!file.exists(macro_count_summ_list_path)) {
  
  # List of count matrices and their metadata
  dat_l <- readRDS(macro_dat_path)
  
  count_summ_hg <- summarize_gene_counts(dat_l[ids_hg])
  count_summ_mm <- summarize_gene_counts(dat_l[ids_mm])
  
  saveRDS(list(Human = count_summ_hg, Mouse = count_summ_mm), macro_count_summ_list_path)

}
