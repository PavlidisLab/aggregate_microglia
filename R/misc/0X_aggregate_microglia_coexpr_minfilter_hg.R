## 
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
# library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/dev_aggregate_functions.R")

# Gene table and dataset meta
pc_df <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
meta <- read.delim(mcg_meta_path) %>% filter(Species == "Human")

# Paths for the Fisher's Z, Rank Sum Rank, and NA tracking matrices
fz_path <- file.path(cmat_dir_hg, "aggregate_cormat_FZ_minfilter_hg.RDS")
rsr_path <- file.path(cmat_dir_hg, "aggregate_cormat_RSR_minfilter_hg.RDS")

# Only consider genes that are minimally measured across the collection

count_summ <- readRDS(count_summ_path)$Human

keep_genes <- count_summ$Summ_df %>% 
  filter(N_msr >= floor(max(N_msr) * 1/3)) %>% 
  pull(Symbol)

pc_df <- filter(pc_df, Symbol %in% keep_genes)



# Generate and save
# ------------------------------------------------------------------------------

# TODO: remove, just ensuring it matches full matrix
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
