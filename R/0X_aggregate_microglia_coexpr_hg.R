## TODO: 
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_df <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
meta <- read.delim(mcg_meta_path) %>% filter(Species == "Human")

# Paths for the lists of aggregate matrices
fz_path <- file.path(cmat_dir_hg, "aggregate_cormat_FZ_hg.RDS")
allrank_path <- file.path(cmat_dir_hg, "aggregate_cormat_allrank_hg.RDS")
colrank_path <- file.path(cmat_dir_hg, "aggregate_cormat_colrank_hg.RDS")


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
  path = allrank_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)



save_function_results(
  path = colrank_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df,
    cor_method = "pearson",
    agg_method = "colrank"
  )
)
