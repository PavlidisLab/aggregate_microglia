## TODO:
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_df <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
meta <- read.delim(macro_meta_path) %>% filter(Species == "Human")

# Paths for the Fisher's Z, Rank Sum Rank, and NA tracking matrices
fz_path <- "/space/scratch/amorin/aggregate_microglia/Cormats/Macrophage_hg/aggregate_cormat_FZ_macrophage_hg.RDS"
rsr_path <- "/space/scratch/amorin/aggregate_microglia/Cormats/Macrophage_hg/aggregate_cormat_RSR_macrophage_hg.RDS"



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
  path = rsr_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)
