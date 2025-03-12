## Saving aggregate macrophage coexpression for human and mouse. Here, just
## allrank and FZ.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
meta <- read.delim(macro_meta_path)
meta_hg <- filter(meta, Species == "Human")
meta_mm <- filter(meta, Species == "Mouse")



# Generate and save
# ------------------------------------------------------------------------------


# Human allrank
save_function_results(
  path = macro_allrank_hg_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_hg,
    pc_df = pc_hg,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)


# Human FZ
save_function_results(
  path = macro_fz_hg_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_hg,
    pc_df = pc_hg,
    cor_method = "pearson",
    agg_method = "FZ"
  )
)



# Mouse allrank
save_function_results(
  path = macro_allrank_mm_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)


# Mouse FZ
save_function_results(
  path = macro_fz_mm_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = "FZ"
  )
)
