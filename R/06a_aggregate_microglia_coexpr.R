## Save aggregate coexpression for microglia datasets in human and mouse. This
## script saves out 3 forms of aggregation for each species: average Fisher's Z
## (FZ), rank sum rank for jointly ranking the full matrix ("allrank"), and
## rank sum rank for column/gene-wise ranking ("colrank")
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
meta <- read.delim(mcg_meta_path)
meta_hg <- filter(meta, Species == "Human")
meta_mm <- filter(meta, Species == "Mouse")



# Generate and save
# ------------------------------------------------------------------------------


# Human FZ
save_function_results(
  path = mcg_fz_hg_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_hg,
    pc_df = pc_hg,
    cor_method = "pearson",
    agg_method = "FZ"
  )
)


# Human allrank
save_function_results(
  path = mcg_allrank_hg_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_hg,
    pc_df = pc_hg,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)


# Human colrank
save_function_results(
  path = mcg_colrank_hg_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_hg,
    pc_df = pc_hg,
    cor_method = "pearson",
    agg_method = "colrank"
  )
)



# Mouse FZ
save_function_results(
  path = mcg_fz_mm_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = "FZ"
  )
)


# Mouse allrank
save_function_results(
  path = mcg_allrank_mm_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)


# Mouse colrank
save_function_results(
  path = mcg_colrank_mm_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = "colrank"
  )
)
