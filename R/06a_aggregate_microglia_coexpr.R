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

# Paths for the lists of aggregate coexpression and NA tracking matrices
mcg_fz_hg_path <- file.path(cmat_dir_hg, "aggregate_cormat_FZ_hg.RDS")
mcg_allrank_hg_path <- file.path(cmat_dir_hg, "aggregate_cormat_allrank_hg.RDS")
mcg_colrank_hg_path <- file.path(cmat_dir_hg, "aggregate_cormat_colrank_hg.RDS")

mcg_fz_mm_path <- file.path(cmat_dir_mm, "aggregate_cormat_FZ_mm.RDS")
mcg_allrank_mm_path <- file.path(cmat_dir_mm, "aggregate_cormat_allrank_mm.RDS")
mcg_colrank_mm_path <- file.path(cmat_dir_mm, "aggregate_cormat_colrank_mm.RDS")


# Generate and save
# ------------------------------------------------------------------------------


save_function_results(
  path = fz_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_hg,
    cor_method = "pearson",
    agg_method = "FZ"
  )
)



save_function_results(
  path = allrank_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_hg,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)



save_function_results(
  path = colrank_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_hg,
    cor_method = "pearson",
    agg_method = "colrank"
  )
)
