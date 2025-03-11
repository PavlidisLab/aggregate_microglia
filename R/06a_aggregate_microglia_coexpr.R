## Save aggregate coexpression for microglia datasets in human and mouse. This
## script saves out 3 forms of aggregation for each species: average Fisher's Z
## (FZ), rank sum rank for jointly ranking the full matrix ("allrank"), and
## rank sum rank for column/gene-wise ranking ("colrank"). In addition, for the
## allrank aggregation, it saves out a version that only considers filtered
## genes measured in microglia. This was done to show effect of measurement on
## ranking procedure. As FZ is not ranked, it was not ran with the subset genes.
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

# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# Filter: require measured in at least a third of datasets
keep_hg <- count_summ$Human$Filter_genes
keep_mm <- count_summ$Mouse$Filter_genes
pc_filt_hg <- filter(pc_hg, Symbol %in% keep_hg)
pc_filt_mm <- filter(pc_mm, Symbol %in% keep_mm)



# Generate and save
# ------------------------------------------------------------------------------


# Human FZ (full protein coding table)
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


# Human allrank (full protein coding table)
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


# Human colrank (full protein coding table)
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


# Human allrank (filtered protein coding table)
save_function_results(
  path = mcg_allrank_filt_hg_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_hg,
    pc_df = pc_filt_hg,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)




# Mouse FZ (full protein coding table)
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


# Mouse allrank (full protein coding table)
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


# Mouse colrank (full protein coding table)
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


# Mouse allrank (filter protein coding table)
save_function_results(
  path = mcg_allrank_filt_mm_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_filt_mm,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)
