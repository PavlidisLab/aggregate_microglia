## Here, generating aggregate coexpression for mouse microglia (using filtered
## gene set) to compare against Pearson's cor.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
source("R/utils/functions.R")
source("R/utils/aggregate_scor.R") # NOTE: patch script to allow spearman cor
source("R/00_config.R")

# Gene table and dataset meta
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
meta_mm <- read.delim(mcg_meta_path) %>% filter(Species == "Mouse")

# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# Filter: require measured in at least a third of datasets
keep_mm <- count_summ$Mouse$Filter_genes
pc_filt_mm <- filter(pc_mm, Symbol %in% keep_mm)



# Generate and save
# ------------------------------------------------------------------------------


save_function_results(
  path = mcg_allrank_scor_mm_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_filt_mm,
    cor_method = "spearman",
    agg_method = "allrank"
  )
)


save_function_results(
  path = mcg_fz_scor_mm_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_filt_mm,
    cor_method = "spearman",
    agg_method = "FZ"
  )
)
