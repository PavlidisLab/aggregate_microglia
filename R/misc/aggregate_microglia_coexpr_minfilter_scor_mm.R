## Generating rank-based aggregate matrices after filtering for genes that
## are measured in a minimal amount of datasets. Ranking procedure (unlike FZ)
## is affected by count of non-measured genes (treated as NA->0->ties)
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
# library(aggtools)
source("R/utils/aggregate_scor.R")
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_df <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
meta <- read.delim(mcg_meta_path) %>% filter(Species == "Mouse")

# Paths for the aggregate matrix lists
allrank_path <- file.path(cmat_dir_mm, "aggregate_cormat_allrank_filter_scor_mm.RDS")
fz_path <- file.path(cmat_dir_mm, "aggregate_cormat_FZ_filter_scor_mm.RDS")


# Only consider genes that are minimally measured across the collection
count_summ <- readRDS(mcg_count_summ_list_path)

# Filter: require measured in at least a third of datasets
keep_genes <- count_summ$Mouse$Summ_df %>% 
  filter(N_msr >= floor(max(N_msr) * 1/3)) %>%
  pull(Symbol)


pc_df <- filter(pc_df, Symbol %in% keep_genes)



# Generate and save
# ------------------------------------------------------------------------------


save_function_results(
  path = allrank_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df,
    cor_method = "spearman",
    agg_method = "allrank"
  )
)


save_function_results(
  path = fz_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df,
    cor_method = "spearman",
    agg_method = "FZ"
  )
)

