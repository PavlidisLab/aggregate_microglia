## Generating rank-based aggregate matrices after filtering for genes that
## are measured in a minimal amount of datasets. Ranking procedure (unlike FZ)
## is affected by count of non-measured genes (treated as NA->0->ties)
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
# library(aggtools)
source("/home/amorin/Projects/aggtools/R/aggregate_functions.R")
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_df <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
meta <- read.delim(mcg_meta_path) %>% filter(Species == "Human")

# Paths for the aggregate matrix lists
allrank_len_path <- file.path(cmat_dir_hg, "aggregate_cormat_allrank_filter_lenient_hg.RDS")
allrank_str_path <- file.path(cmat_dir_hg, "aggregate_cormat_allrank_filter_stringent_hg.RDS")
colrank_len_path <- file.path(cmat_dir_hg, "aggregate_cormat_colrank_filter_lenient_hg.RDS")
colrank_str_path <- file.path(cmat_dir_hg, "aggregate_cormat_colrank_filter_stringent_hg.RDS")


# Only consider genes that are minimally measured across the collection
count_summ <- readRDS(mcg_count_summ_list_path)

# Lenient: require measured in at least a third of datasets
keep_len <- count_summ$Human$Summ_df %>% 
  filter(N_msr >= floor(max(N_msr) * 1/3)) %>%
  pull(Symbol)

# Stringent: At least 90% of datasets
keep_str <- count_summ$Human$Summ_df %>% 
  filter(N_msr >= floor(max(N_msr) * 0.9)) %>% 
  pull(Symbol)


pc_df_len <- filter(pc_df, Symbol %in% keep_len)
pc_df_str <- filter(pc_df, Symbol %in% keep_str)



# Generate and save
# ------------------------------------------------------------------------------


save_function_results(
  path = allrank_len_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df_len,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)


save_function_results(
  path = allrank_str_path,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta,
    pc_df = pc_df_str,
    cor_method = "pearson",
    agg_method = "allrank"
  )
)

