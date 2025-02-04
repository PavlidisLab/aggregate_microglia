## GRNBoost2 was ran 100 times on each microglia dataset. For each dataset ID,
## load each iterations and average into one matrix. Then, generate one global
## matrix that averages all dataset matrices.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("/home/amorin/Projects/aggtools/R/aggregate_functions.R") # TODO: replace
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/grnboost2_utils.R")

# Directory of GRNBoost2 output
grn_dir <- file.path(data_out_dir, "GRNBoost2")

# Output paths
grn_avg_hg_path <- file.path(data_out_dir, "GRNBoost2_average_mat_hg.RDS")
grn_avg_mm_path <- file.path(data_out_dir, "GRNBoost2_average_mat_mm.RDS")
grn_allrank_hg_path <- file.path(data_out_dir, "GRNBoost2_allrank_mat_hg.RDS")
grn_allrank_mm_path <- file.path(data_out_dir, "GRNBoost2_allrank_mat_mm.RDS")

# Microglia datasets
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_mm <- filter(mcg_meta, Species == "Mouse")$ID
ids_hg <- filter(mcg_meta, Species == "Human")$ID

# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# Load genes/TFs
pc_hg <- read.delim(ref_hg_path)
pc_mm <- read.delim(ref_mm_path)
tfs_hg <- read.delim(tfs_hg_path)
tfs_mm <- read.delim(tfs_mm_path)


# Min genes to keep
keep_hg <- count_summ$Human$Summ_df %>% 
  filter(N_msr >= floor(max(N_msr) * 1/3)) %>%
  select(Symbol)

keep_mm <- count_summ$Mouse$Summ_df %>% 
  filter(N_msr >= floor(max(N_msr) * 1/3)) %>%
  select(Symbol)

keep_tfs_hg <- intersect(tfs_hg$Symbol, keep_hg$Symbol)
keep_tfs_mm <- intersect(tfs_mm$Symbol, keep_mm$Symbol)



# For each dataset, average the GRNBoost2 iterations into one matrix
grn_l_hg <- average_and_save_each_network(ids = ids_hg,
                                          grn_dir = grn_dir,
                                          ncore = ncore)

grn_l_mm <- average_and_save_each_network(ids = ids_mm,
                                          grn_dir = grn_dir,
                                          ncore = ncore)


# Average each dataset's matrix into one averaged matrix of importance scores
grn_avg_hg <- average_all_networks(grn_l_hg, keep_hg$Symbol, keep_tfs_hg)
grn_avg_mm <- average_all_networks(grn_l_mm, keep_mm$Symbol, keep_tfs_mm)


# Aggregate ranks for each dataset
grn_allrank_hg <- rank_aggregate_grn(grn_l_hg, keep_hg$Symbol, keep_tfs_hg)
grn_allrank_mm <- rank_aggregate_grn(grn_l_mm, keep_mm$Symbol, keep_tfs_mm)


# Save out
saveRDS(grn_avg_hg, grn_avg_hg_path)
saveRDS(grn_avg_mm, grn_avg_mm_path)
saveRDS(grn_allrank_hg, grn_allrank_hg_path)
saveRDS(grn_allrank_mm, grn_allrank_mm_path)
