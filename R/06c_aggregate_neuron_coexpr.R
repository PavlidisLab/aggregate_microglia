## Saving out neuron aggregates
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table and dataset meta
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

meta <- read.delim(neuron_meta_path)
meta_hg <- filter(meta, Species == "Human")
meta_mm <- filter(meta, Species == "Mouse")

# Paths for the Fisher's Z, Rank Sum Rank, and NA tracking matrices
fz_path_hg <- "/space/scratch/amorin/aggregate_microglia/Cormats/Neuron_hg/aggregate_cormat_FZ_neuron_hg.RDS"
fz_path_mm <- "/space/scratch/amorin/aggregate_microglia/Cormats/Neuron_mm/aggregate_cormat_FZ_neuron_mm.RDS"

allrank_path_hg <- "/space/scratch/amorin/aggregate_microglia/Cormats/Neuron_hg/aggregate_cormat_allrank_neuron_hg.RDS"
allrank_path_mm <- "/space/scratch/amorin/aggregate_microglia/Cormats/Neuron_mm/aggregate_cormat_allrank_neuron_mm.RDS"



# Generate and save
# ------------------------------------------------------------------------------


# Human allrank
save_function_results(
  path = allrank_path_hg,
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
  path = fz_path_hg,
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
  path = allrank_path_mm,
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
  path = fz_path_mm,
  fun = aggr_coexpr_multi_dataset,
  args = list(
    input_df = meta_mm,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = "FZ"
  )
)
