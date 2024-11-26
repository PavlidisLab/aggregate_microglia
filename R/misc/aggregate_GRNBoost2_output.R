## GRNBoost2 was ran 100 times on each microglia dataset. For each dataset ID,
## load each iterations and average into one matrix. Then, generate one global
## matrix that averages all dataset matrices.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/grnboost2_utils.R")

grn_dir <- file.path(data_out_dir, "GRNBoost2")

mcg_meta <- read.delim(mcg_meta_path)
ids_hg <- filter(mcg_meta, Species == "Human")$ID

# Load genes/TFs
tf_hg <- read.delim("/home/amorin/Data/Metadata/TFs_human.tsv")
pc_hg <- read.delim(ref_hg_path)


# TODO: remove/reduce when everything is ran
# GRN output for each dataset
grn_files_list <- lapply(ids_hg, function(x) {
  list.files(file.path(grn_dir, x), full.names = TRUE, pattern = ".*_iter.*")
})

names(grn_files_list) <- ids_hg
grn_files_list <- grn_files_list[sapply(grn_files_list, length) > 0]
ids_hg <- names(grn_files_list)


avg_grn_l <- average_and_save_each_network(ids = ids_hg,
                                           grn_dir = grn_dir,
                                           ncore = ncore)

all_avg_mat <- average_all_networks(avg_grn_l)
