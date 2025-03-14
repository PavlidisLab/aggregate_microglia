## Basic summarization of individual coexpr mats. Original idea was to explore
## the most variable gene coexpression pairs across datasets; didn't get far
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
library(pheatmap)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table
pc_df_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_df_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# Paths of individual coexpression matrices
cor_paths_hg <- list.files(cmat_dir_mcg_hg, full.names = TRUE, pattern = "cormat.tsv")
cor_paths_mm <- list.files(cmat_dir_mcg_mm, full.names = TRUE, pattern = "cormat.tsv")

# Output summary objects
outfile_hg <- file.path(cmat_dir_mcg_hg, "raw_cor_avg_sd_msr_cv_list_hg.RDS")
outfile_mm <- file.path(cmat_dir_mcg_mm, "raw_cor_avg_sd_msr_cv_list_mm.RDS")

# Microglia meta with file paths
mcg_meta <- read.delim(mcg_meta_path) %>% distinct(ID, .keep_all = TRUE)
meta_hg <- filter(mcg_meta, Species == "Human")
meta_mm <- filter(mcg_meta, Species == "Mouse")


# Functions
# ------------------------------------------------------------------------------


# Iteratively load each coexpr mat for averaging. Same idea as aggregation, minus
# FZ or ranking

calc_avg_and_na <- function(paths, pc_df, verbose = TRUE) {
  
  avg_mat <- msr_mat <- init_agg_mat(pc_df)
  
  for (file in paths) {
    message(paste(file, Sys.time()))
    cmat <- fread_to_mat(file, genes = pc_df$Symbol)
    non_na <- !is.na(cmat)
    avg_mat[non_na] <- avg_mat[non_na] + cmat[non_na]
    msr_mat[non_na] <- msr_mat[non_na] + 1
    rm(cmat)
    gc(verbose = FALSE)
  }
  
  avg_mat <- avg_mat / msr_mat
  
  return(list(Avg_mat = avg_mat, Msr_mat = msr_mat))
  
}


# Iteratively load each coexpr mat for calculating std. dev.

calc_sd <- function(paths, pc_df, avg_l, verbose = TRUE) {
  
  var_mat <- init_agg_mat(pc_df)
  avg_mat <- avg_l$Avg_mat
  msr_mat <- avg_l$Msr_mat
  
  for (file in paths) {
    message(paste(file, Sys.time()))
    cmat <- fread_to_mat(file, genes = pc_df$Symbol)
    non_na <- !is.na(cmat)
    diff <- (cmat[non_na] - avg_mat[non_na])^2
    var_mat[non_na] <- var_mat[non_na] + diff
    rm(cmat)
    gc(verbose = FALSE)
  }
  
  sd_mat <- sqrt(var_mat / msr_mat)
  
  return(sd_mat)
}


# Calc avg, SD, and CV for each gene pair across datasets

prepare_summary_l <- function(paths, pc_df, verbose = TRUE) {
  
  avg_l <- calc_avg_and_na(paths = paths, pc_df = pc_df)
  sd_mat <- calc_sd(paths = paths, pc_df = pc_df, avg_l = avg_l)
  summ_l <- avg_l
  summ_l$SD_mat <- sd_mat
  summ_l$CV_mat <- summ_l$SD_mat / summ_l$Avg_mat
  
  return(summ_l)
}



# Run/save/load
# ------------------------------------------------------------------------------


save_function_results(
  path = outfile_hg,
  fun = prepare_summary_l,
  args = list(paths = cor_paths_hg, pc_df = pc_df_hg)
)
  


save_function_results(
  path = outfile_mm,
  fun = prepare_summary_l,
  args = list(paths = cor_paths_mm, pc_df = pc_df_mm)
)



summ_hg <- readRDS(outfile_hg)
summ_mm <- readRDS(outfile_mm)



# Basic exploration
# ------------------------------------------------------------------------------


summ_df_hg <- cbind(
  mat_to_df(summ_hg$Avg_mat, value_name = "Avg"),
  SD = mat_to_df(summ_hg$SD_mat, value_name = "SD")[["SD"]],
  CV =  mat_to_df(summ_hg$CV_mat, value_name = "CV")[["CV"]],
  N_msr = mat_to_df(summ_hg$Msr_mat, value_name = "N_msr")[["N_msr"]]
) %>%
  arrange(desc(Avg))


summ_df_mm <- cbind(
  mat_to_df(summ_mm$Avg_mat, value_name = "Avg"),
  SD = mat_to_df(summ_mm$SD_mat, value_name = "SD")[["SD"]],
  CV =  mat_to_df(summ_mm$CV_mat, value_name = "CV")[["CV"]],
  N_msr = mat_to_df(summ_mm$Msr_mat, value_name = "N_msr")[["N_msr"]]
) %>%
  arrange(desc(Avg))




top_avg_hg <- summ_df_hg %>% 
  filter(N_msr >= floor(nrow(meta_hg) * min_msr_frac)) %>% 
  slice_max(abs(Avg), n = 10e3)

top_sd_hg <- summ_df_hg %>% 
  filter(N_msr >= floor(nrow(meta_hg) * min_msr_frac)) %>% 
  slice_max(SD, n = 10e3)

top_cv_hg <- summ_df_hg %>% 
  filter(N_msr >= floor(nrow(meta_hg) * min_msr_frac)) %>% 
  slice_max(CV, n = 10e3)



top_avg_mm <- summ_df_mm %>% 
  filter(N_msr >= floor(nrow(meta_mm) * min_msr_frac)) %>% 
  slice_max(abs(Avg), n = 10e3)

top_sd_mm <- summ_df_mm %>% 
  filter(N_msr >= floor(nrow(meta_mm) * min_msr_frac)) %>% 
  slice_max(SD, n = 10e3)

top_cv_mm <- summ_df_mm %>% 
  filter(N_msr >= floor(nrow(meta_mm) * min_msr_frac)) %>% 
  slice_max(CV, n = 10e3)
