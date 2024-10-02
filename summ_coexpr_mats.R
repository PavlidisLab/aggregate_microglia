library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Gene table
pc_df_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_df_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# 
dir_hg <- "/space/scratch/amorin/TR_singlecell/Microglia/Hg_pcor_test/"
dir_mm <- "/space/scratch/amorin/TR_singlecell/Microglia/Mm_pcor_test/"

cor_paths_hg <- list.files(dir_hg, full.names = TRUE, pattern = "cormat.tsv")
cor_paths_mm <- list.files(dir_mm, full.names = TRUE, pattern = "cormat.tsv")

outfile_hg <- file.path(dir_hg, "raw_cor_avg_sd_msr_list_hg.RDS")
outfile_mm <- file.path(dir_mm, "raw_cor_avg_sd_msr_list_mm.RDS")


# Microglia meta with file paths
mcg_meta <- read.delim(mcg_meta_path) %>% distinct(ID, .keep_all = TRUE)
meta_hg <- filter(mcg_meta, Species == "Human")
meta_mm <- filter(mcg_meta, Species == "Mouse")


# hg_l <- lapply(paths_hg[1], fread_to_mat, genes = pc_df_hg$Symbol)
# # hg_l[[2]] <- hg_l[[1]]
# msr_hg <- Reduce("+", lapply(hg_l, function(x) !is.na(x)))
# avg_hg <- Reduce("+", hg_l) / msr_hg
# sd_hg <- apply(simplify2array(hg_l), 1:2, sd, na.rm = TRUE)



calc_avg_and_na <- function(paths, pc_df, verbose = TRUE) {
  
  amat <- msr_mat <- init_agg_mat(pc_df)
  
  for (file in paths) {
    message(paste(file, Sys.time()))
    cmat <- fread_to_mat(file, genes = pc_df$Symbol)
    non_na <- !is.na(cmat)
    amat[non_na] <- amat[non_na] + cmat[non_na]
    msr_mat[non_na] <- msr_mat[non_na] + 1
    rm(cmat)
    gc(verbose = FALSE)
  }
  
  amat <- amat / msr_mat
  
  return(list(Avg_mat = amat, Msr_mat = msr_mat))
  
}



calc_sd <- function(paths, pc_df, amat, verbose = TRUE) {
  
  var_mat <- msr_mat <- init_agg_mat(pc_df)
  
  for (file in paths) {
    message(paste(file, Sys.time()))
    cmat <- fread_to_mat(file, genes = pc_df$Symbol)
    non_na <- !is.na(cmat)
    diff <- (cmat[non_na] - amat[non_na])^2
    var_mat[non_na] <- var_mat[non_na] + diff
    msr_mat[non_na] <- msr_mat[non_na] + 1
    rm(cmat)
    gc(verbose = FALSE)
  }
  
  sd_mat <- sqrt(var_mat / msr_mat)
  
  return(sd_mat)
}



if (!file.exists(outfile_hg)) {
  
  l_hg <- calc_avg_and_na(paths = cor_paths_hg, pc_df = pc_df_hg)
  sd_hg <- calc_sd(paths = cor_paths_hg, pc_df = pc_df_hg, amat = l_hg$Avg_mat)
  
  saveRDS(list(Avg_mat = l_hg$Avg_mat, SD_mat = sd_hg, Msr_mat = l_hg$Msr_mat),
          file = outfile_hg)
  
}


if (!file.exists(outfile_mm)) {
  
  l_mm <- calc_avg_and_na(paths = cor_paths_mm, pc_df = pc_df_mm)
  sd_mm <- calc_sd(paths = cor_paths_mm, pc_df = pc_df_mm, amat = l_mm$Avg_mat)
  
  saveRDS(list(Avg_mat = l_mm$Avg_mat, SD_mat = sd_mm, Msr_mat = l_mm$Msr_mat),
          file = outfile_mm)
  
}

