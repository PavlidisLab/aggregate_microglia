library(tidyverse)
library(data.table)
library(aggtools)
library(pheatmap)
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



##


l_hg <- readRDS(outfile_hg)
l_hg$CV_mat <- l_hg$SD_mat / l_hg$Avg_mat


l_mm <- readRDS(outfile_mm)
l_mm$CV_mat <- l_mm$SD_mat / l_mm$Avg_mat



summ_df_hg <- cbind(
  mat_to_df(l_hg$Avg_mat, value_name = "Avg"),
  SD = mat_to_df(l_hg$SD_mat, value_name = "SD")[["SD"]],
  CV =  mat_to_df(l_hg$CV_mat, value_name = "CV")[["CV"]],
  N_msr = mat_to_df(l_hg$Msr_mat, value_name = "N_msr")[["N_msr"]]
) %>%
  arrange(desc(Avg))


summ_df_mm <- cbind(
  mat_to_df(l_mm$Avg_mat, value_name = "Avg"),
  SD = mat_to_df(l_mm$SD_mat, value_name = "SD")[["SD"]],
  CV =  mat_to_df(l_mm$CV_mat, value_name = "CV")[["CV"]],
  N_msr = mat_to_df(l_mm$Msr_mat, value_name = "N_msr")[["N_msr"]]
) %>%
  arrange(desc(Avg))




top_avg_hg <- summ_df_hg %>% 
  filter(N_msr > 5) %>% 
  slice_max(abs(Avg), n = 10e3)

top_sd_hg <- summ_df_hg %>% 
  filter(N_msr > 5) %>% 
  slice_max(SD, n = 10e3)

top_cv_hg <- summ_df_hg %>% 
  filter(N_msr > 5) %>% 
  slice_max(CV, n = 10e3)



top_avg_mm <- summ_df_mm %>% 
  filter(N_msr > 5) %>% 
  slice_max(abs(Avg), n = 10e3)

top_sd_mm <- summ_df_mm %>% 
  filter(N_msr > 5) %>% 
  slice_max(SD, n = 10e3)

top_cv_mm <- summ_df_mm %>% 
  filter(N_msr > 5) %>% 
  slice_max(CV, n = 10e3)





cor.test(
  filter(summ_df_hg, N_msr > 5)[["Avg"]],
  filter(summ_df_hg, N_msr > 5)[["SD"]],
  method = "spearman"
)



cor.test(
  filter(summ_df_mm, N_msr > 5)[["Avg"]],
  filter(summ_df_mm, N_msr > 5)[["SD"]],
  method = "spearman"
)



ggplot(top_avg_hg, aes(x = Avg, y = SD)) +
  geom_point(shape = 21, alpha = 0.5) +
  theme_classic() +
  theme(text = element_text(size = 25))



ggplot(top_avg_mm, aes(x = Avg, y = SD)) +
  geom_point(shape = 21, alpha = 0.5) +
  theme_classic() +
  theme(text = element_text(size = 25))
