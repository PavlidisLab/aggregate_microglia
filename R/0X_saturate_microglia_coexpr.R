library(tidyverse)
library(data.table)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")



# Dataset meta and human/mouse data IDs
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)

pc_mm <- read.delim(ref_mm_path)
tfs_mm <- read.delim("/home/amorin/Data/Metadata/AnimalTFDB_mouse_V4.tsv")
agg_mm <- readRDS("/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/aggregate_cormat_FZ_mm.RDS")



path_mm <- "/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/"
pattern <- "_cormat.tsv"
agg_l_mm_path <- "/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/TF_cor_list.RDS"



# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)


# Loads aggregate correlation matrices or NA count matrices for the given dataset
# ids into a list. 
# Pattern: "_NA_mat.tsv" for NA counts

load_mat_to_list  <- function(ids,
                               dir,
                               pattern,  
                               genes,
                               sub_genes = NULL,
                               verbose = TRUE) {
  
  mat_l <- lapply(ids, function(id) {
    if (verbose) message(paste(id, Sys.time()))
    path <- file.path(dir, paste0(id, pattern))
    fread_to_mat(path, genes, sub_genes)
  })
  
  names(mat_l) <- ids
  gc(verbose = FALSE)
  
  return(mat_l)
}



# Call load_agg_mat_list to load/save a local copy

load_or_generate_agg <- function(path, 
                                 ids, 
                                 dir,
                                 pattern,
                                 genes, 
                                 sub_genes = NULL) {
  
  if (!file.exists(path)) {
    
    agg_l <- load_mat_to_list(ids = ids, 
                              dir = dir,
                              genes = genes, 
                              pattern = pattern,
                              sub_genes = sub_genes)
    
    saveRDS(agg_l, path)
  } else {
    agg_l <- readRDS(path)
  }
  
  return(agg_l)
}


# Subset cormat to min thresholded genes, NA->0 and FZ transform

ready_cmat <- function(cmat, keep_genes, keep_tfs) {
  
  cmat <- cmat[keep_mm, keep_tfs_mm]
  cmat[is.na(cmat)] <- 0
  cmat[cmat > 1] <- 1
  cmat[cmat < -1] <- -1
  fisherz(cmat)
  
}



keep_mm <- count_summ$Mouse$Filter_genes
keep_tfs_mm <- intersect(keep_mm, tfs_mm$Symbol)
agg_mat <- agg_mm$Agg_mat[keep_mm, keep_tfs_mm]


cmat_l <- load_mat_to_list(ids = ids_mm[1:2],
                           dir = path_mm,
                           pattern = pattern,
                           genes = pc_mm$Symbol,
                           sub_genes = tfs_mm$Symbol)






cmat_l <- lapply(cmat_l, ready_cmat, keep_mm, keep_tfs_mm)




ind_topk_l <- lapply(cmat_l, function(mat) {
  pair_colwise_topk(mat1 = mat,
                    mat2 = agg_mat,
                    k = 200,
                    ncores = ncore)
})




