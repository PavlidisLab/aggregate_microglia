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




# Loads aggregate correlation matrices or NA count matrices for the given dataset
# ids into a list. 
# Pattern: "_NA_mat.tsv" for NA counts

load_agg_mat_list  <- function(ids,
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
    
    agg_l <- load_agg_mat_list(ids = ids, 
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



tt1 <- load_agg_mat_list(ids = ids_mm[1],
                         dir = path_mm,
                         pattern = pattern,
                         genes = pc_mm$Symbol,
                         sub_genes = tfs_mm$Symbol)









mat1 <- tt1$HypoMap
mat2 <- agg_mm$Agg_mat[, tfs_mm$Symbol]

mat1[is.na(mat1)] <- 0
mat2[is.na(mat2)] <- 0


topk <- pair_colwise_topk(mat1 = mat1,
                          mat2 = mat2,
                          k = 200,
                          ncores = ncore)


check_gene <- "Rpl11"

df <- data.frame(Symbol = pc_mm$Symbol,
                 Agg = agg_mm$Agg_mat[, check_gene], 
                 Msr = length(ids_mm) - agg_mm$NA_mat[, check_gene],
                 Is_TF = pc_mm$Symbol %in% tfs_mm$Symbol)

df <- filter(df, Msr >= 5)

hist(df$Agg, breaks = 1000 )


keep_genes <- names(diag(agg_mm$NA_mat) >= 5)
keep_tfs <- intersect(tfs_mm$Symbol, keep_genes)

ph <- pheatmap::pheatmap(mat2[keep_genes, keep_tfs],
                         cluster_rows = TRUE,
                         cluster_cols = TRUE)
