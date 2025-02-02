## Process count matrix and get aggregate correlation for Velmeshev 2019
## -----------------------------------------------------------------------------

library(CSCORE)
library(tidyverse)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")

id <- "Velmeshev"
sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/jules_garreau_sc_datasets/", id)
dat_path <- file.path(sc_dir, "matrix.mtx")
meta_path <- file.path(sc_dir, "meta.txt")
genes_path <- file.path(sc_dir, "genes.tsv")
outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))
mcg_dat_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)
ct <- mcg_dat_meta %>% filter(ID == id) %>% pull(Cell_type)

pc <- read.delim(ref_hg_path, stringsAsFactors = FALSE)



if (!file.exists(outfile)) {
  
  # Count matrix, metadata, and genes need to be loaded individually
  
  mat <- Matrix::readMM(dat_path)
  meta <- read.delim(meta_path, stringsAsFactors = FALSE)
  genes <- read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)
  
  colnames(mat) <- meta$cell
  rownames(mat) <- genes$V2
  
  
  # Ready metadata
  
  change_colnames <- c(Cell_type = "cluster", ID = "cell")
  
  meta <- meta %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(assay = "10x 3' v1") %>% 
    add_count_info(mat = mat)
  
  # Remove cells failing QC and keep only protein coding genes
  mat <- rm_low_qc_cells(mat, meta) %>% get_pcoding_only(pcoding_df = pc)
  meta <- filter(meta, ID %in% colnames(mat))

  # Subset to microglia cells
  mcg_meta <- meta %>% filter(Cell_type %in% ct)
  mcg_ids <- mcg_meta$ID
  mcg_mat <- mat[, mcg_ids]
  
  # Keeping genes detected in at least 20 cells. Arb. threshold for min genes
  keep_genes <- names(which(rowSums(mcg_mat > 0) >= 20))

  if (length(keep_genes) <= 100) {
    stop("Too few genes")
  }
  mcg_mat <- mcg_mat[keep_genes, ]

  # Make Seurat object for input to CSCORE and ensure counts are integer
  mcg_dat <- CreateSeuratObject(counts = mcg_mat, meta.data = mcg_meta)
  DefaultAssay(mcg_dat) <- "RNA"
  stopifnot(is.integer(mcg_dat[["RNA"]]$counts@i)) 
    
  # Run CSCORE and save output
  res <- CSCORE(mcg_dat, genes = keep_genes) 
  saveRDS(res, outfile)

}
  
