## Process count matrix and get aggregate correlation for ROSMAP
## -----------------------------------------------------------------------------

library(CSCORE)
library(tidyverse)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")

id <- "ROSMAP"
sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/jules_garreau_sc_datasets/", id)
dat_path <- file.path(sc_dir, "sc_exprmats_rosmap.rds")
outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))
mcg_dat_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)
ct <- mcg_dat_meta %>% filter(ID == id) %>% pull(Cell_type)

pc <- read.delim(ens_hg_path, stringsAsFactors = FALSE)



if (!file.exists(outfile)) {
  
  dat <- readRDS(dat_path)
  
  mat <- dat$raw$ctmat
  
  
  # Ready metadata
  
  change_colnames <- c(Cell_type = "cell_type", ID = "sample")
  
  meta <- dat$raw$samples %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(assay = "10x 3' v2") %>% 
    add_count_info(mat = mat)
  
  # Remove cells failing QC and keep only protein coding genes
  mat <- rm_low_qc_cells(mat, meta) %>% 
    ensembl_to_symbol(ensembl_df = pc) %>% 
    get_pcoding_only(pcoding_df = pc)
  
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
  
