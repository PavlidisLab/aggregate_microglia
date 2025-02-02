## GSE199317
## -----------------------------------------------------------------------------

library(CSCORE)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")

id <- "GSE199317"
species <- "Mouse"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)
mcg_dat_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)

ct <- mcg_dat_meta %>% filter(ID == id) %>% pull(Cell_type)

outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))


pc <- read.delim(ref_mm_path, stringsAsFactors = FALSE)


# Files were directly downloaded from GEO, see GSE199317_download.sh in dat_dir
meta_path <- list.files(dat_dir, pattern = "metadata", full.names = TRUE)
mat_path <- list.files(dat_dir, pattern = "mtx", full.names = TRUE)
features_path <- list.files(dat_dir, pattern = "features", full.names = TRUE)
barcodes_path <- list.files(dat_dir, pattern = "barcodes", full.names = TRUE)



if (!file.exists(outfile)) {
  
  # Load metadata and the count matrix
  
  meta_l <- lapply(meta_path, read.delim)
  
  meta <- rbind(
    meta_l[[1]], 
    dplyr::rename(meta_l[[2]], cell_type = "cell_type_major") %>% select(cell, cell_type),
    meta_l[[3]]
  )
  
  features_l <- lapply(features_path, read.delim, header = FALSE)
  barcodes_l <- lapply(barcodes_path, read.delim, header = FALSE)
  
  mat_l <- lapply(1:length(mat_path), function(x) {
    mat <- Matrix::readMM(mat_path[[x]])
    rownames(mat) <- features_l[[x]]$V1
    colnames(mat) <- barcodes_l[[x]]$V1
    return(mat)
  })
  
  mat <- do.call(cbind, mat_l)
  common <- intersect(colnames(mat), meta$cell)
  mat <- mat[, common]
  meta <- filter(meta, cell %in% common)
  
  stopifnot(identical(colnames(mat), meta$cell))
  
  
  # Ready metadata
  # "GSE199317" clean up cell types and remove ambiguous
  
  change_colnames <- c(Cell_type = "cell_type", ID = "cell")
  
  meta <- meta %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(
      assay = "10x 3' v2",
      Cell_type = str_replace(Cell_type, "Neutrophils", "Neutrophil")
      ) %>% 
    add_count_info(mat = mat)
  
  
  meta <- filter(meta, !str_detect(Cell_type, "[:digit:]+|low UMI|doublet")) 
  mat <- mat[, meta$ID]
  
  
  # Remove cells failing QC and keep only protein coding genes
  mat <- rm_low_qc_cells(mat, meta) %>% get_pcoding_only(pcoding_df = pc)
  meta <- filter(meta, ID %in% colnames(mat))

  # Subset to microglia cells
  mcg_meta <- meta %>% filter(Cell_type %in% unlist(str_split(ct, ";")))
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
  
