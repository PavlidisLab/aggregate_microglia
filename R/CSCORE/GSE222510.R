## GSE222510
## -----------------------------------------------------------------------------

library(CSCORE)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")

id <- "GSE222510"
species <- "Mouse"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)
mcg_dat_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)

ct <- mcg_dat_meta %>% filter(ID == id) %>% pull(Cell_type)

outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))


pc <- read.delim(ref_mm_path, stringsAsFactors = FALSE)


# Files were directly downloaded from GEO, see GSE222510_download.sh in dat_dir
meta_path <- file.path(dat_dir, paste0(id, "_metadata.txt"))
mat_path <- file.path(dat_dir, paste0(id, "_counts.mtx"))
features_path <- file.path(dat_dir, paste0(id, "_features.tsv"))
barcodes_path <- file.path(dat_dir, paste0(id, "_barcodes.tsv"))



if (!file.exists(outfile)) {
  
  # Load metadata and the count matrix
  
  meta <- read.delim(meta_path)
  mat <- Matrix::readMM(mat_path)
  features <- read.delim(features_path, header = FALSE)
  barcodes <- read.delim(barcodes_path, header = FALSE)
  
  rownames(mat) <- features$V1
  colnames(mat) <- barcodes$V1
  
  stopifnot(identical(colnames(mat), meta$barcode))

  
  # Ready metadata
  
  change_colnames <- c(Cell_type = "cell_type", ID = "barcode")
  
  meta <- meta %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(assay = "10x 3' v2") %>% 
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
  
