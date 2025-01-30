## GSE180928
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(Seurat)
library(CSCORE)
source("R/00_config.R")
source("R/utils/functions.R")

id <- "GSE180928"
mcg_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)
ct <- mcg_meta %>% filter(ID == id) %>% pull(Cell_type)

# Input count/meta and output CSCORE paths
sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata", id)
dat_path <- file.path(sc_dir, paste0(id, "_filtered_cell_counts.csv"))
meta_path <- file.path(sc_dir, paste0(id, "_metadata.csv"))
outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))

# Protein coding genes
pc <- read.delim(ref_hg_path, stringsAsFactors = FALSE)

# Get microglia raw count matrix and run CSCORE
if (!file.exists(outfile)) {
  
  # Read meta and count matrix and match IDs
  count_meta <- read.delim(meta_path, sep = ",")
  mat <- read_count_mat(dat_path)
  colnames(mat) <- str_replace_all(colnames(mat), "\\.", "-")
  mat <- mat[, count_meta$X]

  stopifnot(identical(colnames(mat), count_meta$X))

  # Ready metadata
  change_colnames <- c(Cell_type = "Lineage", ID = "X")
  
  count_meta <- count_meta %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(assay = "10x 3' v2/v3") %>% 
    add_count_info(mat = mat)

  # Remove cells failing QC and keep only protein coding genes
  mat <- rm_low_qc_cells(mat, count_meta) %>% get_pcoding_only(pcoding_df = pc)
  count_meta <- filter(count_meta, ID %in% colnames(mat))

  # Subset to microglia cells
  mcg_count_meta <- count_meta %>% filter(Cell_type %in% ct)
  mcg_ids <- mcg_count_meta$ID
  mcg_mat <- mat[, mcg_ids]
  
  # Keeping genes detected in at least 20 cells. Arb. threshold for min genes
  keep_genes <- names(which(rowSums(mcg_mat > 0) >= 20))

  if (length(keep_genes) <= 100) {
    stop("Too few genes")
  }
  mcg_mat <- mcg_mat[keep_genes, ]

  # Make Seurat object for input to CSCORE and ensure counts are integer
  mcg_dat <- CreateSeuratObject(counts = mcg_mat, meta.data = mcg_count_meta)
  DefaultAssay(mcg_dat) <- "RNA"
  stopifnot(is.integer(mcg_dat[["RNA"]]$counts@i)) 
    
  # Run CSCORE and save output
  res <- CSCORE(mcg_dat, genes = keep_genes) 
  saveRDS(res, outfile)

}