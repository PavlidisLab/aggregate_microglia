## This script downloads, preprocesses, and generates the aggregate correlation
## network for the Tabula Muris dataset
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6642641/
## https://figshare.com/articles/dataset/Robject_files_for_tissues_processed_by_Seurat/5821263
## -----------------------------------------------------------------------------

library(CSCORE)
library(tidyverse)
library(Seurat)
library(data.table)
source("R/00_config.R")
source("R/utils/functions.R")

id <- "TabulaMuris"
species <- "Mouse"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)
mcg_dat_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)

ct <- mcg_dat_meta %>% filter(ID == id) %>% pull(Cell_type)

outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))

dl_table <- read.delim(file.path(dat_dir, "TabulaMuris_download.tsv"), stringsAsFactors = FALSE, header = FALSE)


# Download each tissue as an Robj file

for (i in 1:nrow(dl_table)) {
  
  dl_path <- file.path(dat_dir, dl_table[i, 1])
  dl_url <- dl_table[i, 2]
  
  if (!file.exists(dl_path)) {
    curl::curl_download(url = dl_url, destfile = dl_path, quiet = FALSE)
  }
}


dat_path <- list.files(dat_dir, pattern = ".Robj", full.names = TRUE)

pc <- read.delim(ref_mm_path, stringsAsFactors = FALSE)



if (!file.exists(outfile)) {
  
  # Load and update Seurat objects, each of which is named "tiss", then merge
  
  dat_l <- lapply(dat_path, function(x) {
    load(x)
    UpdateSeuratObject(tiss)
  })
  
  dat <- reduce(dat_l, merge)
  gc(verbose = FALSE)

  # Extract count matrix: default counts slot, but use data slot if counts empty
  
  mat <- GetAssayData(dat, slot = "counts")
  
  if (length(mat) == 0 || all(rowSums(mat) == 0)) {
    mat <- GetAssayData(dat, slot = "data")
  }
  
  
  # Ready metadata
  # "TabulaMuris" remove NA cell types
  
  meta <- dat[[]] %>% 
    dplyr::rename(Cell_type = cell_ontology_class) %>% 
    rownames_to_column(var = "ID") %>% 
    filter(!is.na(Cell_type))
  
  mat <- mat[, meta$ID]
  meta <- add_count_info(mat, meta)
  
  meta <- mutate(meta, assay = "Smart-seq")
  
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
  
