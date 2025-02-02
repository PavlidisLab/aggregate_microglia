## Cell by Gene single
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(Seurat)
library(CSCORE)
source("R/00_config.R")
source("R/utils/functions.R")


args <- commandArgs(trailingOnly = TRUE)
id <- args[1]

mcg_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)
ct <- mcg_meta %>% filter(ID == id) %>% pull(Cell_type)
species <- mcg_meta %>%  filter(ID == id) %>% pull(Species)

dat_dir <- file.path(sc_dir, id)
dat_path <- list.files(dat_dir, pattern = "^.*_seurat.*.RDS", full.names = TRUE)
outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))


ens <- if (species == "Human") {
  read.delim(ens_hg_path, stringsAsFactors = FALSE)
} else if (species == "Mouse") {
  read.delim(ens_mm_path, stringsAsFactors = FALSE)
} else {
  stop("Species not recognized")
}


if (!file.exists(outfile)) {
  
  message(paste(id, Sys.time()))
  
  # Load Seurat object -- if multiple, merge into one
  dat <- lapply(dat_path, readRDS)
  dat <- reduce(dat, merge)
  
  # Load existing processed mcg data to ensure exact match in cells
  mcg_dat <- readRDS(mcg_dat_path)
  keep_ids <- mcg_dat[[id]]$Meta$ID
  dat <- dat[, keep_ids]
  
  # Convert Ensembl IDs to symbols... need to extract count matrix and rebuild
  counts <- suppressWarnings(dat@assays$RNA$counts)

  if (length(counts) == 0 || sum(counts == 0)) {
    counts <- dat@assays$RNA$data
    stopifnot(is.integer(counts@i))
  }
  
  keep_ens <- intersect(rownames(counts), ens$Gene_ID)
  counts <- counts[keep_ens, ]
  rownames(counts) <- ens$Symbol[match(keep_ens, ens$Gene_ID)]
  
  # Keeping genes detected in at least 20 cells. Arb. threshold for min genes
  keep_genes <- names(which(rowSums(counts > 0) >= 20))

  if (length(keep_genes) <= 100) {
    stop("Too few genes")
  }
  
  counts <- counts[keep_genes, ]

  # Remake Seurat object with refseq named count matrix
  dat <- CreateSeuratObject(counts = counts, meta.data = dat@meta.data)
  
  # Set raw counts
  DefaultAssay(dat) <- "RNA"
  stopifnot(is.integer(dat[["RNA"]]$counts@i))

  # Run CSCORE and save output
  res <- CSCORE(dat, genes = keep_genes) 
  saveRDS(res, outfile)
  
} 
