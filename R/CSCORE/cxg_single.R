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
dat_path <- file.path(dat_dir, paste0(id, "_cellxgene_seurat.RDS"))
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
  
  # Load Seurat object and subset to microglia
  dat <- readRDS(dat_path)
  dat <- dat[, dat$cell_type == ct]
  
  # Ensembl -> refseq
  counts <- dat@assays$RNA$counts
  keep_ens <- intersect(rownames(counts), ens$Gene_ID)
  counts <- counts[keep_ens, ]
  rownames(counts) <- ens$Symbol[match(keep_ens, ens$Gene_ID)]
  
  # Remake Seurat object with refseq named count matrix
  dat <- CreateSeuratObject(counts = counts, meta.data = dat@meta.data)
  
  # Set raw counts
  DefaultAssay(dat) <- "RNA"
  stopifnot(is.integer(dat[["RNA"]]$counts@i))

  # Run CSCORE and save output
  res <- CSCORE(dat, genes = keep_genes) 
  saveRDS(res, outfile)
  
} 