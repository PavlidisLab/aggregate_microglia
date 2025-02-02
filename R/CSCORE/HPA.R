## HPA
## -----------------------------------------------------------------------------

library(CSCORE)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")

id <- "HPA"
species <- "Human"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)
mcg_dat_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)

ct <- mcg_dat_meta %>% filter(ID == id) %>% pull(Cell_type)

outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))


pc <- read.delim(ens_hg_path, stringsAsFactors = FALSE)


# Files were directly downloaded from human protein atlas, see HPA_download.sh in dat_dir
tissue_dir <- list.dirs(dat_dir, full.names = TRUE, recursive = FALSE)
cluster_meta_path <- file.path(dat_dir, "rna_single_cell_cluster_description.tsv")


# Load each tissue's count matrix and cell IDs into a list, then match the 
# cell ID/cluster numbers with the respective cell type in cluster_meta.

load_and_match_data <- function(cluster_meta, tissue_dir) {
  
  
  tissues <- unique(cluster_meta$Tissue)
  
  
  dat_l <- lapply(1:length(tissues), function(x) {
    
    # Load mat and meta
    
    tissue <- tissues[x]
    cell_ids <- read.delim(file.path(tissue_dir[x], "cell_data.tsv"))
    mat <- fread(file.path(tissue_dir[x], "read_count.tsv"), header = TRUE)
    
    # Match cluster number to cell type, cleaning up names and setting unique ID
    
    ct <- cluster_meta %>% 
      filter(Tissue == tissue) %>% 
      mutate(
        Cell_type = paste0(tissue, "_", Cell_type),
        Cell_type = str_replace_all(Cell_type, " ", "_"))
    
    cell_meta <- cell_ids %>%
      mutate(cluster = paste0("c-", cluster)) %>%
      dplyr::rename(Cluster = cluster) %>%
      left_join(.,
                ct[, c("Tissue", "Cluster", "Cell_type", "Cell_type_group")],
                by = "Cluster") %>% 
      dplyr::rename(ID = cell_id) %>% 
      mutate(ID = paste0(ID, "_", Cell_type))
    
    # Format matrix
    
    genes <- mat[[1]]
    mat <- as.matrix(mat[, -1])
    rownames(mat) <- genes
    colnames(mat) <- cell_meta$ID
    
    gc(verbose = FALSE)
    list(Mat = mat, Meta = cell_meta)
    
  })
  
  names(dat_l) <- tissues
  
  return(dat_l)
  
}



if (!file.exists(outfile)) {
  
  # For HPA must load each tissue individually before binding into one mat
  
  cluster_meta <- read.delim(cluster_meta_path)
  colnames(cluster_meta) <- str_replace_all(colnames(cluster_meta), "\\.", "_")
  
  
  # Tissue names slightly differ between paths and meta, inspect to ensure match
  
  tissues <- data.frame(
    Dir = list.dirs(dat_dir, full.names = FALSE, recursive = FALSE),
    Meta = unique(cluster_meta$Tissue)
  )
  
  stopifnot(identical(
    str_to_lower(str_replace(tissues$Dir, "_", "")),
    str_to_lower(str_replace(tissues$Meta, " ", ""))
  ))

  
  dat <- load_and_match_data(cluster_meta, tissue_dir)
  
  
  mat <- do.call(cbind, lapply(dat, `[[`, "Mat"))
  meta <- do.call(rbind, lapply(dat, `[[`, "Meta"))

  
  # Ready metadata
  
  meta <- meta %>% 
    mutate(assay = "NA") %>% 
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
  
