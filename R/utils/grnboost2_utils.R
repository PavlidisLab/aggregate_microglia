## Code to load and format the TF-target-importance lists produced by GRNBoost2

library(parallel)
library(data.table)
source("R/utils/functions.R")


# GRNBoost2 output saved as table of [TF, Gene, Importance_score]

read_network_list <- function(path) {
  fread(path, sep = "\t", header = FALSE, col.names = c("TF", "Gene", "Importance"))
}



# GRNBoost2 gives a list of TF-target-importance scores. Convert to a gene x TF
# matrix of importance scores

network_list_to_mat <- function(network_list) {
  
  tfs <- unique(network_list$TF)
  genes <- unique(network_list$Gene)
  
  imp_mat <- matrix(0, nrow = length(genes), ncol = length(tfs))
  rownames(imp_mat) <- genes
  colnames(imp_mat) <- tfs
  
  tf_ix <- match(network_list$TF, tfs)
  gene_ix <- match(network_list$Gene, genes)
  
  for (i in 1:nrow(network_list)) {
    imp_mat[gene_ix[i], tf_ix[i]] <- network_list$Importance[i]
  }
  
  return(imp_mat)
}



# Load each network into a list of matrices and consistently order names

ready_networks <- function(paths, ncore) {
  
  mat_l <- mclapply(paths, function(x) {
    network_list <- read_network_list(x)
    network_mat <- network_list_to_mat(network_list)
  }, mc.cores = ncore)
  
  tfs <- Reduce(union, lapply(mat_l, colnames))
  genes <- Reduce(union, lapply(mat_l, rownames))
  mat_l <- lapply(mat_l, function(mat) mat[genes, tfs])
  
  return(mat_l)
}



# Element-wise average of all matrices in mat_l. Assumes identical dimensions.

average_mat_list <- function(mat_l)  Reduce("+", mat_l) / length(mat_l)



# Use fread:: to load .tsv into matrix, making the first column into rownames
# TODO: consider if necessary or can be made general with functions.R fread

fread_mat <- function(path) {
  
  mat <- fread(path, sep = "\t")
  rownames <- mat[["V1"]]
  mat <- as.matrix(mat[, -1, drop = FALSE])
  rownames(mat) <- rownames
  
  return(mat)
}




# For each dataset ID, average the GRNBoost2 iterations into one matrix
# TODO: consider if load logic can be used with fread_to_mat()
# TODO: n_iter is pretty hacky, consider replacing when all is ran

average_and_save_each_network <- function(ids, 
                                          grn_dir, 
                                          ncore, 
                                          n_iter = 100) {
  
  avg_grn_l <- lapply(ids, function(id) {
    
    id_dir <- file.path(grn_dir, id)
    iter_paths <- list.files(id_dir, full.names = TRUE, pattern = ".*_iter.*")
    avg_grn_path <- file.path(id_dir, paste0(id, "_average_GRN.tsv"))
    
    if (!file.exists(avg_grn_path) || length(iter_paths) < n_iter) {

      mat_l <- ready_networks(iter_paths, ncore)
      avg_grn <- average_mat_list(mat_l)
      fwrite_mat(avg_grn, avg_grn_path)
      avg_grn
      
    } else {
      
      avg_grn <- fread_mat(avg_grn_path)
      
    }
  })
  names(avg_grn_l) <- ids
  
  return(avg_grn_l)
}



# Compile all averaged GRNBoost2 networks into one global average matrix

average_all_networks <- function(avg_grn_l, keep_genes, keep_tfs) {
  
  # Init tracking matrix that includes final genes to keep
  track_mat <- matrix(0, nrow = length(keep_genes), ncol = length(keep_tfs))
  rownames(track_mat) <- keep_genes
  colnames(track_mat) <- keep_tfs
  
  # Make each matrix have equal dimensions/ordering
  mat_l <- lapply(avg_grn_l, function(x) {
    
    common_gene <- intersect(keep_genes, rownames(x))
    common_tf <- intersect(keep_tfs, colnames(x))
    track_mat[common_gene, common_tf] <- x[common_gene, common_tf]
    track_mat
    
  })
  
  avg_mat <- average_mat_list(mat_l)
  
  return(avg_mat)
}



# All rank the averaged importance scores

rank_aggregate_grn <- function(avg_grn_l, keep_genes, keep_tfs) {
  
  # Init aggregate matrix and tracking matrix to hold averaged importance scores
  amat <- matrix(0, nrow = length(keep_genes), ncol = length(keep_tfs))
  rownames(amat) <- keep_genes
  colnames(amat) <- keep_tfs
  keep_mat <- amat
  
  for (mat in avg_grn_l) {
    
    # Fill tracking matrix with measured values
    common_genes <- intersect(rownames(mat), keep_genes)
    common_tfs <- intersect(colnames(mat), keep_tfs)
    keep_mat[common_genes, common_tfs] <- mat[common_genes, common_tfs]
    
    # Rank standardize
    rmat <- allrank_mat(-keep_mat, ties_arg = "min")
    rmat <- rmat / max(rmat, na.rm = TRUE)
    amat <- amat + rmat
    
  }
  
  # Finalize aggregate matrix
  amat <- allrank_mat(-amat, ties_arg = "min")
  amat <- amat / max(amat, na.rm = TRUE)

  return(amat)
}
