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



# TODO:
# TODO: consider if necessary or can be made general with functions.R fread

fread_mat <- function(path) {
  
  mat <- fread(path, sep = "\t")
  rownames <- mat[["V1"]]
  mat <- as.matrix(mat[, -1, drop = FALSE])
  rownames(mat) <- rownames
  
  return(mat)
}




# TODO
# TODO: consider if load logic can be used with fread_to_mat()
# TODO: n_iter is pretty hacky, consider replacing when all is ran

average_and_save_each_network <- function(ids, 
                                          grn_dir, 
                                          ncore, 
                                          # force_resave = FALSE) {
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

