## Code to load and format the TF-target-importance lists produced by GRNBoost2

library(parallel)
library(data.table)


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



# Generate a network/importance score matrix for each output network list

ready_networks <- function(paths, ncore) {
  
  mat_l <- mclapply(paths, function(x) {
    network_list <- read_network_list(x)
    network_mat <- network_list_to_mat(network_list)
  }, mc.cores = ncore)
  
  return(mat_l)
}



# Element-wise average of all matrices in mat_l. Assumes identical dimensions.

average_mat_list <- function(mat_l)  Reduce("+", mat_l) / length(mat_l)

