## Code to load and format the TF-target-importance lists produced by GRNBoost2

library(parallel)



# GRNBoost2 gives a list of TF-target-importance scores. Convert to a gene x TF
# matrix of importance scores

network_list_to_mat <- function(network_list, tfs, genes, ncore) {
  
  imp_vec <- setNames(rep(0, length(genes)), genes) # init importance vector
  
  imp_l <- mclapply(tfs, function(tf) { # extract gene scores for each TF
    
    tf_df <- filter(network_list, V1 == tf)
    gene_match <- match(tf_df$V2, genes)
    imp_vec[gene_match] <- tf_df$V3
    imp_vec
    
  }, mc.cores = ncore)
  
  # Bind importance vectors into a matrix
  imp_mat <- do.call(cbind, imp_l)
  colnames(imp_mat) <- tfs
  
  return(imp_mat)
}



# Generate a network/importance score matrix for each output network list

ready_networks <- function(paths, tfs, genes, ncore) {
  
  mat_l <- lapply(paths, function(x) {
    
    network_list <- read.delim(x, header = FALSE)
    if (!all(network_list$V1 %in% tfs)) stop(paste(x, "doesn't have all TFs"))
    network_mat <- network_list_to_mat(network_list, tfs, genes, ncore)
    
  })
  
  return(mat_l)
}



# Element-wise average of all matrices in mat_l. Assumes identical dimensions.

average_mat_list <- function(mat_l)  Reduce("+", mat_l) / length(mat_l)

