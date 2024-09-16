# Testing implementations of sparse matrix correlation 
# https://github.com/koheiw/proxyC
# https://github.com/cysouw/qlcMatrix
# https://slowkow.com/notes/sparse-matrix/
# https://saket-choudhary.me/blog/2022/03/09/sparsespearman/
# https://github.com/saketkc/blog/blob/main/2022-03-10/SparseSpearmanCorrelation2.ipynb
# ------------------------------------------------------------------------------

# install.packages("proxyC")
# install.packages("qlcMatrix")

library(parallel)
library(qlcMatrix)
library(proxyC)
library(tidyverse)
library(data.table)
source("R/00_config.R")
source("R/utils/dev_functions.R")
source("R/utils/plot_functions.R")


# Testing with a real dataset
id <- "GSE180928"
ct <- "Microglia"
out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
dat <- readRDS(processed_path)
meta <- dat$Meta
mat <- dat$Mat

# For sampling genes
set.seed(56)
n_samp <- 1000

# Get sparse cell type count matrix: this auto filters low count genes to 0.
mat_sparse <- subset_celltype_and_filter(mat, meta, cell_type = ct)

# For speed subset matrix to n_samp genes that are measured as well as all zero
which_nonzero <- colSums(mat_sparse > 0) > 0
nonzero <- names(which_nonzero[which_nonzero])
allzero <- names(which_nonzero[!which_nonzero])
sample_nonzero <- sample(nonzero, n_samp)
sample_allzero <- sample(allzero, n_samp)
mat_sparse_sub <- mat_sparse[, c(sample_nonzero, sample_allzero)]

# Keep only genes that are not all zero
mat_sparse_keep <- mat_sparse[, nonzero]

# Dense matrices
mat_dense <- as.matrix(mat_sparse)
mat_dense_sub <- as.matrix(mat_sparse_sub)
mat_dense_keep <- as.matrix(mat_sparse_keep)


# Functions
# ------------------------------------------------------------------------------


# Original implementation of sparse ranks followed by sparse cor
# https://saket-choudhary.me/blog/2022/03/09/sparsespearman/

SparsifiedRanks <- function(X) {
  if (class(X)[1] != "dgCMatrix") {
    X <- as(object = X, Class = "dgCMatrix")
  }
  non_zeros_per_col <- diff(x = X@p)
  n_zeros_per_col <- nrow(x = X) - non_zeros_per_col
  offsets <- (n_zeros_per_col - 1) / 2
  x <- X@x
  ## split entries to columns
  col_lst <- split(x = x, f = rep.int(1:ncol(X), non_zeros_per_col))
  ## calculate sparsified ranks and do shifting
  sparsified_ranks <- unlist(x = lapply(X = seq_along(col_lst), 
                                        FUN = function(i) rank(x = col_lst[[i]]) + offsets[i]))
  ## Create template rank matrix
  X.ranks <- X
  X.ranks@x <- sparsified_ranks
  return(X.ranks)
}


SparseSpearmanCor <- function(X, Y = NULL, cov = FALSE) {
  
  # Get sparsified ranks
  rankX <- SparsifiedRanks(X)
  if (is.null(Y)){
    # Calculate pearson correlation on rank matrices
    return (corSparse(X=rankX, cov=cov))
  }
  rankY <- SparsifiedRanks2(Y)
  return(corSparse( X = rankX, Y = rankY, cov = cov))
}



# Re-impl: sparse ranks just adds parallel, identical otherwise

sparse_rank <- function(mat, ncores = 1) {
  
  stopifnot(inherits(mat, "dgCMatrix"))
  
  nonzero_per_col <- diff(mat@p)
  zero_per_col <- nrow(mat) - nonzero_per_col
  offsets <- (zero_per_col - 1) / 2
  x <- mat@x
  
  ## split entries to columns
  col_l <- split(x, f = rep.int(1:ncol(mat), nonzero_per_col))
  
  # Rank and offset
  sparse_ranks <- mclapply(seq_along(col_l), function(i) {
    rank(col_l[[i]]) + offsets[i]
  }, mc.cores = ncores)
  
  mat@x <- unlist(sparse_ranks)
  
  return(mat)
}



# Sparse pcor wrappers for proxyC or qlcMatrix, +/- filter matrix of all nonzero

sparse_pcor1 <- function(mat) {
  cmat <- qlcMatrix::corSparse(mat)
  colnames(cmat) <- rownames(cmat) <- colnames(mat)
  return(cmat)
}


sparse_pcor2 <- function(mat) {
  cmat <- proxyC::simil(mat, margin = 2, method = "correlation", use_nan = TRUE)
  return(as.matrix(cmat))
}


sparse_pcor3 <- function(mat) {
  
  cmat_all <- matrix(NA, ncol(mat), ncol(mat))
  colnames(cmat_all) <- rownames(cmat_all) <- colnames(mat)
  which_nonzero <- colSums(mat > 0) > 0
  mat <- mat[, which_nonzero]
  cmat <- qlcMatrix::corSparse(mat)
  cmat_all[which_nonzero, which_nonzero] <- cmat
  
  return(cmat_all)
}


sparse_pcor4 <- function(mat) {
  
  cmat_all <- matrix(NA, ncol(mat), ncol(mat))
  colnames(cmat_all) <- rownames(cmat_all) <- colnames(mat)
  which_nonzero <- colSums(mat > 0) > 0
  mat <- mat[, which_nonzero]
  cmat <- proxyC::simil(mat, margin = 2, method = "correlation", use_nan = TRUE)
  cmat_all[which_nonzero, which_nonzero] <- as.matrix(cmat)
  
  return(cmat_all)
}


# Scor calculated over nonzero columns, as zero columns affect the offset logic.
# Fill out all zero columns with NAs. Inspecting both proxyC and qlcMatrix 

sparse_scor1 <- function(mat, ncores = 1) {
  
  cmat_all <- matrix(NA, ncol(mat), ncol(mat))
  colnames(cmat_all) <- rownames(cmat_all) <- colnames(mat)
  which_nonzero <- colSums(mat > 0) > 0
  mat <- mat[, which_nonzero]
  
  rmat <- sparse_rank(mat, ncores)
  cmat <- qlcMatrix::corSparse(rmat)
  colnames(cmat) <- rownames(cmat) <- colnames(mat)
  cmat_all[which_nonzero, which_nonzero] <- cmat
  cmat_all <- diag_to_one(cmat_all)
  
  return(cmat_all)
}


sparse_scor2 <- function(mat, ncores = 1) {
  
  cmat_all <- matrix(NA, ncol(mat), ncol(mat))
  colnames(cmat_all) <- rownames(cmat_all) <- colnames(mat)
  which_nonzero <- colSums(mat > 0) > 0
  mat <- mat[, which_nonzero]
  
  rmat <- sparse_rank(mat, ncores)
  cmat <- proxyC::simil(rmat, margin = 2, method = "correlation", use_nan = TRUE)
  cmat <- as.matrix(cmat)
  cmat_all[which_nonzero, which_nonzero] <- cmat
  cmat_all <- diag_to_one(cmat_all)
  
  return(cmat_all)
}




# Run and check
# ------------------------------------------------------------------------------


# Original sparse implementation
scor_sparse <- SparseSpearmanCor(mat_sparse)
colnames(scor_sparse) <- rownames(scor_sparse) <- colnames(mat_sparse)
scor_sparse <- diag_to_one(scor_sparse)

scor_sparse_sub <- SparseSpearmanCor(mat_sparse_sub)
colnames(scor_sparse_sub) <- rownames(scor_sparse_sub) <- colnames(mat_sparse_sub)
scor_sparse_sub <- diag_to_one(scor_sparse_sub)

scor_sparse_keep <- SparseSpearmanCor(mat_sparse_keep)
colnames(scor_sparse_keep) <- rownames(scor_sparse_keep) <- colnames(mat_sparse_keep)
scor_sparse_keep <- diag_to_one(scor_sparse_keep)


# Using updated sparse that ignores nonzero columns: qlcMatrix
scor_sparse1 <- sparse_scor1(mat_sparse, ncores = ncore)
scor_sparse_sub1 <- sparse_scor1(mat_sparse_sub, ncores = ncore)
scor_sparse_keep1 <- sparse_scor1(mat_sparse_keep, ncores = ncore)


# Using updated sparse that ignores nonzero columns: proxyC
scor_sparse2 <- sparse_scor2(mat_sparse, ncores = ncore)
scor_sparse_sub2 <- sparse_scor2(mat_sparse_sub, ncores = ncore)
scor_sparse_keep2 <- sparse_scor2(mat_sparse_keep, ncores = ncore)


# Dense spearman ("gold standard")
scor_dense <- WGCNA::cor(mat_dense, method = "spearman", use = "everything", nThreads = ncores)
scor_dense_sub <- WGCNA::cor(mat_dense_sub, method = "spearman", use = "everything", nThreads = ncores)
scor_dense_keep <- WGCNA::cor(mat_dense_keep, method = "spearman", use = "everything", nThreads = ncores)




check_gene <- "SRGAP2B"




# Dense and original sparse sub are giving identical answers
all.equal(scor_dense_sub, scor_sparse_sub, tolerance = 1e-6)

# As are re-impls
all.equal(scor_dense_sub, scor_sparse_sub1, tolerance = 1e-6)
all.equal(scor_dense_sub, scor_sparse_sub2, tolerance = 1e-6)


# Dense and original sparse keeping only nonzero genes gives identical answers
all.equal(scor_dense_keep, scor_sparse_keep, tolerance = 1e-6)

# Ditto for re-impls
all.equal(scor_dense_keep, scor_sparse_keep1, tolerance = 1e-6)
all.equal(scor_dense_keep, scor_sparse_keep2, tolerance = 1e-6)


# Dense and original sparse with all genes does not give identical
plot(scor_dense[, check_gene], scor_sparse[, check_gene])
all.equal(scor_dense, scor_sparse, tolerance = 1e-6)

# Re-impls however to ignore nonzero genes gives identical
all.equal(scor_dense, scor_sparse1, tolerance = 1e-6)
all.equal(scor_dense, scor_sparse2, tolerance = 1e-6)






# testing sparse pcor time


a1 <- sparse_pcor1(mat_sparse_sub)
a2 <- sparse_pcor2(mat_sparse_sub)
a3 <- sparse_pcor3(mat_sparse_sub)
a4 <- sparse_pcor4(mat_sparse_sub)

all.equal(a1, a2, tolerance = 1e-6)
all.equal(a3, a4, tolerance = 1e-6)
all.equal(a1, a4, tolerance = 1e-6)



res1 <- microbenchmark::microbenchmark(
  qlcmat_nofilt = sparse_pcor1(mat_sparse_keep),
  proxyc_nofilt = sparse_pcor2(mat_sparse_keep),
  qlcmat_filt = sparse_pcor3(mat_sparse_keep),
  proxyc_filt = sparse_pcor4(mat_sparse_keep),
  times = 10
)


# qlcMatrix:: faster and pre-filtering is marginally slower for both
boxplot(res1)



# testing scor time


res2 <- microbenchmark::microbenchmark(
  qlcmat_1core = sparse_scor1(mat_sparse_sub, ncores = 1),
  qlcmat_8core = sparse_scor1(mat_sparse_sub, ncores = 8),
  proxyc_1core = sparse_scor2(mat_sparse_sub, ncores = 1),
  proxyc_8core = sparse_scor2(mat_sparse_sub, ncores = 8),
  dense_1core = WGCNA::cor(mat_dense_sub, method = "spearman", use = "everything", nThreads = 1),
  dense_8core = WGCNA::cor(mat_dense_sub, method = "spearman", use = "everything", nThreads = 8),
  times = 10
)

# qlcMatrix:: faster (as expected) and both faster than dense. adding cores
# slows down on sub
boxplot(res2)



# last check of parallel on larger matrix requiring more filtering (ignore dense)

res3 <- microbenchmark::microbenchmark(
  qlcmat_1core = sparse_scor1(mat_sparse, ncores = 1),
  qlcmat_8core = sparse_scor1(mat_sparse, ncores = 8),
  proxyc_1core = sparse_scor2(mat_sparse, ncores = 1),
  proxyc_8core = sparse_scor2(mat_sparse, ncores = 8),
  times = 5
)

# parallel not worth it
boxplot(res3)
