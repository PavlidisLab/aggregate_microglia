## This script is essentially a copy of aggtools, save that it allows for
## calculation of Spearman's correlation. This requires (and thus relaxes checks)
## dense matrices.


#' Rank the columns of a matrix such that rank=1 is the highest positive value
#'
#' @param mat A numeric matrix
#' @param ties_arg A character indicating how ties are assigned, as fed into
#' the ties argument of the base rank() function. "min" is the assumed default.
#' @param na_arg A character indicating how NAs are handled, as fed into
#' the na.last argument of the base rank() function. "keep" is the assumed default.
#'
#' @return Return a matrix of integer ranks the same dimension as input matrix
#' @export
#'
#' @examples
colrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {

  rank_mat <- apply(-mat, 2, rank, ties.method = ties_arg, na.last = na_arg)
  return(rank_mat)
}



#' Rank the entire matrix jointly such that rank=1 is the highest positive value
#'
#' @param mat A numeric matrix
#' @param ties_arg A character indicating how ties are assigned, as fed into
#' the ties argument of the base rank() function. "min" is the assumed default.
#' @param na_arg A character indicating how NAs are handled, as fed into
#' the na.last argument of the base rank() function. "keep" is the assumed default.
#'
#' @return Return a matrix of integer ranks the same dimension as input matrix
#' @export
#'
#' @examples
allrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {

  rank_mat <- rank(-mat, ties.method = ties_arg, na.last = na_arg)
  rank_mat <- matrix(rank_mat, nrow = nrow(mat), ncol = ncol(mat))
  rownames(rank_mat) <- rownames(mat)
  colnames(rank_mat) <- colnames(mat)

  return(rank_mat)
}


#' Fisher Z transformation
#'
#' @param cor A numeric correlation value (single element, array, matrix)
#'
#' @return The transformed correlations
#' @seealso DescTools::FisherZ()
#' @export
#'
#' @examples
fisherz <- function(cor) {

  stopifnot(is.numeric(cor), all(cor >= -1 & cor <= 1))

  fz <- 0.5 * log((1 + cor) / (1 - cor))
  return(fz)
}




# Wraps to qlcMatrix::sparseCor() to make resulting matrix inherit input names

#' Column-wise Pearson's correlation for sparse matrices
#'
#' @importFrom qlcMatrix corSparse
#' @import Matrix
#' @param mat A sparse numeric m by n Matrix
#' @return A dense n by n matrix of the Pearson's correlation between the
#' columns of mat
#' @export
#'
#' @examples
sparse_pcor <- function(mat) {

  stopifnot(inherits(mat, "dgCMatrix"))

  cmat <- qlcMatrix::corSparse(mat)
  colnames(cmat) <- rownames(cmat) <- colnames(mat)

  return(cmat)
}



#' Wrapper for calling sparse Pearson or Spearman correlation
#'
#' @param mat A sparse numeric m by n Matrix
#' @param cor_method One of "pearson" or "spearman"
#'
#' @return A dense n by n matrix of the correlations between the columns of mat
#' @export
#'
#' @examples
calc_sparse_correlation <- function(mat, cor_method) {

  # stopifnot(inherits(mat, "dgCMatrix"))
  stopifnot(cor_method %in% c("pearson", "spearman"))

  cmat <- if (cor_method == "pearson") {
    sparse_pcor(mat)
  } else {
    mat <- as.matrix(mat)
    WGCNA::cor(mat, method = "spearman")
  }

  return(cmat)
}



#' Set NAs in matrix to 0
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all values are preserved save
#' for NAs assigned to 0
#' @export
#'
#' @examples
na_to_zero <- function(mat) {

  stopifnot(is.matrix(mat))
  mat[is.na(mat)] <- 0
  return(mat)
}



#' Set diag of a square matrix to 1
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all values are preserved save
#' for the diagonals assigned to 1
#' @export
#'
#' @examples
diag_to_one <- function(mat) {

  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  diag(mat) <- 1
  return(mat)
}



#' Set upper triangle of a square matrix to NA
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all values are preserved save
#' for the values above the diagonal set to NA and the diagonal is preserved
#' @export
#'
#' @examples
uppertri_to_na <- function(mat) {

  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  mat[upper.tri(mat)] <- NA
  return(mat)
}



#' Replace the upper triangle of square matrix with its lower triangle
#'
#' @param mat A matrix
#'
#' @return A matrix of the same dimensions where all upper triangular values are
#' replaced by its lower triangle values
#' @export
#'
#' @examples
lowertri_to_symm <- function(mat) {

  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  mat[upper.tri(mat)] <-  t(mat)[upper.tri(mat)]
  return(mat)
}



#' Set sparse columns to 0
#'
#' If a col of mat has fewer non-zero elements than min_cell, set that col to
#' 0. This is done to produce an NA during correlation, instead of allowing cors
#' derived from overly sparse genes/columns.
#' n=20 default used from https://doi.org/10.1093/bioinformatics/btv118
#'
#' @import Matrix
#' @param mat A sparse matrix
#' @param min_cell A non-negative integer corresponding to how many cells of a
#' cell type must have at least one count for a gene to be considered "measured"
#'
#' @return A matrix of the same dimensions where all sparse columns are set to 0
#' @export
#'
#' @examples
zero_sparse_cols <- function(mat, min_cell = 20) {

  # stopifnot(inherits(mat, "dgCMatrix"),
  stopifnot(is.numeric(min_cell),
            min_cell >= 0)

  nonzero_cells <- colSums(mat != 0)
  filt_genes <- nonzero_cells < min_cell
  if (any(filt_genes)) mat[, filt_genes] <- 0

  return(mat)
}



#' Prepare a cell type matrix for coexpression
#'
#' Subset a sparse gene by cell count matrix to cells annotated to the input
#' cell_type and genes contained in pc_df, transpose the resulting matrix to
#' cells by genes, and set sparse genes/columns to 0

#' @param mat A sparse gene by cell matrix
#' @param meta A data frame that maps cell IDs to cell types
#' @param pc_df A data frame of unique protein coding gene symbols
#' @param cell_type A character of the cell type to subset mat on
#' @param min_cell A non-negative integer corresponding to how many cells of a
#' cell type must have at least one count for a gene to be considered "measured"
#'
#' @return A sparse cell by gene matrix
#' @export
#'
#' @examples
prepare_celltype_mat <- function(mat, meta, pc_df, cell_type, min_cell = 20) {

  # stopifnot(inherits(mat, "dgCMatrix"),
  stopifnot(c("ID", "Cell_type") %in% colnames(meta),
            "Symbol" %in% colnames(pc_df),
            cell_type %in% meta[["Cell_type"]])

  ct_meta <- meta[meta$Cell_type %in% cell_type, "ID", drop = FALSE]
  ids <- ct_meta[["ID"]]

  common_ids <- intersect(ids, colnames(mat))
  diff_ids <- setdiff(ids, colnames(mat))
  common_genes <- intersect(rownames(mat), pc_df$Symbol)
  diff_genes <- setdiff(pc_df$Symbol, rownames(mat))

  if (length(common_ids) == 0) stop("No common IDs between metadata and count matrix")
  if (length(diff_ids) > 0) stop("Not all IDs in metadata are in count matrix")
  if (length(common_genes) == 0) stop("No common genes between pc_df and count matrix")
  if (length(diff_genes) > 0) stop("Not all genes in pc_df are found in count matrix")

  ct_mat <- t(mat[pc_df$Symbol, ids])
  ct_mat <- zero_sparse_cols(ct_mat, min_cell)

  return(ct_mat)
}



#' Initiate aggregate matrix
#'
#' Initates a gene by gene matrix of 0s for holding aggregate correlations.
#' pc_df assumed to be protein coding table with "Symbol" as a column and unique
#' genes within Symbol.
#'
#' @param pc_df A data frame of unique protein coding gene symbols
#'
#' @return A gene by gene matrix of 0s with dimensions equal to the count of
#' symbols in pc_df
#' @export
#'
#' @examples
init_agg_mat <- function(pc_df) {

  stopifnot("Symbol" %in% colnames(pc_df),
            anyDuplicated(pc_df[["Symbol"]]) == 0)

  amat <- matrix(0, nrow = nrow(pc_df), ncol = nrow(pc_df))
  rownames(amat) <- colnames(amat) <- pc_df[["Symbol"]]
  return(amat)
}



#' Count the NAs in cmat and increment the corresponding indices of na_mat
#'
#' @param cmat A gene by gene correlation matrix
#' @param na_mat A gene by gene matrix that tracks NA presence across multiple
#' correlation matrices
#'
#' @return na_mat with the indices corresponding to NAs in cmat increased by 1
#' @export
#'
#' @examples
increment_na_mat <- function(cmat, na_mat) {

  na_ix <- which(is.na(cmat), arr.ind = TRUE)
  na_mat[na_ix] <- na_mat[na_ix] + 1
  return(na_mat)
}



#' Transform correlation matrix
#'
#' Take an input correlation matrix and transform it to prepare for aggregation.
#' Transformation entails all ranking, column ranking, or Fisher's Z.
#' Each methods impute NA cors to 0 and ensure that the diagonal equals 1.
#'
#' allrank: make the matrix tri. to prevent double ranking symmetric elements
#' and then jointly rank the matrix (lower rank = positive cor).
#'
#' colrank: rank each column separately (lower rank = positive cor).
#'
#' FZ: perform Fisher's Z transform on the raw correlations

#' @param cmat A gene by gene correlation matrix
#' @param agg_method One of "allrank", "colrank", or "FZ"
#'
#' @return A gene by gene matrix of transformed correlations
#' @export
#'
#' @examples
transform_correlation_mat <- function(cmat, agg_method) {

  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"),
            is.matrix(cmat), identical(rownames(cmat), colnames(cmat)))

  cmat <- na_to_zero(cmat)
  cmat <- diag_to_one(cmat)

  # Coerce values slightly above 1 or below -1 due to floating-point precision
  cmat[cmat > 1] <- 1
  cmat[cmat < -1] <- -1

  if (agg_method == "allrank") {

    # Rank standardizing: larger values more important then divide by max rank
    cmat <- uppertri_to_na(cmat)  # Prevent ranking symmetric elements
    cmat <- allrank_mat(-cmat, ties_arg = "min")
    cmat <- cmat / max(cmat, na.rm = TRUE)

  } else if (agg_method == "colrank") {

    cmat <- colrank_mat(-cmat, ties_arg = "min")  # Larger values more important
    cmat <- apply(cmat, 2, function(x) x / max(x)) # Rank std.

  } else if (agg_method == "FZ") {

    cmat <- fisherz(cmat)
  }

  return(cmat)
}



#' Finalize aggregate matrix
#'
#' Format the final summed aggregate correlations depending on the aggregate
#' strategy.
#'
#' all_rank: jointly re-rank the summed ranks and standardize into \[0, 1\], then
#' convert back to a symmetric matrix for ease of downstream operations
#'
#' col_rank: re-rank the summed ranks for each column separately and standardize
#' into \[0, 1\]
#'
#' FZ: divide each element by its count of its measured (i.e., non-NA)
#' observations
#'
#' @param amat A gene by gene matrix tracking the aggregate correlations
#' @param agg_method One of "allrank", "colrank", or "FZ"
#' @param n_celltypes The count of cell types that went into the aggregated
#' correlation
#' @param na_mat A gene by gene matrix tracking the count of NAs across the
#' correlation matrices
#'
#' @return A gene by gene matrix of the formatted aggregated correlations
#' @export
#'
#' @examples
finalize_agg_mat <- function(amat, agg_method, n_celltypes, na_mat) {

  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"),
            is.matrix(amat), identical(rownames(amat), colnames(amat)))

  if (agg_method == "allrank") {

    # Rank standardizing: larger values more important then divide by max rank
    amat <- allrank_mat(-amat, ties_arg = "min")
    amat <- amat / max(amat, na.rm = TRUE)
    amat <- diag_to_one(amat)  # Make aggregate symmetric and self cor = 1
    amat <- lowertri_to_symm(amat)

  } else if (agg_method == "colrank") {

    amat <- colrank_mat(-amat, ties_arg = "min")  # Larger values more important
    amat <- apply(amat, 2, function(x) x / max(x)) # Rank std.

  } else if (agg_method == "FZ") {

    amat <- amat / (n_celltypes - na_mat)

  }

  return(amat)
}



#' Aggregate cell type correlation matrices within a single dataset
#'
#' Iterates through each unique cell type in meta, subsetting mat to cells of
#' the given type and generating a cell type coexpression matrix for each,
#' which are aggregated into a single matrix.
#'
#' https://pubmed.ncbi.nlm.nih.gov/27165153/
#' https://pubmed.ncbi.nlm.nih.gov/25717192/
#' https://pubmed.ncbi.nlm.nih.gov/34015329/
#'
#' all_rank: each cell type coexpression matrix is ranked across the entire
#' matrix, then summed and rank standardized into \[0, 1\].
#'
#' col_rank: each cell type coexpression matrix is ranked column-wise, then
#' summed and rank standardized into \[0, 1\].
#'
#' FZ: Fisher's Z transformation is applied to each cell type coexpression
#' matrix, which are then summed and divided element-wise by the count of times
#' a given gene-gene pair was co-measured across all cell types.
#'
#' @param mat A sparse gene by cell matrix
#' @param meta A data frame that maps cell IDs to cell types
#' @param pc_df A data frame of unique protein coding gene symbols
#' @param cor_method One of "pearson" or "spearman"
#' @param agg_method One of "allrank", "colrank", or "FZ"
#' @param min_cell A non-negative integer corresponding to how many cells of a
#' cell type must have at least one count for a gene to be considered "measured"
#' and thus viable for calculating correlation
#' @param verbose Logical determining if the current cell type and time should
#' be displayed as a message during iteration
#'
#' @return A list of two gene by gene matrices, the first of which is the
#' aggregate coexpression matrix and the second of which is the NA tracking
#' matrix
#' @export
#'
#' @examples
aggr_coexpr_single_dataset <- function(mat,
                                       meta,
                                       pc_df,
                                       cor_method,
                                       agg_method,
                                       min_cell = 20,
                                       verbose = TRUE) {

  # stopifnot(inherits(mat, "dgCMatrix"),
  stopifnot(c("ID", "Cell_type") %in% colnames(meta),
            cor_method %in% c("pearson", "spearman"),
            agg_method %in% c("allrank", "colrank", "FZ"))

  cts <- unique(meta[["Cell_type"]])
  n_cts <- length(cts)

  # Matrices of 0s for tracking aggregate correlation and count of NAs
  amat <- init_agg_mat(pc_df)
  na_mat <- amat

  for (ct in cts) {

    if (verbose) message(paste(ct, Sys.time()))

    # Get count matrix for current cell type, coercing low count genes to 0
    ct_mat <- prepare_celltype_mat(mat = mat,
                                   meta = meta,
                                   pc_df = pc_df,
                                   cell_type = ct,
                                   min_cell = min_cell)

    no_msr <- all(ct_mat == 0)

    if (no_msr) {
      message(paste(ct, "skipped due to insufficient counts"))
      na_mat <- na_mat + 1
      next()
    }

    # Get cell-type cor matrix and increment count of NAs before imputing to 0
    cmat <- calc_sparse_correlation(ct_mat, cor_method)
    na_mat <- increment_na_mat(cmat, na_mat)

    # Transform raw correlation matrix, add to aggregate and clean up
    cmat <- transform_correlation_mat(cmat, agg_method)
    amat <- amat + cmat
    rm(cmat, ct_mat)
    gc(verbose = FALSE)

  }

  # Final format of aggregate matrix and return along with the tracked NA mat
  amat <- finalize_agg_mat(amat, agg_method, n_cts, na_mat)

  return(list(Agg_mat = amat,
              NA_mat = na_mat))

}




#' Load a single cell dataset
#'
#' Assumes path leads to an RDS of a list with 2 named elements: "Mat", the
#' sparse count matrix, and "Meta", a metadata data frame mapping cell IDs to an
#' annotated cell type
#'
#' @param path A character corresponding to the path of the dataset
#'
#' @return A list of the single cell count matrix and cell type metadata
#' @export
#'
#' @examples
load_scdat <- function(path) {

  stopifnot(is.character(path), length(path) == 1, file.exists(path))

  dat <- readRDS(path)
  meta <- dat[["Meta"]]
  mat <- dat[["Mat"]]
  stopifnot(all(colnames(mat) %in% meta[["ID"]]))

  return(list(Mat = mat, Meta = meta))
}



#' Aggregate cell type correlation matrices across datasets
#'
#' Iteratively loads datasets in input_dat, subsetting the count matrix of each
#' to the specified cell type and generating a cell type coexpression matrix,
#' which are then all aggregated into a single matrix. Note that if multiple
#' cell types are listed for a single dataset ID, these cell types are
#' collapsed into one matrix for calculating coexpression and tracking measurement.
#'
#' https://pubmed.ncbi.nlm.nih.gov/27165153/
#' https://pubmed.ncbi.nlm.nih.gov/25717192/
#' https://pubmed.ncbi.nlm.nih.gov/34015329/
#'
#' all_rank: each cell type coexpression matrix is ranked across the entire
#' matrix, then summed and rank standardized into \[0, 1\].
#'
#' col_rank: each cell type coexpression matrix is ranked column-wise, then
#' summed and rank standardized into \[0, 1\].
#'
#' FZ: Fisher's Z transformation is applied to each cell type coexpression
#' matrix, which are then summed and divided element-wise by the count of times
#' a given gene-gene pair was co-measured across all cell types.
#'
#' @param input_df A dataframe of the dataset paths and cell types to consider
#' @param pc_df A data frame of unique protein coding gene symbols
#' @param cor_method One of "pearson" or "spearman"
#' @param agg_method One of "allrank", "colrank", or "FZ"
#' @param min_cell A non-negative integer corresponding to how many cells of a
#' cell type must have at least one count for a gene to be considered "measured"
#' and thus viable for calculating correlation
#' @param verbose Logical determining if the current cell type and time should
#' be displayed as a message during iteration
#'
#' @return A list of two gene by gene matrices, the first of which is the
#' aggregate coexpression matrix and the second of which is the NA tracking
#' matrix
#' @export
#'
#' @examples
aggr_coexpr_multi_dataset <- function(input_df,
                                      pc_df,
                                      cor_method,
                                      agg_method,
                                      min_cell = 20,
                                      verbose = TRUE) {

  stopifnot(all(c("ID", "Cell_type", "Path") %in% colnames(input_df)),
            cor_method %in% c("pearson", "spearman"),
            agg_method %in% c("allrank", "colrank", "FZ"))

  data_ids <- unique(input_df[["ID"]])
  n_dat <- length(data_ids) # All cell types for a dataset are collapsed

  # Matrices of 0s for tracking aggregate correlation and count of NAs
  amat <- init_agg_mat(pc_df)
  na_mat <- amat

  for (id in data_ids) {

    if (verbose) message(paste(id, Sys.time()))

    ct <- input_df[input_df$ID == id, "Cell_type"]
    dat_path <- unique(input_df[input_df$ID == id, "Path"])

    # Load dataset and get count matrix for current cell type
    dat <- load_scdat(dat_path)

    ct_mat <- prepare_celltype_mat(mat = dat[["Mat"]],
                                   meta = dat[["Meta"]],
                                   pc_df = pc_df,
                                   cell_type = ct,
                                   min_cell = min_cell)

    # Check if filtering removed all genes
    no_msr <- all(ct_mat == 0)

    if (no_msr) {
      message(paste(ct, "skipped due to insufficient counts"))
      na_mat <- na_mat + 1
      next()
    }

    # Get cell-type cor matrix and increment count of NAs before imputing to 0
    cmat <- calc_sparse_correlation(ct_mat, cor_method)
    na_mat <- increment_na_mat(cmat, na_mat)

    # Transform raw correlation matrix, add to aggregate and clean up
    cmat <- transform_correlation_mat(cmat, agg_method)
    amat <- amat + cmat
    rm(cmat, ct_mat)
    gc(verbose = FALSE)

  }

  # Final format of aggregate matrix and return along with the tracked NA mat

  amat <- finalize_agg_mat(amat, agg_method, n_dat, na_mat)

  return(list(Agg_mat = amat,
              NA_mat = na_mat))
}
