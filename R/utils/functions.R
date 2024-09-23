## Project functions
## -----------------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(parallel)
library(Matrix)
library(qlcMatrix)
library(DescTools)



# Misc/general



# Not in

'%!in%' <- function(x, y) !('%in%'(x, y))



# Execute a function with provided args and save an RDS to path

save_function_results <- function(path, 
                                  fun,
                                  args,
                                  force_resave = FALSE) {
  
  if (!file.exists(path) || force_resave) {
    
    result <- do.call(fun, args)
    
    if (!is.null(result)) {
      saveRDS(result, path)
    }
    
    return(invisible(NULL))
  }
}



# Convert a matrix into a long and skinny df. If symmetric, only return the
# unique values.

mat_to_df <- function(mat, symmetric = TRUE, value_name = NULL) {
  
  if (symmetric) {
    df <- data.frame(
      Row = rownames(mat)[row(mat)[lower.tri(mat)]],
      Col = colnames(mat)[col(mat)[lower.tri(mat)]],
      Value = mat[lower.tri(mat)],
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Row = rownames(mat)[row(mat)],
      Col = colnames(mat)[col(mat)],
      Value = c(mat),
      stringsAsFactors = FALSE
    )
  }
  
  if (!is.null(value_name)) colnames(df)[colnames(df) == "Value"] <- value_name
  
  return(df)
}








# Loading and saving data objects
# ------------------------------------------------------------------------------



# Wrapper around data.table::fwrite() to save a matrix as a tab delim data.frame
# so that matrix rownames are saved as the first column

fwrite_mat <- function(mat, path) {
  
  fwrite(
    data.frame(mat, check.names = FALSE),
    sep = "\t",
    row.names = TRUE,
    quote = FALSE,
    verbose = FALSE,
    showProgress = FALSE,
    file = path
  )
}



# Calls data.table::fread() to read in a table where the first column corresponds
# to gene names, and convert this table to a gene x gene matrix. 
# sub_genes controls if only a subset of column genes should be loaded 

fread_to_mat <- function(path, genes, sub_genes = NULL) {
  
  if (!is.null(sub_genes)) {
    
    dat <- fread(path, sep = "\t", select = c("V1", sub_genes))
    mat <- as.matrix(dat[, -1, drop = FALSE])
    rownames(mat) <- dat[["V1"]]
    colnames(mat) <- sub_genes
    mat <- mat[genes, sub_genes, drop = FALSE]
    
  } else {
    
    dat <- fread(path, sep = "\t")
    mat <- as.matrix(dat[, -1, drop = FALSE])
    rownames(mat) <- colnames(mat) <- dat[["V1"]]
    mat <- mat[genes, genes, drop = FALSE]
  }
  
  return(mat)
}



# Uses fread() to read from path, assuming that the introduced V1 (rownames)
# column is for genes. Coerces to matrix and assigns gene rownames. Columns
# are assumed to be cell IDs (thus distinct from fread_to_mat())

read_count_mat <- function(dat_path) {
  
  dat <- suppressWarnings(fread(dat_path))
  genes <- dat[["V1"]]
  dat[["V1"]] <- NULL
  mat <- as.matrix(dat)
  rownames(mat) <- genes
  
  return(mat)
}

