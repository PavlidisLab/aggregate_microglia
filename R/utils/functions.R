## Project functions
## -----------------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(parallel)
library(Matrix)
library(qlcMatrix)
library(aggtools)



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



# Preparing data
# ------------------------------------------------------------------------------


# Create a filtered data frame of the cell types matching the input pattern

# pattern: character vector of terms to search for within cell types
# ct_list: list of dataset data.frames containing the cell types and their counts
# sc_meta: data frame containing metadata information of datasets in ct_list
# return: a data frame of the specified cell types and their counts, joined to sc_meta

create_celltype_df <- function(pattern, ct_list, sc_meta) {
  
  keep_meta_cols <- c("ID", "Species", "Platform", "Design", "GEO_link", "Data_link")
  
  # Return list of relevant cell types matching string search
  filt_ct_list <- lapply(ct_list, function(x) {
    ct_vec <- str_to_lower(x$Ct_count$Cell_type)
    match_indices <- str_detect(ct_vec, pattern)
    x$Ct_count[match_indices, , drop = FALSE]
  }) %>% 
    discard(~nrow(.) == 0)
  
  # Bind into a df and join with relevant details from meta
  ct_df <- bind_rows(filt_ct_list, .id = "ID") %>% 
    left_join(., sc_meta[, keep_meta_cols], by = "ID")

  return(ct_df)
}



# Return a list of single cell data subset to a cell type based on provided input_df
# input_df: data.frame with dataset ID, data path, and cell type to filter
# TODO: consider if aggtools::load_scdat() should be isolated

load_and_filter_celltype <- function(input_df, verbose = TRUE) {
  
  ids <- unique(input_df$ID)
  
  dat_l <- lapply(ids, function(id) {
    
    if (verbose) message(paste(id, Sys.time()))
    
    df <- filter(input_df, ID == id)
    path <- unique(df$Path)
    dat <- aggtools::load_scdat(path)
    meta <- filter(dat$Meta, Cell_type %in% df$Cell_type)
    mat <- dat$Mat[, meta$ID, drop = FALSE]
    
    list(Mat = mat, Meta = meta)
    
  })
  names(dat_l) <- ids
  
  return(dat_l)
}




# Add additional columns to meta on the gene measurement of data in dat_l
# dat_l: list of single cell datasets (count matrix and metadata)
# sc_meta: metadata table tracking the single cell experiments
# -- Added columns:
# N_msr_prefilt: Genes that had at least one count in any cell.
# N_msr_posfilt: Genes that had at least one count in 20+ cells (for coexpression)
# Median_UMI: Median UMI (or UMI-like) per cell -- measure of seq. depth

add_counts_to_meta <- function(dat_l, sc_meta) {
  
  # Getting gene counts across all count matrices in dat_l
  count_l <- lapply(dat_l, function(x) {
    
    prefilt <- t(x$Mat)
    postfilt <- zero_sparse_cols(t(x$Mat))
    
    data.frame(
      N_msr_prefilt = sum(colSums(prefilt) != 0),
      N_msr_postfilt = sum(colSums(postfilt) != 0),
      Median_UMI = median(x$Meta$UMI_counts)
    )
  })
  
  # Join counts with sc_meta, removing N_genes (as calc'd over all cell types)
  sc_meta <- bind_rows(count_l) %>% 
    mutate(ID = names(dat_l)) %>% 
    left_join(., sc_meta, by = "ID") %>% 
    relocate(c(N_msr_prefilt, N_msr_postfilt, Median_UMI), .after = N_cells)
  
  return(sc_meta)
}



# A single dataset may have multiple relevant cell types, encoded as distinct
# rows in the metadata. Collapse these rows into one such that each row of the 
# metadata is a unique data ID

collapse_dupl_ids <- function(sc_meta) {
  
  dupl_ids <- unique(sc_meta$ID[duplicated(sc_meta$ID)])
  
  dedupl_l <- lapply(dupl_ids, function(x) {
    
    df <- filter(sc_meta, ID == x)
    ct <- paste(df$Cell_type, collapse = ";")
    platform <- paste(unique(df$Platform), collapse = ";")
    ncells <- sum(df$N_cells)
    
    df$Cell_type[1] <- ct
    df$Platform[1] <- platform
    df$N_cells[1] <- ncells
    
    df[1, ]
    
  })
  
  dedupl_meta <- bind_rows(dedupl_l) %>%
    rbind(., filter(sc_meta, ID %!in% dupl_ids)) %>% 
    arrange(match(ID, sc_meta$ID))
  
  return(dedupl_meta)
}
