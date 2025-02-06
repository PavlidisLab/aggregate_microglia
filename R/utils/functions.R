## Project functions
## -----------------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(parallel)
library(Matrix)
library(qlcMatrix)
# library(aggtools)
library(preprocessCore)
library(matrixStats)



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
  mat <- Matrix(as.matrix(dat), sparse = TRUE)
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



# Create a list that summarizes gene count levels across all datasets in dat_l
# dat_l: list of single cell datasets (count matrix and metadata)
# Assumes dat_l contains a gene x cell CPM matrix called "Mat"
# Assumes all matrices have the same genes and ordering
# Note: matrixStats requires coercing sparse count matrices to dense.
# TODO: consider checks

summarize_gene_counts <- function(dat_l) {
  
  genes <- rownames(dat_l[[1]]$Mat)
  ids <- names(dat_l)
  stopifnot(length(genes) > 0, length(ids) > 0) 
  
  # Log transform all CPM matrices first
  dat_l <- lapply(dat_l, function(x) as.matrix(log2(x$Mat + 1)))  
  gc(verbose = FALSE)
  
  # Get average, SD, and CV of log CPM counts and bind into a matrix
  avg <- do.call(cbind, lapply(dat_l, rowMeans))
  sd <- do.call(cbind, lapply(dat_l, rowSds))
  cv <- sd / avg
  
  # Quantile normalize the averaged profiles
  qn_avg <- preprocessCore::normalize.quantiles(avg, keep.names = TRUE)
  
  # Rank product of the averaged profiles
  rank_avg <- aggtools::colrank_mat(avg)
  rp_avg <- rowSums(log(rank_avg)) / length(genes)
  
  # Get binary gene measurement status (min 20 cells with at least one count)
  is_measured <- function(mat) rowSums(mat > 0) >= 20
  msr <- do.call(cbind, lapply(dat_l, is_measured))
  
  # Summary dataset of point estimates for each gene collapsed across datasets
  
  summ_df <- data.frame(
    Symbol = genes,
    Avg = rowMeans(avg, na.rm = TRUE),
    QN_avg = rowMeans(qn_avg, na.rm = TRUE),
    SD = rowMeans(sd, na.rm = TRUE),
    CV = rowMeans(cv, na.rm = TRUE),
    RP = rank(rp_avg),
    N_msr = rowSums(msr)
  )
  
  # Genes that are measured in a minimum proportion of datasets
  cutoff <- floor(length(ids) * (1/3))
  filter_genes <- filter(summ_df, N_msr >= cutoff) %>% pull(Symbol)
  
  return(list(
    Avg = avg,
    QN_Avg = qn_avg,
    SD = sd,
    CV = cv,
    Msr = msr,
    Summ_df = summ_df,
    Filter_genes = filter_genes
  ))
}



# Functions for QC/preprocessing count matrices and metadata
# ------------------------------------------------------------------------------



# Gives a matrix with ENSEMBL IDs as rownames, return the matrix with the 
# corresponding gene symbols as rownames. Blank gene symbols are removed

ensembl_to_symbol <- function(mat, ensembl_df) {
  
  stopifnot(c("Gene_ID", "Symbol") %in% colnames(ensembl_df))
  
  ids <- intersect(pc$Gene_ID, rownames(mat))
  
  if (length(ids) == 0) stop("No common ENSEMBL IDs in rownames of matrix")
  
  common_genes <- data.frame(
    ID = ids,
    Symbol = pc$Symbol[match(ids, pc$Gene_ID)]) %>% 
    filter(Symbol != "")
  
  # dupl_genes <- common_genes$Symbol[which(duplicated(common_genes$Symbol))]
  
  mat <- mat[common_genes$ID, ]
  rownames(mat) <- common_genes$Symbol
  
  return(mat)
}



# Assumes that mat is sparse gene x cell count matrix. Filters the matrix for 
# unique gene symbols in pcoding_df, and fills the missing genes as 0s. 

get_pcoding_only <- function(mat, pcoding_df) {
  
  stopifnot("Symbol" %in% colnames(pcoding_df))
  
  genes <- unique(pcoding_df$Symbol)
  common <- intersect(rownames(mat), genes)
  missing <- setdiff(genes, rownames(mat))
  
  if (length(common) == 0) stop("No common symbols in rownames of mat")
  
  pc_mat <- mat[common, ]
  pc_mat <- rbind(pc_mat, Matrix(0, nrow = length(missing), ncol = ncol(mat)))
  rownames(pc_mat) <- c(common, missing)
  pc_mat <- pc_mat[genes, ]
  
  return(pc_mat)
}



# Given a vector of genes that have either common gene symbols or ensembl
# ids, return a subst of gene_vec only containing the mitochondrial genes. 
# Assumes gene_vec has only mouse or human symbols/ensembl IDs.

get_mt_genes <- function(gene_vec,
                         mt_path = "/home/amorin/Data/Metadata/mitochondrial_genes_all.tsv") {
  
  mt_table <- read.delim(mt_path, stringsAsFactors = FALSE)
  mt_genes <- gene_vec[gene_vec %in% c(mt_table$Gene_stable_ID, mt_table$Gene_name)]
  
  return(mt_genes)
}



# This adds columns to metadata: the number of total UMI counts for each 
# cell/column of mat, the number of non-zero expressing genes, and the RNA
# novelty/compexity, which is the ratio of the log10 gene counts to log10 umi 
# counts. It additionally adds ratio of mitochondrial if available. 
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

add_count_info <- function(mat, meta) {
  
  mt_genes <- get_mt_genes(rownames(mat))
  
  meta <- meta %>% 
    mutate(
      UMI_counts = colSums(mat),
      Gene_counts = colSums(mat > 0),
      RNA_novelty = log10(Gene_counts) / log10(UMI_counts)
    )
  
  if (length(mt_genes) > 0) {
    mt_ratio <- colSums(mat[mt_genes, , drop = FALSE]) / meta$UMI_counts
    meta$MT_ratio = mt_ratio
  }
  
  # Remove cell x gene features of this type, if present
  meta <- meta[, !(colnames(meta) %in% c("nFeature_RNA", "nCount_RNA"))]
  
  return(meta)
}



# This subsets mat to remove cells that fail any of the filters laid out in: 
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html
# The RNA novelty filter is relaxed for Smart-seq runs as more reads can go to
# genes like mitochondrial https://pubmed.ncbi.nlm.nih.gov/33662621/

rm_low_qc_cells <- function(mat, 
                            meta,
                            min_counts = 500,
                            min_genes = 250,
                            min_novelty = NULL,
                            max_mt_ratio = 0.2) {
  
  keep_id <- meta %>%
    mutate(Is_smartseq = str_detect(str_to_lower(assay), "smart-seq")) %>%
    filter(
      UMI_counts >= min_counts,
      Gene_counts >= min_genes,
      ifelse(Is_smartseq, RNA_novelty > 0.5, RNA_novelty > 0.8)) %>%
    pull(ID)
  
  if ("MT_ratio" %in% colnames(meta)) {
    keep_id <- intersect(keep_id, filter(meta, MT_ratio < max_mt_ratio)$ID)
  }
  
  if (length(keep_id) == 0) stop("No remaining cells after filtering")
  
  return(mat[, keep_id])
}

