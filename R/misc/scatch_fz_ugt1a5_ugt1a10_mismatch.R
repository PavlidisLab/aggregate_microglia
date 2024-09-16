# Inspecting aggregate FZ gene pairs that had unexpected (non-diag) Inf values

source("R/utils/dev_functions.R")
source("R/00_config.R")

#

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# List of cell types
ct_l <- readRDS(celltype_list_path)

pc_df <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# Save out
out_dir <- "/space/scratch/amorin/TR_singlecell/Microglia/Mm_pcor"
dir.create(out_dir, showWarnings = FALSE)
out_path <- file.path(out_dir, "check_ugt1a5.RDS")


keep_genes <- c("Ugt1a5",
                "Ugt1a10",
                "Ugt1a9",
                "Ugt1a2",
                "Ugt1a1",
                "Ugt1a6b",
                "Ugt1a6a",
                "Ugt1a7c",
                "H3c13",
                "H3c14")


# Extract IDs of specific cell type of interest
# ------------------------------------------------------------------------------


ct <- "microglia|mcg|mgl|microg"  


ct_df <- lapply(ct_l, function(x) {
  ct_vec <- str_to_lower(x$Ct_count$Cell_type)
  ct_which <- str_detect(ct_vec, ct)
  
  if (sum(ct_which) == 0) {
    return(NA)
  }
  
  x$Ct_count[ct_which, ]
  
})


ct_df <- ct_df[!is.na(ct_df)]



ct_df <- data.frame(
  do.call(rbind, ct_df)
) %>% 
  rownames_to_column(var = "ID") %>% 
  mutate(ID = str_replace(ID, "\\.[:digit:]+$", "")) %>% 
  left_join(., sc_meta[, c("ID", "Species")], by = "ID")


ct_df <- filter(ct_df, Species == "Mouse")

# Add # GSE118020: Micro (micro grep grabs non-microglia) 
ct_df <- rbind(ct_df, c("GSE118020", "Micro", 3230, "Mouse"))


ct_df$Path <- paste0(amat_dir, ct_df$ID, "/",  ct_df$ID, "_clean_mat_and_meta_CPM.RDS")


#

# TODO: check before merge

init_agg_mat <- function(pc_df) {
  
  amat <- matrix(0, nrow = nrow(pc_df), ncol = nrow(pc_df))
  rownames(amat) <- colnames(amat) <- pc_df$Symbol
  return(amat)
}



# Assumes path leads to an RDS of a list with 2 elements: the sparse count
# matrix and the metadata mapping cell IDs to cell type

load_dataset <- function(path) {
  
  dat <- readRDS(path)
  meta <- dat$Meta
  mat <- dat$Mat
  stopifnot(identical(colnames(mat), meta$ID))
  
  return(list(Mat = mat, Meta = meta))
}



cor_method = "pearson"
agg_method = "FZ"
min_cell = 20
verbose = TRUE



# TODO:

aggr_coexpr_across_datasets <- function(ct_df,
                                        pc_df,
                                        cor_method = "pearson",
                                        agg_method = "FZ",
                                        min_cell = 20,
                                        verbose = TRUE) {
  
  stopifnot(cor_method %in% c("pearson", "spearman"))
  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"))
  
  data_ids <- unique(ct_df[["ID"]])
  n_cts <- length(data_ids) # All cell types for a dataset are collapsed
  
  # Matrices of 0s for tracking aggregate correlation and count of NAs
  
  amat <- init_agg_mat(pc_df)
  na_mat <- amat  
  
  cor_l <- lapply(data_ids, function(id) {
    
    if (verbose) message(paste(id, Sys.time()))
    
    ct <- filter(ct_df, ID == id)[["Cell_type"]]
    dat_path <- unique(filter(ct_df, ID == id)[["Path"]])
    
    # Load dataset and get count matrix for current cell type
    
    dat = load_dataset(dat_path)
    
    ct_mat <- prepare_celltype_mat(mat = dat$Mat, 
                                   meta = dat$Meta, 
                                   cell_type = ct, 
                                   min_count = min_cell)
    
    ct_mat <- ct_mat[, keep_genes]
    
    # Check if filtering removed all genes
    
    # no_msr <- all(ct_mat == 0)
    
    # if (no_msr) {
    #   message(paste(ct, "skipped due to insufficient counts"))
    #   na_mat <- na_mat + 1
    #   next()
    # }

    # Get cell-type cor matrix and increment count of NAs before imputing to 0
    
    cmat <- calc_sparse_correlation(ct_mat, cor_method)
    # na_mat <- increment_na_mat(cmat, na_mat)
    
    # Transform raw correlation matrix, add to aggregate and clean up
    
    fzmat <- transform_correlation_mat(cmat, agg_method)

    list(Cor = cmat, FZ = fzmat)
    
    
  })
  
  # Final format of aggregate matrix and return along with the tracked NA mat
  
  amat <- finalize_agg_mat(amat, agg_method, n_cts, na_mat)
  return(list(Agg_mat = amat, NA_mat = na_mat))
}




cor_l
check1 <- "Ugt1a5"
check2 <- "Ugt1a10"
cor_vec <- unlist(lapply(cor_l, function(x) x$Cor[check1, check2]))
fz_vec <- unlist(lapply(cor_l, function(x) x$FZ[check1, check2]))
n_na <- sum(is.na(cor_vec))


mean(cor_vec)
mean(fz_vec)



# cell_cor <- calc_sparse_correlation(t(ct_mat), cor_method = "pearson")
