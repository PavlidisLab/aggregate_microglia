# Generating aggregate correlation for a specific cell type
# ------------------------------------------------------------------------------

source("R/00_config.R")
library(aggtools)
library(tidyverse)


# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# List of cell types
ct_l <- readRDS(celltype_list_path)

pc_df <- read.delim(ref_mm_path, stringsAsFactors = FALSE)





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


# aggr_coexpr_multi_dataset()
# ------------------------------------------------------------------------------


input_df <- ct_df
cor_method <- "pearson"
agg_method <- "FZ"
min_cell = 20
verbose <- TRUE


data_ids <- unique(input_df[["ID"]])
n_dat <- length(data_ids)

amat <- init_agg_mat(pc_df)
na_mat <- amat

# Within one iteration

id <- "GSE160512"
if (verbose) message(paste(id, Sys.time()))

ct <- input_df[input_df$ID == id, "Cell_type"]
dat_path <- unique(input_df[input_df$ID == id, "Path"])


dat <- load_scdat(dat_path)


ct_mat <- prepare_celltype_mat(mat = dat[["Mat"]],
                               meta = dat[["Meta"]],
                               cell_type = ct,
                               min_cell = min_cell)

no_msr <- all(ct_mat == 0)


if (no_msr) {
  message(paste(ct, "skipped due to insufficient counts"))
  na_mat <- na_mat + 1
  # next()
}


cmat <- calc_sparse_correlation(ct_mat, cor_method)
na_mat <- increment_na_mat(cmat, na_mat)


cmat <- transform_correlation_mat(cmat, agg_method)
amat <- amat + cmat
# rm(cmat, ct_mat)
# gc(verbose = FALSE)


# Outside iteration

amat <- finalize_agg_mat(amat, agg_method, n_dat, na_mat)
