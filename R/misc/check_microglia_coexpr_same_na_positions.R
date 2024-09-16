# Ensure P/S cor matrices have same NA positions 


library(tidyverse)
library(parallel)
source("R/utils/functions.R")
source("R/00_config.R")


pc_df <- read.delim(ref_mm_path, stringsAsFactors = FALSE)


mcg_p_ids <- list.files(file.path(amat_dir, "Microglia", "Mm_pcor"),
                        full.names = TRUE, pattern = "cormat")

mcg_s_ids <- list.files(file.path(amat_dir, "Microglia", "Mm_scor"),
                        full.names = TRUE, pattern = "cormat")


stopifnot(length(mcg_p_ids) == length(mcg_s_ids))




# NA positions are the same
check_l <- lapply(1:length(mcg_p_ids), function(x) {

  pmat <- fread_to_mat(mcg_p_ids[x], genes = pc_df$Symbol)
  smat <- fread_to_mat(mcg_s_ids[x], genes = pc_df$Symbol)

  # Because cor argument can change whether diag coerced to 1 for all NA vecs
  pmat <- diag_to_one(pmat)
  smat <- diag_to_one(smat)

  identical(is.na(pmat), is.na(smat))

})



# Ensuring that no matrix has all NAs
# check_l <- lapply(1:length(mcg_p_ids), function(x) {
#   pmat <- fread_to_mat(mcg_p_ids[x], genes = pc_df$Symbol)
#   !all(is.na(pmat))
# })




stopifnot(all(unlist(check_l)))


message("All good")
