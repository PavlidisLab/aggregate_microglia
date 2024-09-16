# Does one of P/S cor have more negative cors


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



pmat_track <- matrix(0, ncol = length(pc_df$Symbol), nrow = length(pc_df$Symbol))
colnames(pmat_track) <- rownames(pmat_track) <- pc_df$Symbol
smat_track <- pmat_track


for (i in seq_along(mcg_p_ids)) {
  
  pmat <- fread_to_mat(mcg_p_ids[i], genes = pc_df$Symbol)
  smat <- fread_to_mat(mcg_s_ids[i], genes = pc_df$Symbol)
  
  # NAs to 0 so doesn't cause NA in summation
  pmat <- pmat %>% na_to_zero() %>% diag_to_one()
  smat <- smat %>% na_to_zero() %>% diag_to_one()

  # Integer tracking of whether the raw cor values were below 0
  pmat_neg <- (pmat < 0) * 1
  smat_neg <- (smat < 0) * 1
  
  pmat_track <- pmat_track + pmat_neg
  smat_track <- smat_track + smat_neg
  
}



fwrite(
  data.frame(pmat_track, check.names = FALSE),
  sep = "\t",
  row.names = TRUE,
  quote = FALSE,
  verbose = FALSE,
  showProgress = FALSE,
  file = file.path(amat_dir, "Microglia", "Mm_pcor", "count_negcor_matrix.tsv")
)


fwrite(
  data.frame(smat_track, check.names = FALSE),
  sep = "\t",
  row.names = TRUE,
  quote = FALSE,
  verbose = FALSE,
  showProgress = FALSE,
  file = file.path(amat_dir, "Microglia", "Mm_scor", "count_negcor_matrix.tsv")
)

