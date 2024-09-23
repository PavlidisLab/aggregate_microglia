## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

mcg_meta <- read.delim(mcg_meta_path)
dat_l <- readRDS(mcg_dat_path)


# Note that I inspected taking cor after filtering away 0 genes, but saw
# trivial change to summary stats

calc_cell_cor <- function(dat_l) {
  
  cor_l <- lapply(dat_l, function(x) {
    
    message(paste(x, Sys.time()))
    
    cell_cor <- aggtools::sparse_pcor(x$Mat)
    cell_cor <- mat_to_df(cell_cor, value_name = "Cor")
    summary(cell_cor$Cor)
  })
  names(cor_l) <- names(dat_l)
  
  return(cor_l)
}


if (!file.exists(cell_cor_path)) {
  cell_cor <- calc_cell_cor(dat_l)
  saveRDS(cell_cor, cell_cor_path)
} else {
  cell_cor <- readRDS(cell_cor_path)
}



print(length(cell_cor))



# id <- "GSE121891"
# cell_cor <- aggtools::sparse_pcor(dat_l[[id]]$Mat)
# cell_cor <- mat_to_df(cell_cor, value_name = "Cor")
# summary(cell_cor$Cor)
# plot(density(cell_cor$Cor))
# 
# keep_genes <- names(which(rowMeans(dat_l[[id]]$Mat) != 0))
# cell_cor_keep <- aggtools::sparse_pcor(dat_l[[id]]$Mat[keep_genes, ])
# cell_cor_keep <- mat_to_df(cell_cor_keep, value_name = "Cor")
# summary(cell_cor_keep$Cor)
