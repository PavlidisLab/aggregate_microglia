## Load single cell data and filter count matrix and meta to just microglia
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

mcg_meta <- read.delim(mcg_meta_path)


load_and_filter_dat <- function(input_df) {
  
  ids <- unique(input_df$ID)
  
  dat_l <- lapply(ids, function(id) {
    
    message(paste(id, Sys.time()))
    df <- filter(input_df, ID == id)
    path <- unique(df$Path)
    
    dat <- aggtools::load_scdat(path)
    meta <- filter(dat$Meta, Cell_type %in% df$Cell_type)
    mat <- dat$Mat[, meta$ID]
    list(Mat = mat, Meta = meta)
    
  })
  names(dat_l) <- ids
  
  return(dat_l)
}



# TODO: run if not saved logic
dat_l <- load_and_filter_dat(mcg_meta)

saveRDS(dat_l, mcg_dat_path)


dat <- aggtools::load_scdat(mcg_meta$Path[1])
meta <- filter(dat$Meta, Cell_type %in% mcg_meta$Cell_type[1])
mat <- dat$Mat[, meta$ID]


mat2 <- t(zero_sparse_cols(t(mat)))

count_before <- rowSums(mat)
count_after <- rowSums(mat2)

sum(count_before != 0)
sum(count_after != 0)


plot(count_before, count_after)
