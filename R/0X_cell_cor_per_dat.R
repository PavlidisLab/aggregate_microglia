## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

mcg_meta <- read.delim(mcg_meta_path)
dat_l <- readRDS(mcg_dat_path)


mcg_meta <- mcg_meta %>% 
  group_by(ID) %>% 
  summarise(N_cells = sum(N_cells)) %>% 
  arrange(N_cells)


# Note that I inspected taking cor after filtering away 0 genes, but saw
# trivial change to summary stats

# TODO: check cell count and subsample at threshold

calc_cell_cor <- function(dat_l) {
  
  cor_l <- lapply(names(dat_l), function(x) {
    
    print(message(paste(x, Sys.time())))
    
    result <- tryCatch({
      cell_cor <- aggtools::sparse_pcor(dat_l[[x]]$Mat)
      cell_cor <- mat_to_df(cell_cor, value_name = "Cor")
      summary(cell_cor$Cor)
    }, 
    error = function(e) NULL)
    
    return(result)
    
  })
  names(cor_l) <- names(dat_l)
  
  return(cor_l)
}



calc_sample_cell_cor <- function(dat_l, 
                                 # mcg_meta, 
                                 n_cell_threshold = 10e3, 
                                 n_sample_cells = 1e3,
                                 n_samples = 10) {
  
  cor_l <- lapply(names(dat_l), function(x) {
    
    mat <- dat_l[[x]]$Mat
    
    sample_l <- lapply(1:n_samples, function(i) {
      sample_cells <-  sample(colnames(mat), n_sample_cells, replace = FALSE)
      submat <- mat[, sample_cells]
      cell_cor <- aggtools::sparse_pcor(submat)
      cell_cor <- mat_to_df(cell_cor, value_name = "Cor")[["Cor"]]
    })
    
    summary(unlist(sample_l))
    
    names(sample_l) <- names(dat_l)
    
    return(cor_l)
  })
}





if (!file.exists(cell_cor_path)) {
  cor_l <- calc_cell_cor(dat_l)
  saveRDS(cor_l, cell_cor_path)
} else {
  cor_l <- readRDS(cell_cor_path)
}



# TODO: some elements failed in lapply call, but can be regen isolated...
# figure out how to do one pass

null_cor <- mcg_meta %>% 
  filter(ID %in% names(which(sapply(cor_l, function(x) is.null(x))))) %>% 
  arrange(N_cells)


rerun_cor_ids <- filter(null_cor, N_cells < 20e3) %>% pull(ID) %>% unique()
sample_cor_ids <- setdiff(null_cor$ID, rerun_cor_ids)


rerun_cor <- calc_cell_cor(dat_l[rerun_cor_ids])



subsample_cor <- calc_sample_cell_cor(dat_l[sample_cor_ids])



test_id <- "GSE160523"
test_cell_cor <- aggtools::sparse_pcor(dat_l[[id]]$Mat)
test_cell_cor <- mat_to_df(cell_cor, value_name = "Cor")
summary(test_cell_cor$Cor)
plot(density(test_cell_cor$Cor))
# 
# keep_genes <- names(which(rowMeans(dat_l[[id]]$Mat) != 0))
# cell_cor_keep <- aggtools::sparse_pcor(dat_l[[id]]$Mat[keep_genes, ])
# cell_cor_keep <- mat_to_df(cell_cor_keep, value_name = "Cor")
# summary(cell_cor_keep$Cor)
