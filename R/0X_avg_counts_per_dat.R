## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
# library(aggtools)
library(pheatmap)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")


mcg_meta <- read.delim(mcg_meta_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)


dat_l <- readRDS(mcg_dat_path)


avg_hg <- do.call(cbind, lapply(dat_l[ids_hg], function(x) rowMeans(x$Mat)))
avg_mm <- do.call(cbind, lapply(dat_l[ids_mm], function(x) rowMeans(x$Mat)))


sd_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) apply(x$Mat, 1, sd), mc.cores = ncore))
sd_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) apply(x$Mat, 1, sd), mc.cores = ncore))


cv_hg <- sd_hg / avg_hg
cv_mm <- sd_mm / avg_mm




# TODO: compare cor using SD/CV, and cor after remove low genes

dat_cor_hg <- cor(avg_hg, method = "spearman")
dat_cor_mm <- cor(avg_mm, method = "spearman")


cor_heatmap <- function(mat) {
  
  min_break <-  min(mat)
  
  pheatmap(mat,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = viridis::inferno(10),
           breaks = seq(min(min_break), 1, length.out = 10),
           border_color = "black",
           fontsize = 20)
  
}


cor_heatmap(dat_cor_hg)
cor_heatmap(dat_cor_mm)


