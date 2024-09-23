## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
library(pheatmap)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")


mcg_meta <- read.delim(mcg_meta_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)


dat_l <- readRDS(mcg_dat_path)


# TODO: summarization into a function
# TODO: add more rank-based summarization; weigh against scaling counts
# TODO: is parallel helping or slowing


avg_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) rowMeans(x$Mat), mc.cores = ncore))
avg_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) rowMeans(x$Mat), mc.cores = ncore))

med_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) apply(x$Mat, 1, median), mc.cores = ncore))
med_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) apply(x$Mat, 1, median), mc.cores = ncore))

sd_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) apply(x$Mat, 1, sd), mc.cores = ncore))
sd_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) apply(x$Mat, 1, sd), mc.cores = ncore))

cv_hg <- sd_hg / avg_hg
cv_mm <- sd_mm / avg_mm

msr_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) rowSums(zero_sparse_cols(x$Mat)) == 0, mc.cores = ncore))
msr_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) rowSums(zero_sparse_cols(x$Mat)) == 0, mc.cores = ncore))




summ_hg <- data.frame(
  Symbol = rownames(avg_hg),
  Avg = rowMeans(avg_hg),
  Med = rowMeans(med_hg),
  SD = rowMeans(sd_hg),
  CV = rowMeans(cv_hg),
  N_msr = rowSums(msr_hg)
)



summ_mm <- data.frame(
  Symbol = rownames(avg_mm),
  Avg = rowMeans(avg_mm),
  Med = rowMeans(med_mm),
  SD = rowMeans(sd_mm),
  CV = rowMeans(cv_mm),
  N_msr = rowSums(msr_mm)
)




all0_hg <- filter(Symbol)





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





# More deviation between avg and median at lower end
plot(summ_hg$Avg, summ_hg$Med)
plot(summ_hg$Avg[summ_hg$Avg < 1000 & summ_hg$Med], summ_hg$Med[summ_hg$Avg < 1000 & summ_hg$Med])

plot(summ_hg$Avg, summ_hg$SD)
plot(summ_hg$Med, summ_hg$CV)
