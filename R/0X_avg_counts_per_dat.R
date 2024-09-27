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


mcg_meta <- distinct(mcg_meta, ID, .keep_all = TRUE)


# TODO: summarization into a function
# TODO: add more rank-based summarization; weigh against scaling counts
# TODO: is parallel helping or slowing


calc_count_summaries <- function(dat_l, ids_hg, ids_mm) {
  
  
  avg_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) rowMeans(x$Mat), mc.cores = ncore))
  avg_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) rowMeans(x$Mat), mc.cores = ncore))
  
  med_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) apply(x$Mat, 1, median), mc.cores = ncore))
  med_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) apply(x$Mat, 1, median), mc.cores = ncore))
  
  sd_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) apply(x$Mat, 1, sd), mc.cores = ncore))
  sd_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) apply(x$Mat, 1, sd), mc.cores = ncore))
  
  cv_hg <- sd_hg / avg_hg
  cv_mm <- sd_mm / avg_mm
  
  msr_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) rowSums(zero_sparse_cols(x$Mat)) != 0, mc.cores = ncore))
  msr_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) rowSums(zero_sparse_cols(x$Mat)) != 0, mc.cores = ncore))
  
  
  summ_hg <- data.frame(
    Symbol = rownames(avg_hg),
    Avg = rowMeans(avg_hg, na.rm = TRUE),
    Med = rowMeans(med_hg, na.rm = TRUE),
    SD = rowMeans(sd_hg, na.rm = TRUE),
    CV = rowMeans(cv_hg, na.rm = TRUE),
    N_msr = rowSums(msr_hg)
  )
  
  
  summ_mm <- data.frame(
    Symbol = rownames(avg_mm),
    Avg = rowMeans(avg_mm, na.rm = TRUE),
    Med = rowMeans(med_mm, na.rm = TRUE),
    SD = rowMeans(sd_mm, na.rm = TRUE),
    CV = rowMeans(cv_mm, na.rm = TRUE),
    N_msr = rowSums(msr_mm)
  )
  
  list(Avg_hg = avg_hg,
       Avg_mm = avg_mm,
       Med_hg = med_hg,
       Med_mm = med_mm,
       SD_hg = sd_hg,
       SD_mm = sd_mm,
       CV_hg = cv_hg,
       CV_mm = cv_mm,
       Msr_hg = msr_hg,
       Msr_mm = msr_mm,
       Summ_hg = summ_hg,
       Summ_mm = summ_mm)
  
  
}



if (!file.exists(count_summ_path)) {
  count_summ <- calc_count_summaries(dat_l, ids_hg, ids_mm)
  saveRDS(count_summ, count_summ_path)
} else {
  count_summ <- readRDS(count_summ_path)
}


rank_avg_hg <- aggtools::colrank_mat(count_summ$Avg_hg)
rank_avg_mm <- aggtools::colrank_mat(count_summ$Avg_mm)

avg_log_hg <- do.call(cbind, mclapply(dat_l[ids_hg], function(x) rowMeans(log2(x$Mat+1)), mc.cores = ncore))
avg_log_mm <- do.call(cbind, mclapply(dat_l[ids_mm], function(x) rowMeans(log2(x$Mat+1)), mc.cores = ncore))


count_summ$Summ_hg <- cbind(
  count_summ$Summ_hg,
  Avg_log_avg = rowMeans(avg_log_hg),
  Avg_rank_avg = rowMeans(rank_avg_hg),
  RP_avg = rank((rowSums(log(rank_avg_hg)) / nrow(rank_avg_hg)))
)


count_summ$Summ_mm <- cbind(
  count_summ$Summ_mm,
  Avg_log_avg = rowMeans(avg_log_mm),
  Avg_rank_avg = rowMeans(rank_avg_mm),
  RP_avg = rank((rowSums(log(rank_avg_mm)) / nrow(rank_avg_mm)))
)








# Trying to figure out genes with acceptable enough expression to keep
# ------------------------------------------------------------------------------

plot(density(log10(count_summ$Summ_hg$Avg + 1)))
hist(log10(count_summ$Summ_hg$Avg + 1), breaks = 1000)
cutoff_hg <- 0.3
abline(v = cutoff_hg, col = "red")

# keep_hg <- filter(count_summ$Summ_hg, N_msr > 0 & Avg > 0 & Med > 0) %>% arrange(desc(Avg))
keep_hg <- filter(count_summ$Summ_hg, log10(Avg + 1) > cutoff_hg) %>% arrange(desc(Avg))
# keep_hg <- filter(count_summ$Summ_hg, N_msr > 2) %>% arrange(desc(N_msr))
# keep_hg <- filter(count_summ$Summ_hg, N_msr > 0 & Avg > 0) %>% arrange(desc(Avg))


plot(density(log10(keep_hg$Avg + 1)))
hist(log10(keep_hg$Avg + 1), breaks = 1000)




plot(density(log10(count_summ$Summ_mm$Avg + 1)))
hist(log10(count_summ$Summ_mm$Avg + 1), breaks = 1000)
cutoff_mm <- 0.5
abline(v = cutoff_mm, col = "red")

# keep_mm <- filter(count_summ$Summ_mm, N_msr > 0 & Avg > 0 & Med > 0) %>% arrange(desc(Avg))
keep_mm <- filter(count_summ$Summ_mm, log10(Avg + 1) > cutoff_mm) %>% arrange(desc(Avg))
# keep_mm <- filter(count_summ$Summ_mm, N_msr > 5) %>% arrange(desc(N_msr))
# keep_mm <- filter(count_summ$Summ_mm, N_msr > 0 & Avg > 0) %>% arrange(desc(Avg))


plot(density(log10(keep_mm$Avg + 1)))
hist(log10(keep_mm$Avg + 1), breaks = 1000)





# Experiment similarity/clustering
# ------------------------------------------------------------------------------



# TODO: compare cor using SD/CV, and cor after remove low genes

# cor_avg_hg <- cor(count_summ$Avg_hg, method = "spearman")
# cor_avg_mm <- cor(count_summ$Avg_mm, method = "spearman")

cor_avg_hg <- cor(count_summ$Avg_hg[keep_hg$Symbol, ], method = "spearman")
cor_avg_mm <- cor(count_summ$Avg_mm[keep_mm$Symbol, ], method = "spearman")



cor_cv_hg <- cor(count_summ$CV_hg, method = "spearman", use = "pairwise.complete.obs")
cor_cv_mm <- cor(count_summ$CV_mm, method = "spearman", use = "pairwise.complete.obs")




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


cor_heatmap(cor_avg_hg)
cor_heatmap(cor_avg_mm)


cor_heatmap(cor_cv_hg)
cor_heatmap(cor_cv_mm)



hclust_hg <- hclust(d = as.dist(1 - cor_avg_hg))
hclust_mm <- hclust(d = as.dist(1 - cor_avg_mm))
plot(hclust_hg)
plot(hclust_mm)


clus_df_hg <- data.frame(Group = cutree(hclust_hg, k = 2)) %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(filter(mcg_meta, Species == "Human"), by = "ID")



clus_df_mm <- data.frame(Group = cutree(hclust_mm, k = 2)) %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(filter(mcg_meta, Species == "Mouse"), by = "ID")




# Count of cells by cluster
boxplot(log10(clus_df_hg$N_cells+1) ~ clus_df_hg$Group)
wilcox.test(clus_df_hg$N_cells ~ clus_df_hg$Group)
# kruskal.test(clus_df_hg$N_cells ~ clus_df_hg$Group)

boxplot(log10(clus_df_mm$N_cells+1) ~ clus_df_mm$Group)
wilcox.test(clus_df_mm$N_cells ~ clus_df_mm$Group)
# kruskal.test(clus_df_mm$N_cells ~ clus_df_mm$Group)


# Median UMI by cluster
boxplot(log10(clus_df_hg$Median_UMI+1) ~ clus_df_hg$Group)
wilcox.test(clus_df_hg$Median_UMI+1 ~ clus_df_hg$Group)
# kruskal.test(clus_df_hg$N_cells ~ clus_df_hg$Group)

boxplot(log10(clus_df_mm$Median_UMI+1) ~ clus_df_mm$Group)
wilcox.test(clus_df_mm$Median_UMI+1 ~ clus_df_mm$Group)
# kruskal.test(clus_df_hg$N_cells ~ clus_df_hg$Group)



# Gene msr by cluster
boxplot(log10(clus_df_hg$N_msr_postfilt+1) ~ clus_df_hg$Group)
wilcox.test(clus_df_hg$N_msr_postfilt+1 ~ clus_df_hg$Group)
# kruskal.test(clus_df_hg$N_cells ~ clus_df_hg$Group)

boxplot(log10(clus_df_mm$N_msr_postfilt+1) ~ clus_df_mm$Group)
wilcox.test(clus_df_mm$N_msr_postfilt+1 ~ clus_df_mm$Group)
# kruskal.test(clus_df_hg$N_cells ~ clus_df_hg$Group)




# Plotting average expression of counts
# ------------------------------------------------------------------------------

# head(count_summ$Summ_hg)

plot_mat <- count_summ$Avg_hg[keep_hg$Symbol, ]
plot_mat <- cbind(plot_mat, Average = count_summ$Summ_hg[keep_hg$Symbol, "Avg"])
plot_mat <- scale(plot_mat)
plot_mat[plot_mat > 3] <- 3

pheatmap(plot_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = RColorBrewer::brewer.pal(9, "Blues"),
         show_rownames = FALSE,
         # breaks = seq(0, 1, length.out = 10),
         border_color = NA,
         fontsize = 20,
         gaps_col = ncol(plot_mat) - 1)



plot_mat <- count_summ$Avg_mm[keep_mm$Symbol, ]
plot_mat <- cbind(plot_mat, Average = count_summ$Summ_mm[keep_mm$Symbol, "Avg"])
plot_mat <- scale(plot_mat)
plot_mat[plot_mat > 3] <- 3

pheatmap(plot_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = RColorBrewer::brewer.pal(9, "Blues"),
         show_rownames = FALSE,
         # breaks = seq(0, 1, length.out = 10),
         border_color = NA,
         fontsize = 20,
         gaps_col = ncol(plot_mat) - 1)





# Plotting binary measurement for coexpression
# ------------------------------------------------------------------------------


plot_mat <- count_summ$Msr_hg[keep_hg$Symbol, ] * 1



pheatmap(plot_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         show_rownames = FALSE,
         # breaks = seq(0, 1, length.out = 10),
         border_color = NA,
         fontsize = 20)



plot_mat <- count_summ$Msr_mm[keep_mm$Symbol, ] * 1


pheatmap(plot_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         show_rownames = FALSE,
         # breaks = seq(0, 1, length.out = 10),
         border_color = NA,
         fontsize = 20)



# Plotting counts of a single experiments
# ------------------------------------------------------------------------------


plot_mat <- dat_l$GSE118020$Mat
plot_mat <- cbind(plot_mat, Average = count_summ$Avg_mm[, "GSE118020"])


pheatmap(plot_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = RColorBrewer::brewer.pal(9, "Reds"),
         show_rownames = FALSE,
         show_colnames = FALSE,
         # breaks = seq(0, 1, length.out = 10),
         border_color = NA,
         fontsize = 20,
         gaps_col = ncol(plot_mat) - 1)









# More deviation between avg and median at lower end
plot(summ_hg$Avg, summ_hg$Med)
plot(summ_hg$Avg[summ_hg$Avg < 1000 & summ_hg$Med], summ_hg$Med[summ_hg$Avg < 1000 & summ_hg$Med])

plot(summ_hg$Avg, summ_hg$SD)
plot(summ_hg$Med, summ_hg$CV)
