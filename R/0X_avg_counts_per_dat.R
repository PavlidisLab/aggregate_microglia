## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
library(pheatmap)
library(parallel)
library(cluster)
library(preprocessCore)
library(matrixStats)
source("R/00_config.R")
source("R/utils/functions.R")


mcg_meta <- read.delim(mcg_meta_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)

dat_l <- readRDS(mcg_dat_path)


mcg_meta <- distinct(mcg_meta, ID, .keep_all = TRUE)


# Assumes dat_l contains a gene x cell CPM matrix called "Mat"
# Assumes all matrices have the same genes and ordering
# Returns a list of the gene x experiment point estimates, as well as a summary
# df collapsing these summarizations across all datasets
# Note: matrixStats requires coercing sparse count matrices to dense.

summarize_count_list <- function(dat_l) {
  
  genes <- rownames(dat_l[[1]]$Mat)
  stopifnot(length(genes) > 0, identical(genes, rownames(dat_l[[length(dat_l)]]$Mat)))

  # Log transform all CPM matrices first
  count_l <- lapply(dat_l, function(x) as.matrix(log2(x$Mat + 1)))  
  
  # Get average, median, SD, and CV of log CPM counts and bind into a matrix
  avg <- do.call(cbind, lapply(count_l, rowMeans))
  med <- do.call(cbind, lapply(count_l, rowMedians))
  sd <- do.call(cbind, lapply(count_l, rowSds))
  cv <- sd / avg
  
  # Quantile normalize the averaged profiles
  qn_avg <- preprocessCore::normalize.quantiles(avg)
  
  # Rank and rank product of the averaged profiles
  rank_avg <- aggtools::colrank_mat(avg)
  rp_avg <- rowSums(log(rank_avg)) / length(genes)
  
  # Get binary gene measurement status (min 20 cells with at least one count)
  is_measured <- function(mat) rowSums(mat > 0) >= 20
  msr <- do.call(cbind, lapply(count_l, is_measured))
  
  # Summary dataset of point estimates for each gene collapsed across datasets

  summ_df <- data.frame(
    Symbol = genes,
    Avg = rowMeans(avg, na.rm = TRUE),
    Med = rowMedians(med, na.rm = TRUE),
    SD = rowMeans(sd, na.rm = TRUE),
    CV = rowMeans(cv, na.rm = TRUE),
    Avg_rank = rowMeans(rank_avg),
    RP = rank(rp_avg),
    N_msr = rowSums(msr)
  )

  return(list(
    Avg = avg,
    QN_Avg = qn_avg,
    Med = med,
    SD = sd,
    CV = cv,
    Msr = msr,
    Summ_df = summ_df
  ))
}



if (!file.exists(count_summ_path)) {
  count_summ_hg <- summarize_count_list(dat_l[ids_hg])
  count_summ_mm <- summarize_count_list(dat_l[ids_mm])
  saveRDS(list(Human = count_summ_hg, Mouse = count_summ_mm), count_summ_path)
  write.table(count_summ$Human$Summ_df, sep = "\t", quote = FALSE, row.names = FALSE, file = count_summ_table_hg)
  write.table(count_summ$Mouse$Summ_df, sep = "\t", quote = FALSE, row.names = FALSE, file = count_summ_table_mm)
} else {
  count_summ <- readRDS(count_summ_path)
}







# Trying to figure out genes with acceptable enough expression to keep
# ------------------------------------------------------------------------------


cutoff_hg <- 0.3
hist(log10(count_summ$Summ_hg$Avg + 1), breaks = 1000)
abline(v = cutoff_hg, col = "red")

# keep_hg <- filter(count_summ$Summ_hg, N_msr > 0 & Avg > 0 & Med > 0) %>% arrange(desc(Avg))
keep_hg <- filter(count_summ$Summ_hg, log10(Avg + 1) > cutoff_hg) %>% arrange(desc(Avg))
# keep_hg <- filter(count_summ$Summ_hg, N_msr > 2) %>% arrange(desc(N_msr))
# keep_hg <- filter(count_summ$Summ_hg, N_msr > 0 & Avg > 0) %>% arrange(desc(Avg))

hist(log10(keep_hg$Avg + 1), breaks = 1000)



cutoff_mm <- 0.5
hist(log10(count_summ$Summ_mm$Avg + 1), breaks = 1000)
abline(v = cutoff_mm, col = "red")

# keep_mm <- filter(count_summ$Summ_mm, N_msr > 0 & Avg > 0 & Med > 0) %>% arrange(desc(Avg))
keep_mm <- filter(count_summ$Summ_mm, log10(Avg + 1) > cutoff_mm) %>% arrange(desc(Avg))
# keep_mm <- filter(count_summ$Summ_mm, N_msr > 5) %>% arrange(desc(N_msr))
# keep_mm <- filter(count_summ$Summ_mm, N_msr > 0 & Avg > 0) %>% arrange(desc(Avg))

hist(log10(keep_mm$Avg + 1), breaks = 1000)





# Experiment similarity/clustering
# ------------------------------------------------------------------------------


cor_avg_hg <- cor(count_summ$Avg_hg[keep_hg$Symbol, ], method = "spearman")
cor_avg_mm <- cor(count_summ$Avg_mm[keep_mm$Symbol, ], method = "spearman")
# cor_avg_hg <- cor(avg_log_hg[keep_hg$Symbol, ], method = "spearman")
# cor_avg_mm <- cor(avg_log_mm[keep_mm$Symbol, ], method = "spearman")

cor_cv_hg <- cor(count_summ$CV_hg[keep_hg$Symbol, ], method = "spearman", use = "pairwise.complete.obs")
cor_cv_mm <- cor(count_summ$CV_mm[keep_mm$Symbol, ], method = "spearman", use = "pairwise.complete.obs")




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



hclust_avg_hg <- hclust(d = as.dist(1 - cor_avg_hg))
hclust_avg_mm <- hclust(d = as.dist(1 - cor_avg_mm))


hclust_cv_hg <- hclust(d = as.dist(1 - cor_cv_hg))
hclust_cv_mm <- hclust(d = as.dist(1 - cor_cv_mm))


sil_avg_hg <- cluster::silhouette(cutree(hclust_avg_hg, k = 2), as.dist(1 - cor_avg_hg))
mean(sil_avg_hg[, "sil_width"])
plot(sil_avg_hg)


sil_avg_mm <- cluster::silhouette(cutree(hclust_avg_mm, k = 2), as.dist(1 - cor_avg_mm))
mean(sil_avg_mm[, "sil_width"])
plot(sil_avg_mm)




# TODO: you have to be more explicit about expected shape in since it calls t()
pca_and_var <- function(mat) {
  
  # Performs PCA with prcomp and returns list of the resulting 
  # object as well as the variance explained
  
  # prcomp expects samples as rows, features (genes) as columns so transpose
  pcmat <- prcomp(t(mat), scale = TRUE)
  
  # variance explained by the PCs
  prc_var <- pcmat$sdev ^ 2
  var_explained <- round(prc_var / sum(prc_var) * 100, 2)
  cumvar_explained <- cumsum(var_explained)/sum(var_explained)
  return(list(PC = pcmat, 
              Var_explained = var_explained, 
              Cumvar_explained = cumvar_explained))
}



pc_hg <- pca_and_var(count_summ$Avg_hg[keep_hg$Symbol, ])
pc_mm <- pca_and_var(count_summ$Avg_mm[keep_mm$Symbol, ])
# pc_hg <- pca_and_var(avg_log_hg[keep_hg$Symbol, ])
# pc_mm <- pca_and_var(avg_log_mm[keep_mm$Symbol, ])


pc_df_hg <- data.frame(pc_hg$PC$x[, 1:5]) %>% rownames_to_column(var = "ID")
pc_df_mm <- data.frame(pc_mm$PC$x[, 1:5]) %>% rownames_to_column(var = "ID")


clus_df_hg <- data.frame(
  Avg_cluster = cutree(hclust_avg_hg, k = 2),
  CV_cluster = cutree(hclust_cv_hg, k = 2)) %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(filter(mcg_meta, Species == "Human"), by = "ID") %>% 
  left_join(pc_df_hg, by = "ID") %>% 
  mutate(Platform = str_replace_all(Platform, " ", ""),
         Platform = ifelse(is.na(Platform), "Mixed", Platform),
         Platform = ifelse(sapply(str_split(Platform, ","), length) != 1, "Mixed", Platform),
         Is_10X = str_detect(Platform, "^10x.*"))



clus_df_mm <- data.frame(
  Avg_cluster = cutree(hclust_avg_mm, k = 2),
  CV_cluster = cutree(hclust_cv_mm, k = 2)) %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(filter(mcg_meta, Species == "Mouse"), by = "ID") %>% 
  left_join(pc_df_mm, by = "ID") %>% 
  mutate(Platform = str_replace_all(Platform, " ", ""),
         Platform = ifelse(is.na(Platform), "Mixed", Platform),
         Platform = ifelse(sapply(str_split(Platform, ","), length) != 1, "Mixed", Platform),
         Is_10X = str_detect(Platform, "^10x.*"))


plot(hclust_avg_hg)
plot(hclust_avg_mm)

plot(hclust_cv_hg)
plot(hclust_cv_mm)


table(clus_df_hg$Avg_cluster, clus_df_hg$CV_cluster)
table(clus_df_mm$Avg_cluster, clus_df_mm$CV_cluster)



# Count of cells by cluster
boxplot(log10(clus_df_hg$N_cells+1) ~ clus_df_hg$Avg_cluster)
wilcox.test(clus_df_hg$N_cells ~ clus_df_hg$Avg_cluster)
# kruskal.test(clus_df_hg$N_cells ~ clus_df_hg$Avg_cluster)

boxplot(log10(clus_df_mm$N_cells+1) ~ clus_df_mm$Avg_cluster)
wilcox.test(clus_df_mm$N_cells ~ clus_df_mm$Avg_cluster)
# kruskal.test(clus_df_mm$N_cells ~ clus_df_mm$Avg_cluster)


# Median UMI by cluster
boxplot(log10(clus_df_hg$Median_UMI+1) ~ clus_df_hg$Avg_cluster)
wilcox.test(clus_df_hg$Median_UMI+1 ~ clus_df_hg$Avg_cluster)
# kruskal.test(clus_df_hg$Median_UMI ~ clus_df_hg$Avg_cluster)

boxplot(log10(clus_df_mm$Median_UMI+1) ~ clus_df_mm$Avg_cluster)
wilcox.test(clus_df_mm$Median_UMI+1 ~ clus_df_mm$Avg_cluster)
# kruskal.test(clus_df_mm$Median_UMI ~ clus_df_mm$Avg_cluster)


# Gene msr by cluster
boxplot(log10(clus_df_hg$N_msr_postfilt+1) ~ clus_df_hg$Avg_cluster)
wilcox.test(clus_df_hg$N_msr_postfilt+1 ~ clus_df_hg$Avg_cluster)
# kruskal.test(clus_df_hg$N_msr_postfilt ~ clus_df_hg$Avg_cluster)

boxplot(log10(clus_df_mm$N_msr_postfilt+1) ~ clus_df_mm$Avg_cluster)
wilcox.test(clus_df_mm$N_msr_postfilt+1 ~ clus_df_mm$Avg_cluster)
# kruskal.test(clus_df_mm$N_msr_postfilt ~ clus_df_mm$Avg_cluster)


fisher.test(table(clus_df_hg$Platform, clus_df_hg$Avg_cluster))
fisher.test(table(clus_df_mm$Platform, clus_df_mm$Avg_cluster))

fisher.test(table(clus_df_hg$Is_10X, clus_df_hg$Avg_cluster))
fisher.test(table(clus_df_mm$Is_10X, clus_df_mm$Avg_cluster))


round(cor(select_if(clus_df_hg, is.numeric)), 4)
round(cor(select_if(clus_df_mm, is.numeric)), 4)


boxplot(clus_df_hg$PC1 ~ clus_df_hg$Platform)
boxplot(clus_df_hg$PC2 ~ clus_df_hg$Platform)
kruskal.test(clus_df_hg$PC1 ~ clus_df_hg$Platform)
kruskal.test(clus_df_hg$PC2 ~ clus_df_hg$Platform)

boxplot(clus_df_mm$PC1 ~ clus_df_mm$Platform)
boxplot(clus_df_mm$PC2 ~ clus_df_mm$Platform)
kruskal.test(clus_df_mm$PC1 ~ clus_df_mm$Platform)
kruskal.test(clus_df_mm$PC2 ~ clus_df_mm$Platform)


boxplot(clus_df_hg$PC1 ~ clus_df_hg$Is_10X)
boxplot(clus_df_hg$PC2 ~ clus_df_hg$Is_10X)
wilcox.test(clus_df_hg$PC1 ~ clus_df_hg$Is_10X)
wilcox.test(clus_df_hg$PC2 ~ clus_df_hg$Is_10X)

boxplot(clus_df_mm$PC1 ~ clus_df_mm$Is_10X)
boxplot(clus_df_mm$PC2 ~ clus_df_mm$Is_10X)
kruskal.test(clus_df_mm$PC2 ~ clus_df_mm$Is_10X)
kruskal.test(clus_df_mm$PC1 ~ clus_df_mm$Is_10X)



ggplot(clus_df_hg, aes(x = PC1, y = PC2, fill = log10(N_cells))) +
  geom_point(size = 4, shape = 21) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  theme_classic()


ggplot(clus_df_mm, aes(x = PC1, y = PC2, fill = log10(N_cells))) +
  geom_point(size = 4, shape = 21) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  theme_classic()



ggplot(clus_df_hg, aes(x = PC1, y = PC2, fill = log10(Median_UMI))) +
  geom_point(size = 4, shape = 21) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  theme_classic()
plot(clus_df_hg$PC1, log10(clus_df_hg$Median_UMI+1))

ggplot(clus_df_mm, aes(x = PC1, y = PC2, fill = log10(Median_UMI))) +
  geom_point(size = 4, shape = 21) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  theme_classic()
plot(clus_df_mm$PC1, log10(clus_df_mm$Median_UMI+1))



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
