## Experiment similarity/clustering
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")

# Dataset meta and data IDs
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)


# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# Gene by dataset average log2 CPM matrices
avg_hg <- count_summ$Human$QN_Avg[count_summ$Human$Filter_genes, ]
avg_mm <- count_summ$Mouse$QN_Avg[count_summ$Mouse$Filter_genes, ]


# Exp-wise similarity
cor_avg_hg <- cor(avg_hg, method = "spearman")
cor_avg_mm <- cor(avg_mm, method = "spearman")


# view(mat_to_df(cor_avg_hg, value = "Cor"))
# view(mat_to_df(cor_avg_mm, value = "Cor"))


# Assumes mat is a dataset-dataset cor matrix

cor_heatmap <- function(mat) {
  
  min_break <- min(mat)
  
  pheatmap(mat,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = as.dist(1 - mat),
           clustering_distance_cols = as.dist(1 - mat),
           color = viridis::inferno(10),
           breaks = seq(min(min_break), 1, length.out = 10),
           border_color = "black",
           fontsize = 15)
}



hm_hg <- cor_heatmap(cor_avg_hg)
hm_mm <- cor_heatmap(cor_avg_mm)


cut_hg <- cutree(hm_hg$tree_row, k = 2)
cut_mm <- cutree(hm_mm$tree_col, k = 2)

png(height = 10, width = 12, units = "in", res = 300,
    file = file.path(plot_dir, "mcg_data_cluster_mm.png"))
hm_mm
dev.off()


png(height = 10, width = 12, units = "in", res = 300,
    file = file.path(plot_dir, "mcg_data_cluster_hg.png"))
hm_hg
dev.off()

# hclust_avg_hg <- hclust(d = as.dist(1 - cor_avg_hg))
# hclust_avg_mm <- hclust(d = as.dist(1 - cor_avg_mm))

# sil_avg_hg <- cluster::silhouette(cut_hg, as.dist(1 - cor_avg_hg))
# mean(sil_avg_hg[, "sil_width"])
# plot(sil_avg_hg)
# 
# 
# sil_avg_mm <- cluster::silhouette(cut_mm, as.dist(1 - cor_avg_mm))
# mean(sil_avg_mm[, "sil_width"])
# plot(sil_avg_mm)




# TODO: you have to be more explicit about expected shape in since it calls t()
# Performs PCA with prcomp and returns list of the resulting 
# object as well as the variance explained

pca_and_var <- function(mat) {
  
  # prcomp expects samples as rows, features (genes) as columns so transpose
  pcmat <- prcomp(t(mat), scale = TRUE)
  
  # variance explained by the PCs
  prc_var <- pcmat$sdev ^ 2
  var_explained <- round(prc_var / sum(prc_var) * 100, 2)
  cumvar_explained <- cumsum(var_explained) / sum(var_explained)
  
  return(list(PC = pcmat, 
              Var_explained = var_explained, 
              Cumvar_explained = cumvar_explained))
}



pca_hg <- pca_and_var(avg_hg)
pca_mm <- pca_and_var(avg_mm)



pca_df_hg <- data.frame(pca_hg$PC$x[, 1:5]) %>% rownames_to_column(var = "ID")
pca_df_mm <- data.frame(pca_mm$PC$x[, 1:5]) %>% rownames_to_column(var = "ID")


clus_df_hg <- data.frame(
  Avg_cut = cut_hg,
  ID = names(cut_hg)
  ) %>% 
  left_join(filter(mcg_meta, Species == "Human"), by = "ID") %>% 
  left_join(pca_df_hg, by = "ID") %>% 
  mutate(Platform = str_replace_all(Platform, " ", ""),
         Platform = ifelse(is.na(Platform), "Mixed", Platform),
         Platform = ifelse(sapply(str_split(Platform, ","), length) != 1, "Mixed", Platform),
         Is_10X = str_detect(Platform, "^10x.*"),
         Avg_cut = factor(Avg_cut, levels = unique(Avg_cut)))



clus_df_mm <- data.frame(
  Avg_cut = cut_mm,
  ID = names(cut_mm)
) %>% 
  left_join(filter(mcg_meta, Species == "Mouse"), by = "ID") %>% 
  left_join(pca_df_mm, by = "ID") %>% 
  mutate(Platform = str_replace_all(Platform, " ", ""),
         Platform = ifelse(is.na(Platform), "Mixed", Platform),
         Platform = ifelse(sapply(str_split(Platform, ","), length) != 1, "Mixed", Platform),
         Is_10X = str_detect(Platform, "^10x.*"),
         Avg_cut = factor(Avg_cut, levels = unique(Avg_cut)))





plot(hm_hg$tree_row)


plot(hm_mm$tree_row)



ggplot(clus_df_hg, aes(x = Avg_cut, y = log10(Median_UMI))) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, colour = "slategrey", size = 3.4, width = 0.05) +
  xlab("Cluster") +
  ggtitle("Human") +
  theme_classic() +
  theme(text = element_text(size = 20))


ggplot(clus_df_mm, aes(x = Avg_cut, y = log10(Median_UMI))) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, colour = "slategrey", size = 3.4, width = 0.05) +
  xlab("Cluster") +
  ggtitle("Mouse") +
  theme_classic() +
  theme(text = element_text(size = 20))


# Count of cells by cluster
boxplot(log10(clus_df_hg$N_cells+1) ~ clus_df_hg$Avg_cut)
wilcox.test(clus_df_hg$N_cells ~ clus_df_hg$Avg_cut)
# kruskal.test(clus_df_hg$N_cells ~ clus_df_hg$Avg_cut)

boxplot(log10(clus_df_mm$N_cells+1) ~ clus_df_mm$Avg_cut)
wilcox.test(clus_df_mm$N_cells ~ clus_df_mm$Avg_cut)
# kruskal.test(clus_df_mm$N_cells ~ clus_df_mm$Avg_cut)


# Median UMI by cluster
boxplot(log10(clus_df_hg$Median_UMI+1) ~ clus_df_hg$Avg_cut)
wilcox.test(clus_df_hg$Median_UMI+1 ~ clus_df_hg$Avg_cut)
# kruskal.test(clus_df_hg$Median_UMI ~ clus_df_hg$Avg_cut)

boxplot(log10(clus_df_mm$Median_UMI+1) ~ clus_df_mm$Avg_cut)
wilcox.test(clus_df_mm$Median_UMI+1 ~ clus_df_mm$Avg_cut)
# kruskal.test(clus_df_mm$Median_UMI ~ clus_df_mm$Avg_cut)


# Gene msr by cluster
boxplot(log10(clus_df_hg$N_msr_postfilt+1) ~ clus_df_hg$Avg_cut)
wilcox.test(clus_df_hg$N_msr_postfilt+1 ~ clus_df_hg$Avg_cut)
# kruskal.test(clus_df_hg$N_msr_postfilt ~ clus_df_hg$Avg_cut)

boxplot(log10(clus_df_mm$N_msr_postfilt+1) ~ clus_df_mm$Avg_cut)
wilcox.test(clus_df_mm$N_msr_postfilt+1 ~ clus_df_mm$Avg_cut)
# kruskal.test(clus_df_mm$N_msr_postfilt ~ clus_df_mm$Avg_cut)


fisher.test(table(clus_df_hg$Platform, clus_df_hg$Avg_cut))
fisher.test(table(clus_df_mm$Platform, clus_df_mm$Avg_cut))

fisher.test(table(clus_df_hg$Is_10X, clus_df_hg$Avg_cut))
fisher.test(table(clus_df_mm$Is_10X, clus_df_mm$Avg_cut))


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

