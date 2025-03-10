## Each microglia dataset can be represented as a pseuodbulked vector. Exploring
## the similarity/clustering these datasets
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(cluster)
library(viridis)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")

# Dataset meta and data IDs
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)


# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)



# Functions
# ------------------------------------------------------------------------------


# Assumes mat is a dataset-dataset correlation matrix; produces a pheatmap. This
# is also used to extract the hclust object/tree

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



# Performs PCA with prcomp and returns list of the resulting object as well as 
# the variance explained. Assumes features as columns and samples as rows

pca_and_var <- function(mat) {
  
  pc_mat <- prcomp(mat, scale = TRUE)
  
  # variance explained by the PCs
  prc_var <- pc_mat$sdev ^ 2
  var_explained <- round(prc_var / sum(prc_var) * 100, 2)
  cumvar_explained <- cumsum(var_explained) / sum(var_explained)
  
  return(list(PC = pc_mat, 
              Var_explained = var_explained, 
              Cumvar_explained = cumvar_explained))
}




# Generating similarity between datasets
# ------------------------------------------------------------------------------



# Isolate gene by dataset average/pseuodbulked log2+quantile norm CPM matrices
mat_hg <- count_summ$Human$QN_Avg[count_summ$Human$Filter_genes, ]
mat_mm <- count_summ$Mouse$QN_Avg[count_summ$Mouse$Filter_genes, ]


# Experiment-wise spearman's correlation
cor_mat_hg <- cor(mat_hg, method = "spearman")
cor_mat_mm <- cor(mat_mm, method = "spearman")


# Correlation of experiment pairs as a table
cor_df_hg <- mat_to_df(cor_mat_hg, value_name = "Cor")
cor_df_mm <- mat_to_df(cor_mat_mm, value_name = "Cor")


# Heatmap object (includes hclust) of experiment correlation
heatmap_hg <- cor_heatmap(cor_mat_hg)
heatmap_mm <- cor_heatmap(cor_mat_mm)


# Cutting the dendogram produced by pheatmap. k=2 chosen by visual exam.
cut_hg <- cutree(heatmap_hg$tree_row, k = 2)
cut_mm <- cutree(heatmap_mm$tree_row, k = 2)


# This produces the same hclust as the above pheatmap object
hclust_mat_hg <- hclust(d = as.dist(1 - cor_mat_hg))
hclust_mat_mm <- hclust(d = as.dist(1 - cor_mat_mm))


# Silhoutte analysis aligns with visual inspection of k=2
sil_mat_hg <- cluster::silhouette(cut_hg, as.dist(1 - cor_mat_hg))
# mean(sil_mat_hg[, "sil_width"])
# plot(sil_mat_hg)

sil_mat_mm <- cluster::silhouette(cut_mm, as.dist(1 - cor_mat_mm))
# mean(sil_mat_mm[, "sil_width"])
# plot(sil_mat_mm)



# PCA: expects samples as rows, features (genes) as columns so must transpose
pca_hg <- pca_and_var(t(mat_hg))
pca_mm <- pca_and_var(t(mat_mm))


# Isolating first 5 PCs
pca_df_hg <- data.frame(pca_hg$PC$x[, 1:5]) %>% rownames_to_column(var = "ID")
pca_df_mm <- data.frame(pca_mm$PC$x[, 1:5]) %>% rownames_to_column(var = "ID")


# Organize hclust cuts and PCs into a df, and simplify platform information
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


# Again for mouse
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



# Correlation between numeric metadata factors and PCs
meta_cor_hg <- round(cor(select_if(clus_df_hg, is.numeric)), 4)
meta_cor_mm <- round(cor(select_if(clus_df_mm, is.numeric)), 4)



# Plots
# ------------------------------------------------------------------------------



# Saving out human correlation heatmap
png(height = 10, width = 12, units = "in", res = 300,
    file = file.path(plot_dir, "mcg_data_cluster_mm.png"))
heatmap_mm
dev.off()

# Saving out mouse correlation heatmap
png(height = 10, width = 12, units = "in", res = 300,
    file = file.path(plot_dir, "mcg_data_cluster_hg.png"))
heatmap_hg
dev.off()



# Closer examination of dendogram 
# plot(heatmap_hg$tree_row)
# plot(heatmap_mm$tree_row)


# The following was used to inspect potential technical factors driving 
# clustering. Median UMI seemed to have the strongest association


p1a <- ggplot(clus_df_hg, aes(x = Avg_cut, y = log10(Median_UMI))) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, colour = "slategrey", size = 3.4, width = 0.05) +
  xlab("Cluster") +
  ggtitle("Human") +
  theme_classic() +
  theme(text = element_text(size = 20))


p1b <- ggplot(clus_df_mm, aes(x = Avg_cut, y = log10(Median_UMI))) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, colour = "slategrey", size = 3.4, width = 0.05) +
  xlab("Cluster") +
  ggtitle("Mouse") +
  theme_classic() +
  theme(text = element_text(size = 20))


# wilcox.test(clus_df_hg$Median_UMI+1 ~ clus_df_hg$Avg_cut)
# wilcox.test(clus_df_mm$Median_UMI+1 ~ clus_df_mm$Avg_cut)


# Count of cells by cluster
# boxplot(log10(clus_df_hg$N_cells+1) ~ clus_df_hg$Avg_cut)
# boxplot(log10(clus_df_mm$N_cells+1) ~ clus_df_mm$Avg_cut)
# wilcox.test(clus_df_hg$N_cells ~ clus_df_hg$Avg_cut)
# wilcox.test(clus_df_mm$N_cells ~ clus_df_mm$Avg_cut)


# Gene msr by cluster
# boxplot(log10(clus_df_hg$N_msr_postfilt+1) ~ clus_df_hg$Avg_cut)
# boxplot(log10(clus_df_mm$N_msr_postfilt+1) ~ clus_df_mm$Avg_cut)
# wilcox.test(clus_df_hg$N_msr_postfilt+1 ~ clus_df_hg$Avg_cut)
# wilcox.test(clus_df_mm$N_msr_postfilt+1 ~ clus_df_mm$Avg_cut)


# Platform
# fisher.test(table(clus_df_hg$Platform, clus_df_hg$Avg_cut))
# fisher.test(table(clus_df_mm$Platform, clus_df_mm$Avg_cut))
# fisher.test(table(clus_df_hg$Is_10X, clus_df_hg$Avg_cut))
# fisher.test(table(clus_df_mm$Is_10X, clus_df_mm$Avg_cut))
