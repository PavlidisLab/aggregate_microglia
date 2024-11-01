# Experiment similarity/clustering
# ------------------------------------------------------------------------------


cor_avg_hg <- cor(count_summ$Human$QN_Avg[keep_hg$Symbol, ], method = "spearman")
cor_avg_mm <- cor(count_summ$Mouse$QN_Avg[keep_mm$Symbol, ], method = "spearman")

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