library(tidyverse)
library(pheatmap)
source("R/00_config.R")
source("R/utils/functions.R")


mcg_summ <- readRDS(mcg_count_summ_list_path)
macro_summ <- readRDS(macro_count_summ_list_path)

species <- "Human"


summ_df <- left_join(
  mcg_summ[[species]]$Summ_df,
  macro_summ[[species]]$Summ_df,
  by = "Symbol",
  suffix = c("_mcg", "_macro")
  ) %>% 
  select(Symbol, QN_avg_mcg, QN_avg_macro, N_msr_mcg, N_msr_macro) %>% 
  mutate(Diff = QN_avg_mcg - QN_avg_macro) %>% 
  relocate(Diff, .after = Symbol)


# TODO: should macro have the same n for filter?
keep_mcg <- filter(summ_df, N_msr_mcg >= (max(N_msr_mcg) * 1/3))
keep_macro <- filter(summ_df, N_msr_macro >= (max(N_msr_macro) * 1/3))
common <- intersect(keep_mcg$Symbol, keep_macro$Symbol)
diff_mcg <- setdiff(keep_mcg$Symbol, keep_macro$Symbol)
diff_macro <- setdiff(keep_macro$Symbol, keep_mcg$Symbol)
keep_union <- union(keep_mcg$Symbol, keep_macro$Symbol)


avg_mat <- cbind(mcg_summ[[species]]$QN_Avg, macro_summ[[species]]$QN_Avg)
avg_mat <- avg_mat[keep_union, ]




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



pca <- pca_and_var(avg_mat)



pc_df <- data.frame(
  ID = colnames(avg_mat),
  Cell_type = c(rep("Microglia", ncol(mcg_summ[[species]]$QN_Avg)),
                rep("Macrophage", ncol(macro_summ[[species]]$QN_Avg))),
  pca$PC$x[, 1:5]
)



summ_df %>% 
  filter(Symbol %in% keep_union) %>% 
  ggplot(., aes(x = QN_avg_mcg, y = QN_avg_macro)) +
  geom_point(shape = 21, size = 2.4, alpha = 0.6) +
  ggtitle(species) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



summ_df %>% 
  filter(Symbol %in% keep_union) %>% 
  ggplot(., aes(x = Diff)) +
  geom_histogram(bins = 100) +
  ggtitle(species) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


filter(summ_df, N_msr_macro == 0) %>% view


ggplot(pc_df, aes(x = PC1, y = PC2, fill = Cell_type)) +
  geom_point(size = 3, shape = 21, colour = "black") +
  scale_fill_manual(values = c("firebrick", "darkgrey")) +
  ggtitle(species) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = c(0.80, 0.95),
        plot.margin = margin(c(10, 20, 10, 10)))


boxplot(pc_df$PC1 ~ pc_df$Cell_type)
boxplot(pc_df$PC2 ~ pc_df$Cell_type)
boxplot(pc_df$PC3 ~ pc_df$Cell_type)
boxplot(pc_df$PC4 ~ pc_df$Cell_type)
boxplot(pc_df$PC5 ~ pc_df$Cell_type)



