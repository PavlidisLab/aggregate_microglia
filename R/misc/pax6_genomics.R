

bind_list <- readRDS(bind_summary_path)
rank_list <- readRDS(rank_path)
dat_list <- readRDS(alldat_path)
de_list <- readRDS(prank_path_group)
"/space/scratch/amorin/R_objects//bind_summary_refseq_Apr2022.RDS"
"/space/scratch/amorin/R_objects//ranked_target_list_Apr2022.RDS"
"/space/scratch/amorin/R_objects//all_data_list_Apr2022.RDS"
"/space/grp/amorin/Expression_files/Gemma/TF_perturb_DE_counts_list_by_TF_FDR01_Apr2022.RDS"



library(tidyverse)
library(pheatmap)
library(viridisLite)
source("R/setup-01_config.R")
source("R/utils/plot_functions.R")

bind_list <- readRDS(bind_summary_path)
rank_list <- readRDS(rank_path)
dat_list <- readRDS(alldat_path)
de_list <- readRDS(prank_path_group)

df <- bind_list$Human_bind_tf$ASCL1 %>% 
  mutate(Top = Symbol %in% slice_max(., Mean_bind, n = 50)$Symbol)

pa <- 
  ggplot() +
  geom_point(data = df, 
             aes(y = Proportion_binary, x = Mean_bind),
             shape = 21, size = 2.3, alpha = 0.4, colour = "black") +
  geom_text(data = filter(df, Top),
            check_overlap = TRUE,
            aes(y = Proportion_binary, x = Mean_bind, label = Symbol)) +
  ylab("Proportion of bound experiments (+/- 25kb)") +
  xlab("Mean binding score") +
  ggtitle("ASCL1") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "none")


ascl1_mat <- dat_list$Binding$Human$QN_log[, filter(dat_list$Binding$Meta, Symbol == "ASCL1")$Experiment_ID]

df <- filter(df, Symbol %in% rownames(ascl1_mat)) %>% arrange(desc(Mean_bind))

top_bound <- slice_max(df, Mean_bind, n = 5)$Symbol

mid_bound <- sample(df[1000:2000,]$Symbol, 5)

bottom_bound <- sample(slice_min(df, Mean_bind)$Symbol, 5)

plot_mat <- ascl1_mat[c(top_bound, mid_bound, bottom_bound), ]

plot_mat <- cbind(plot_mat, Mean = rowMeans(plot_mat))


# pal <- rev(magma(n = 11))
# pal <- rev(viridis(n = 11))
# pal <- c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b')
pal <- c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b')
color_breaks <- seq(0, 0.7, length.out = length(pal))



png(width = 9, height = 6, units = "in", res = 300,
    filename = "~/Plots/Chipseq/Binding_summary/ASCL1_bind_heatmap.png")

pheatmap(
  plot_mat,
  color = pal,
  breaks = color_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  border_col = "black",
  cellheight = 20,
  cellwidth = 20,
  fontsize_row = 20,
  gaps_col = ncol(plot_mat) - 1,
  gaps_row = c(5, 10),
  height = 6,
  width = 9
)

graphics.off()



##


perturb_df <- de_list$Human$ASCL1 %>% 
  filter(Count_NA == 0) %>% 
  arrange(desc(Count_DE), desc(Avg_abs_FC))

perturb_top <- head(perturb_df, 5)

perturb_mid  <- perturb_df %>% filter(Count_DE == 2) %>% slice_sample(n = 5)

perturb_low <- tail(perturb_df, 5)

perturb_plot <- do.call(rbind, list(perturb_top, perturb_mid, perturb_low))


fc_heatmap <- function(tf,
                       tf_df,
                       fc_mat,
                       meta,
                       FC_min = -2.5,
                       FC_max = 2.5,
                       legend_arg = TRUE,
                       anno_legend = TRUE) {
  
  
  symbols <- tf_df$Symbol
  tf_meta <- filter(meta, Symbol == tf) %>% arrange(Perturbation)
  fc_mat <- fc_mat[symbols, tf_meta$Experiment_ID]
  
  anno <- tf_meta[, "Perturbation", drop = FALSE]
  anno$Perturbation <- factor(anno$Perturbation, levels = unique(anno$Perturbation))
  droplevels(anno)
  rownames(anno) <- colnames(fc_mat)
  
  
  color_breaks <- seq(FC_min, FC_max, length.out = pal_length)
  
  pheatmap(
    fc_mat,
    color = bluered_pal,
    breaks = color_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    border_col = "black",
    annotation_col = anno,
    annotation_colors = pert_anno,
    annotation_names_col = FALSE,
    annotation_legend = anno_legend,
    cellheight = 20,
    cellwidth = 20,
    fontsize_row = 20,
    legend = legend_arg,
    gaps_row = c(5, 10),
    height = 6,
    width = 9
  )
}


png(width = 6, height = 6, units = "in", res = 300,
    filename = "~/Plots/TF_perturb/Describe_FDR_counts/segmented_ASCL1_FC_heatmap.png")

fc_heatmap(tf = "ASCL1", tf_df = perturb_plot, fc_mat = dat_list$Perturbation$Human$FC_mat, meta = dat_list$Perturbation$Meta)

graphics.off()
