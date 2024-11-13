## TODO: Keep genes as stored vector (instead of filtered df?)
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
library(pheatmap)
library(preprocessCore)
library(matrixStats)
source("R/00_config.R")
source("R/utils/functions.R")

# Dataset meta and human/mouse data IDs
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)

# List of count matrices and their metadata
dat_l <- readRDS(mcg_dat_path)



if (!file.exists(mcg_count_summ_list_path)) {
  count_summ_hg <- summarize_gene_counts(dat_l[ids_hg])
  count_summ_mm <- summarize_gene_counts(dat_l[ids_mm])
  saveRDS(list(Human = count_summ_hg, Mouse = count_summ_mm), mcg_count_summ_list_path)
  write.table(count_summ_hg$Summ_df, sep = "\t", quote = FALSE, row.names = FALSE, file = mcg_count_summ_table_hg)
  write.table(count_summ_mm$Summ_df, sep = "\t", quote = FALSE, row.names = FALSE, file = mcg_count_summ_table_mm)
} else {
  count_summ <- readRDS(mcg_count_summ_list_path)
}


stop()




# Identifying genes to filter
# NOTE: Median is wildly stringent (median of medians too harsh?)
# Heuristic: gene must be measured (min 20 cells) in at least a third of datasets
# ------------------------------------------------------------------------------

summ_df_hg <- count_summ$Human$Summ_df
cutoff_hg <- floor(length(ids_hg) * (1/3))
keep_hg <- filter(summ_df_hg, N_msr >= cutoff_hg) %>% arrange(desc(N_msr))


summ_df_mm <- count_summ$Mouse$Summ_df
cutoff_mm <- floor(length(ids_mm) * (1/3))
keep_mm <- filter(summ_df_mm, N_msr >= cutoff_mm) %>% arrange(desc(Avg))






# Plots
# ------------------------------------------------------------------------------


# Hists before/after filter 
hist(summ_df_hg$Avg, breaks = 1000)
hist(keep_hg$Avg, breaks = 1000)

hist(summ_df_mm$Avg, breaks = 1000)
hist(keep_mm$Avg, breaks = 1000)



# Vbplot of counts for given gene


vbplot <- function(dat_l, gene) {
  
  ids <- names(dat_l)
  
  gene_l <- lapply(ids, function(x) {
    data.frame(Log2_count = log2(dat_l[[x]]$Mat[gene, ] + 1), ID = x)
  })
  
  gene_df <- do.call(rbind, gene_l)
  
  ggplot(gene_df, aes(x = ID, y = Log2_count)) +
    geom_violin(fill = "slategrey") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    ylab("Log[[2]] CPM") +
    ggtitle(gene) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}



# The max average gene
vbplot(dat_l[ids_hg], slice_max(keep_hg, QN_avg, n = 1)$Symbol)
vbplot(dat_l[ids_mm], slice_max(keep_mm, QN_avg, n = 1)$Symbol)

# The lowest average gene that was kept
vbplot(dat_l[ids_hg], slice_min(keep_hg, QN_avg, n = 1)$Symbol)
vbplot(dat_l[ids_mm], slice_min(keep_mm, QN_avg, n = 1)$Symbol)

# The highest average gene that was removed
gene_hg <- filter(summ_df_hg, Symbol %!in% keep_hg$Symbol) %>% slice_max(QN_avg, n = 1) %>% pull(Symbol)
gene_mm <- filter(summ_df_mm, Symbol %!in% keep_mm$Symbol) %>% slice_max(QN_avg, n = 1) %>% pull(Symbol)
vbplot(dat_l[ids_hg], gene_hg)
vbplot(dat_l[ids_mm], gene_mm)


# Plotting average expression of counts

# TODO: w/o keep_df
plot_avg_heatmap <- function(summ_l, keep_df) {
  
  # gene_order <- arrange(keep_df, desc(QN_avg))$Symbol
  # gene_order <- filter(summ_l$Summ_df, N_msr == max(N_msr)) %>% arrange(desc(QN_avg)) %>% pull(Symbol)
  gene_order <- slice_min(summ_l$Summ_df, RP, n = 500) %>% arrange(desc(QN_avg)) %>% pull(Symbol)
  
  plot_mat <- summ_l$QN_Avg
  plot_mat <- plot_mat[gene_order, ]
  plot_mat <- cbind(plot_mat, Average = summ_l$Summ_df[gene_order, "QN_avg"])
  
  # plot_mat <- scale(plot_mat)
  # plot_mat[plot_mat > 3] <- 3
  
  pheatmap(plot_mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = RColorBrewer::brewer.pal(9, "Blues"),
           show_rownames = FALSE,
           border_color = NA,
           fontsize = 20,
           gaps_col = ncol(plot_mat) - 1)
  
}


plot_avg_heatmap(count_summ$Human, keep_hg)
plot_avg_heatmap(count_summ$Mouse, keep_mm)





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
plot_mat <- cbind(plot_mat, Average = count_summ$Mouse$QN_Avg[, "GSE118020"])


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

