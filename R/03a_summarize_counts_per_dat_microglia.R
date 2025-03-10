## Summarize the CPM-normalized count data across microglia datasets, saving
## this info out as a list. Also perform basic exploration of this information.
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
library(ComplexHeatmap)
library(preprocessCore)
library(matrixStats)
library(ggrepel)
source("R/00_config.R")
source("R/utils/functions.R")

# Dataset meta and human/mouse data IDs
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)

tfs_hg <- read.delim(tfs_hg_path)
tfs_mm <- read.delim(tfs_mm_path)
pc_ortho <- read.delim(pc_ortho_path)

# List of count matrices and their metadata
dat_l <- readRDS(mcg_dat_path)


if (!file.exists(mcg_count_summ_list_path)) {
  
  # Generate summaries for human and mouse separately
  count_summ_hg <- summarize_gene_counts(dat_l[ids_hg])
  count_summ_mm <- summarize_gene_counts(dat_l[ids_mm])
  
  # Save out
  saveRDS(list(Human = count_summ_hg, Mouse = count_summ_mm), mcg_count_summ_list_path)
  write.table(count_summ_hg$Summ_df, sep = "\t", quote = FALSE, row.names = FALSE, file = mcg_count_summ_table_hg)
  write.table(count_summ_mm$Summ_df, sep = "\t", quote = FALSE, row.names = FALSE, file = mcg_count_summ_table_mm)

} else {
  
  count_summ <- readRDS(mcg_count_summ_list_path)

}



summ_hg <- count_summ$Human$Summ_df
summ_mm <- count_summ$Mouse$Summ_df


# Min number of experiments a gene needs to be measured in to be kept
min_exp_hg <- floor(length(ids_hg) * min_msr_frac)
min_exp_mm <- floor(length(ids_mm) * min_msr_frac)



# Organize TF and ortho expression
# ------------------------------------------------------------------------------


summ_hg$Is_TF <- summ_hg$Symbol %in% tfs_hg$Symbol
summ_mm$Is_TF <- summ_mm$Symbol %in% tfs_mm$Symbol


ortho_df <- pc_ortho %>% 
  left_join(summ_hg, by = c("Symbol_hg" = "Symbol")) %>% 
  left_join(summ_mm, by = c("Symbol_mm" = "Symbol"),
            suffix = c("_hg", "_mm"))


# Identifying the TFs that are highly expressed in either species  
label_tfs_hg <- ortho_df %>% 
  filter(Is_TF_hg) %>% 
  slice_max(QN_avg_hg, n = 15)

label_tfs_mm <- ortho_df %>% 
  filter(Is_TF_mm) %>% 
  slice_max(QN_avg_mm, n = 15)


label_tfs <- union(label_tfs_hg$Symbol_hg, label_tfs_mm$Symbol_hg)


ortho_df_filt <- ortho_df %>% 
  filter(N_msr_hg >= min_exp_hg | N_msr_mm >= min_exp_mm) %>% 
  mutate(Is_TF = Is_TF_hg & Is_TF_mm,
         Is_top_TF = Symbol_hg %in% label_tfs)


ortho_expr_cor <- cor(ortho_df_filt$QN_avg_mm, ortho_df_filt$QN_avg_hg)



# Identifying candidate species-specific expression
label_diff_hg <- ortho_df_filt %>% 
  slice_min(RP_hg, n = 1000) %>% 
  slice_max(RP_mm, n = 100) %>% 
  arrange(desc(QN_avg_hg))


label_diff_mm <- ortho_df_filt %>% 
  slice_min(RP_mm, n = 1000) %>% 
  slice_max(RP_hg, n = 100) %>% 
  arrange(desc(QN_avg_mm))


# Add label columns for plotting
ortho_df_filt <- ortho_df_filt %>% 
  mutate(
    Diff_hg = Symbol_hg %in% label_diff_hg$Symbol_hg,
    Diff_label_hg = Symbol_hg %in% label_diff_hg$Symbol_hg[1:5],
    Diff_mm = Symbol_hg %in% label_diff_mm$Symbol_hg,
    Diff_label_mm = Symbol_hg %in% label_diff_mm$Symbol_hg[1:5]
  )




# Plots
# ------------------------------------------------------------------------------


# Scatter of ortho expression, with trend line and TFs emphasized
p1a <- ggplot(ortho_df_filt, aes(x = QN_avg_mm, y = QN_avg_hg)) +
  geom_point(data = filter(ortho_df_filt, !Is_TF),
             shape = 21, size = 2.4, alpha = 0.2) +
  geom_point(data = filter(ortho_df_filt, Is_TF),
  shape = 21, size = 3.4) +
  geom_smooth(method = "lm", colour = "royalblue") +
  xlab("Mouse mean log2 CPM") +
  ylab("Human mean log2 CPM") +
  theme_classic() +
  theme(text = element_text(size = 25),
        plot.margin = margin(c(10, 20, 10, 10)))


# Scatter of ortho expression, with top TFs labeled
p1b <- ggplot(ortho_df_filt, aes(x = QN_avg_mm, y = QN_avg_hg)) +
  geom_point(data = filter(ortho_df_filt, !Is_TF), 
             shape = 21, size = 2.4, alpha = 0.2) +
  geom_point(data = filter(ortho_df_filt, Is_TF), 
             shape = 21, size = 3.4, colour = "firebrick") +
  geom_text_repel(data = filter(ortho_df_filt, Is_top_TF),
                  aes(label = Symbol_hg),
                  size = 6, max.overlaps = length(label_tfs)) +
  xlab("Mouse mean log2 CPM") +
  ylab("Human mean log2 CPM") +
  theme_classic() +
  theme(text = element_text(size = 25),
        plot.margin = margin(c(10, 20, 10, 10)))


# Scatter of ortho expression labelling species-specific expression
p1c <- ggplot(ortho_df_filt, aes(x = QN_avg_mm, y = QN_avg_hg)) +
  geom_point(data = filter(ortho_df_filt, !Diff_hg & !Diff_mm), 
             shape = 21, size = 2.4, alpha = 0.2) +
  geom_point(data = filter(ortho_df_filt, Diff_hg), 
             shape = 21, size = 3.4, colour = "royalblue") +
  geom_point(data = filter(ortho_df_filt, Diff_mm), 
             shape = 21, size = 3.4, colour = "goldenrod") +
  geom_text_repel(data = filter(ortho_df_filt, Diff_label_hg | Diff_label_mm),
                  aes(label = Symbol_hg),
                  size = 6, max.overlaps = 10) +
  xlab("Mouse mean log2 CPM") +
  ylab("Human mean log2 CPM") +
  theme_classic() +
  theme(text = element_text(size = 25),
        plot.margin = margin(c(10, 20, 10, 10)))



# Plotting the averaged+pseudobulked+quantile norm'd expression vs count msr'd,
# overlaid with line indicating min msr'd count filter
p2a <- ggplot(summ_hg, aes(x = N_msr, y = QN_avg)) +
  geom_point(shape = 21, size = 2.4) +
  geom_vline(xintercept = min_exp_hg - 0.5, 
             colour = "firebrick", 
             linetype = "dashed") +
  ggtitle("Human") +
  ylab("Average Log2 CPM") +
  xlab("Count measured") +
  theme_classic() +
  theme(text = element_text(size = 25))


# Mouse
p2b <- ggplot(summ_mm, aes(x = N_msr, y = QN_avg)) +
  geom_point(shape = 21, size = 2.4) +
  geom_vline(xintercept = min_exp_mm - 0.5, 
             colour = "firebrick", 
             linetype = "dashed") +
  ggtitle("Mouse") +
  ylab("Average Log2 CPM") +
  xlab("Count measured") +
  theme_classic() +
  theme(text = element_text(size = 25))



# Violin+box plot of CPM counts for given gene


vbplot <- function(dat_l, gene) {
  
  ids <- names(dat_l)
  
  gene_l <- lapply(ids, function(x) {
    data.frame(Log2_count = log2(dat_l[[x]]$Mat[gene, ] + 1), ID = x)
  })
  
  gene_df <- do.call(rbind, gene_l)
  
  ggplot(gene_df, aes(x = ID, y = Log2_count)) +
    geom_violin(fill = "slategrey") +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    ylab(bquote(~Log[2]~ "CPM")) +
    ggtitle(gene) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 25),
          plot.margin = margin(c(10, 20, 10, 10)))
}



# The max average gene
max_hg <- slice_max(summ_hg, QN_avg, n = 1) %>% pull(Symbol)
max_mm <- slice_max(summ_mm, QN_avg, n = 1) %>% pull(Symbol)
p3a <- vbplot(dat_l[ids_hg], max_hg)
p3b <- vbplot(dat_l[ids_mm], max_mm)

# The lowest average gene that was kept
low_hg <- summ_hg %>% filter(Symbol %in% count_summ$Human$Filter_genes) %>% slice_min(QN_avg) %>% pull(Symbol)
low_mm <- summ_mm %>% filter(Symbol %in% count_summ$Mouse$Filter_genes) %>% slice_min(QN_avg) %>% pull(Symbol)
p3c <- vbplot(dat_l[ids_hg], low_hg)
p3d <- vbplot(dat_l[ids_mm], low_mm)

# The highest average gene that was removed
rm_hg <- summ_hg %>% filter(Symbol %!in% count_summ$Human$Filter_genes) %>% slice_max(QN_avg, n = 1) %>% pull(Symbol)
rm_mm <- summ_mm %>% filter(Symbol %!in% count_summ$Mouse$Filter_genes) %>% slice_max(QN_avg, n = 1) %>% pull(Symbol)
p3e <- vbplot(dat_l[ids_hg], rm_hg)
p3f <- vbplot(dat_l[ids_mm], rm_mm)



# Heatmap of average expression of counts

plot_avg_heatmap <- function(summ_l) {
  
  keep_genes <- summ_l$Filter_genes
  
  gene_order <- filter(summ_l$Summ_df, Symbol %in% keep_genes) %>% arrange(desc(QN_avg)) %>% pull(Symbol)
  # gene_order <- filter(summ_l$Summ_df, N_msr == max(N_msr)) %>% arrange(desc(QN_avg)) %>% pull(Symbol)
  # gene_order <- slice_min(summ_l$Summ_df, RP, n = 100) %>% arrange(desc(QN_avg)) %>% pull(Symbol)
  
  plot_mat <- summ_l$QN_Avg
  plot_mat <- cbind(plot_mat, Average = rowMeans(plot_mat, na.rm = TRUE))
  plot_mat <- plot_mat[gene_order, ]
  
  
  Heatmap(plot_mat,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d'),
          show_row_names = FALSE,
          column_split = c(rep("Datasets", ncol(plot_mat) - 1), rep("Mean", 1)),
          heatmap_legend_param = list(title = "Log2 CPM")
          )
  
}


# Human
png(height = 9, width = 7.5, units = "in", res = 300,
    file = file.path(plot_dir, "avg_measurement_heatmap_human.png"))

plot_avg_heatmap(count_summ$Human)

dev.off()

# Mouse
png(height = 6, width = 9, units = "in", res = 300,
    file = file.path(plot_dir, "avg_measurement_heatmap_mouse.png"))

plot_avg_heatmap(count_summ$Mouse)

dev.off()


# Heatmap of binary measurement, with bar plot counting measurement across data


plot_msr_heatmap <- function(summ_l) {
  
  plot_mat <- summ_l$Msr * 1  # logical -> integer
  order_genes <- summ_l$Summ_df %>% arrange(desc(N_msr)) %>% pull(Symbol)
  plot_mat <- plot_mat[order_genes, ]
  
  bar <- rowAnnotation(N_msr = anno_barplot(summ_l$Summ_df[order_genes, "N_msr"]))
  cols <- c("white", "black")
  
  Heatmap(plot_mat,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          name = "Measurement", 
          col = cols,
          right_annotation = bar,
          use_raster = TRUE)
}



# Human
png(height = 9, width = 7.5, units = "in", res = 300,
    file = file.path(plot_dir, "binary_measurement_heatmap_human.png"))

plot_msr_heatmap(count_summ$Human)

dev.off()


# Mouse
png(height = 9, width = 9, units = "in", res = 300,
    file = file.path(plot_dir, "binary_measurement_heatmap_mouse.png"))

plot_msr_heatmap(count_summ$Mouse)

dev.off()
