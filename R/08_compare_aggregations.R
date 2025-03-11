## This script collects various coexpression aggregation in mouse microglia (as
## there is more data than human) to compare/perform interactive exploration
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
source("R/00_config.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")

k <- 200  # Cutoff for TopK overlap
set.seed(5)

# Data meta and IDs
mcg_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)
ids_mm <- filter(mcg_meta, Species == "Mouse")$ID

# Gene tables
pc_mm <- read.delim(ref_mm_path)
tfs_mm <- read.delim(tfs_mm_path)

# Load aggregates
# --
# R = All rank
# FZ = Fisher's Z transformation
# CSC = CSCORE  ** Only has threshold genes
# GRN = GRNBoost2  ** Only has TFs
# Scor = Spearman's cor  ** Only has threshold genes
# *Sub = Using threshold genes
# --
allrank_sub_mm <- readRDS(mcg_allrank_filt_mm_path)
allrank_mm <- readRDS(mcg_allrank_mm_path)
scor_allrank_mm <- readRDS(mcg_allrank_scor_mm_path)
scor_fz_mm <- readRDS(mcg_fz_scor_mm_path)
fz_mm <- readRDS(mcg_fz_mm_path)
csc_mm <- readRDS(cscore_mm_path)
grn_avg_mm <- readRDS(grn_avg_mm_path)
grn_allrank_mm <- readRDS(grn_allrank_mm_path)

# List of microglia count matrices and their meta
mcg_dat <- readRDS(mcg_dat_path)

# List of gene count measurement summaries
count_summ <- readRDS(mcg_count_summ_list_path)


# Filter genes: require measured in at least a third of datasets
keep_mm <- count_summ$Mouse$Filter_genes 
keep_tfs_mm <- intersect(tfs_mm$Symbol, keep_mm)
msr_mm <- count_summ$Mouse$Summ_df


stopifnot(identical(rownames(csc_mm), keep_mm))
stopifnot(identical(rownames(allrank_sub_mm$Agg_mat), keep_mm))
stopifnot(identical(rownames(grn_avg_mm), keep_mm))
stopifnot(identical(rownames(grn_allrank_mm), keep_mm))
stopifnot(identical(colnames(grn_avg_mm), keep_tfs_mm))



# Load specific comparisons of aggregates into a list
comparisons <- list(
  R_vs_RSub = list(allrank_mm$Agg_mat[keep_mm, keep_mm], allrank_sub_mm$Agg_mat[keep_mm, keep_mm]), 
  FZ_vs_R = list(fz_mm$Agg_mat, allrank_mm$Agg_mat),
  FZSub_vs_RSub = list(fz_mm$Agg_mat[keep_mm, keep_mm], allrank_sub_mm$Agg_mat[keep_mm, keep_mm]),
  FZSub_vs_CSC = list(fz_mm$Agg_mat[keep_mm, keep_mm], csc_mm),
  RSub_vs_CSC = list(allrank_sub_mm$Agg_mat, csc_mm),
  RSub_vs_GRN = list(allrank_sub_mm$Agg_mat[, keep_tfs_mm], grn_avg_mm),
  FZSub_vs_GRN = list(fz_mm$Agg_mat[keep_mm, keep_tfs_mm], grn_avg_mm),
  CSC_vs_GRN = list(csc_mm[, keep_tfs_mm], grn_avg_mm),
  ScorR_vs_Rsub = list(scor_allrank_mm$Agg_mat, allrank_sub_mm$Agg_mat),
  ScorFZ_vs_FZ = list(scor_fz_mm$Agg_mat, fz_mm$Agg_mat[keep_mm, keep_mm])
)


# Perform TopK overlap of all comparisons
res <- lapply(comparisons, function(mats) pair_colwise_topk(mats[[1]], mats[[2]], k, ncore))


# Helper to fill out gene vector with NAs for GRN
fill_gene_vec <- function(vec, keep_genes) {
  fill_vec <- setNames(rep(NA, length(keep_genes)), keep_genes)
  fill_vec[names(vec)] <- vec
  return(fill_vec)
}


# Organize topk overlap into df
compare_df <- data.frame(
  Symbol = keep_mm,
  R_vs_RSub = res$R_vs_RSub,
  FZ_vs_R = res$FZ_vs_R[keep_mm],  # This uses all genes for rank, subset to threshold
  FZSub_vs_RSub = res$FZSub_vs_RSub,
  FZSub_vs_CSC = res$FZSub_vs_CSC,
  RSub_vs_CSC = res$RSub_vs_CSC,
  RSub_vs_GRN = fill_gene_vec(res$RSub_vs_GRN, keep_mm),
  FZSub_vs_GRN = fill_gene_vec(res$FZSub_vs_GRN, keep_mm),
  CSC_vs_GRN = fill_gene_vec(res$CSC_vs_GRN, keep_mm),
  ScorR_vs_Rsub = res$ScorR_vs_Rsub,
  ScorFZ_vs_FZ = res$ScorFZ_vs_FZ
)


# Adding average overlap and TF status
compare_df <- compare_df %>% 
  mutate(
    Avg_overlap = rowMeans(select_if(., is.numeric), na.rm = TRUE),
    Is_TF = Symbol %in% keep_tfs_mm
    ) %>% 
  left_join(msr_mm, by = "Symbol")


# Which aggregate gene profiles are most consistent? Pull for downstream ordering
top_overlap <- compare_df %>% 
  arrange(desc(Avg_overlap)) %>% 
  pull(Symbol)

top_tfs_overlap <- top_overlap[top_overlap %in% keep_tfs_mm]


# Generating null overlap between all rank Scor and Pcor
null_topk1 <- pair_shuffle_topk(
  comparisons$ScorR_vs_Rsub[[1]][, keep_tfs_mm],
  comparisons$ScorR_vs_Rsub[[2]][, keep_tfs_mm],
  k = k,
  ncores = ncore
)


# This is used to isolate a single TF to explore its gene rankings and merge
# with the gene measurement info

check_gene <- "Tcf4"

check_df <- data.frame(
  Symbol = keep_mm,
  Allrank = allrank_mm$Agg_mat[keep_mm, check_gene],
  Allrank_sub = allrank_sub_mm$Agg_mat[keep_mm, check_gene],
  Scor_allrank = scor_allrank_mm$Agg_mat[keep_mm, check_gene],
  Scor_FZ = scor_fz_mm$Agg_mat[keep_mm, check_gene],
  FZ = fz_mm$Agg_mat[keep_mm, check_gene],
  CSCORE = csc_mm[keep_mm, check_gene],
  GRN_avg = grn_avg_mm[, check_gene],
  GRN_allrank = grn_allrank_mm[, check_gene]
  ) %>% 
  left_join(msr_mm, by = "Symbol")
rownames(check_df) <- keep_mm
check_df <- filter(check_df, Symbol != check_gene)  # Remove self (perfect cor)


check_cor <- cor(select_if(check_df, is.numeric), method = "spearman")
check_topk <- colwise_topk_intersect(as.matrix(select_if(check_df, is.numeric)), k = 200)
check_bottomk <- colwise_topk_intersect(-as.matrix(select_if(check_df, is.numeric)), k = 200)




# Plots
# ------------------------------------------------------------------------------


# Hist showing the overlap between allrank +/- filtering genes
p1 <- ggplot(compare_df, aes(x = R_vs_RSub)) +
  geom_histogram(bins = 30, fill = "#38761dff", alpha = 0.8) +
  xlab("Top200") +
  ylab("Count of genes") +
  ggtitle("Ranking Pcor with all genes vs. thresholded genes") +
  theme_classic() +
  theme(text = element_text(size = 20))


# Hist showing the overlap between allrank Pcor vs Scor (+filtering) overlaid with null
p2 <- ggplot(filter(compare_df, Is_TF), aes(x = ScorR_vs_Rsub)) +
  geom_histogram(bins = 30, fill = "#54278f", alpha = 0.8) +
  geom_vline(xintercept = mean(null_topk1), linetype = "dashed", colour = "darkgrey") +
  xlab("Top200") +
  ylab("Count of genes") +
  ggtitle("Ranked Pcor vs. Ranked Scor (TFs only)") +
  theme_classic() +
  theme(text = element_text(size = 20))


# Heatmap of overlaps
plot_mat <- plot_mat[top_tfs_overlap, ]
pal <- colorRampPalette(c("white", "#ca0020"))(10)

pheatmap(plot_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         color = pal,
         border_color = NA,
         gaps_col = rep((ncol(plot_mat) - 1), 3),
         fontsize = 20)


# Scatter plots exploring TF expression versus agreement across aggregates
p3a <- ggplot(filter(compare_df, Is_TF), aes(x = QN_avg, y = Avg_overlap)) +
  geom_point(shape = 21) +
  geom_smooth(method = "lm") +
  xlab("Mean log2 CPM") +
  ylab("Mean Top200 across strategies") +
  theme_classic() +
  theme(text = element_text(size = 20))


p3b <- ggplot(filter(compare_df, Is_TF), aes(x = N_msr, y = Avg_overlap)) +
  geom_point(shape = 21) +
  geom_smooth(method = "lm") +
  xlab("Count measured") +
  ylab("Mean Top200 across strategies") +
  theme_classic() +
  theme(text = element_text(size = 20))



# Scatter plots comparing selected aggregates for check TF

plot_scatter <- function(check_df, check_tf, x, y, xlab, ylab) {
  
  ggplot(check_df, aes(x = !!sym(x), y = !!sym(y))) +
  geom_point(shape = 21) +
  xlab(xlab) +
  ylab(ylab) +
  ggtitle(check_gene) +
  theme_classic() +
  theme(text = element_text(size = 30),
        plot.margin = margin(10, 20, 10, 10))

}


p4a <- plot_scatter(check_df, check_tf, x = "Allrank_sub", y = "Scor_allrank", 
                    xlab = "Ranked Pcor", ylab = "Ranked Scor")

p4b <- plot_scatter(check_df, check_tf, x = "FZ", y = "Scor_FZ", 
                    xlab = "Average Fisher's Z Pcor", ylab = "Average Fisher's Z Scor")

p4c <- plot_scatter(check_df, check_tf, x = "Allrank_sub", y = "CSCORE", 
                    xlab = "Ranked Pcor", ylab = "Ranked CSCORE")

p4d <- plot_scatter(check_df, check_tf, x = "Allrank_sub", y = "FZ", 
                    xlab = "Ranked Pcor", ylab = "Average Fisher's Z Pcor")

p4e <- plot_scatter(check_df, check_tf, x = "Allrank_sub", y = "GRN_avg", 
                    xlab = "Ranked Pcor", ylab = "Average GRNBoost2")

p4f <- plot_scatter(check_df, check_tf, x = "FZ", y = "GRN_avg", 
                    xlab = "Average Fisher's Z Pcor", ylab = "Average GRNBoost2")
