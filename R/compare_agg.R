# Collecting aggregates across microglia...

library(tidyverse)
library(pheatmap)
source("R/00_config.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")

k <- 200
set.seed(5)

dat_dir <- file.path(data_out_dir, "CSCORE")
mcg_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)

pc_hg <- read.delim(ref_hg_path)
pc_mm <- read.delim(ref_mm_path)
tfs_hg <- read.delim(tfs_hg_path)
tfs_mm <- read.delim(tfs_mm_path)

allrank_sub_mm <- readRDS("~/mcgdat/Cormats/Mm_pcor/aggregate_cormat_allrank_filter_lenient_mm.RDS")
allrank_mm <- readRDS("~/mcgdat/Cormats/Mm_pcor/aggregate_cormat_allrank_mm.RDS")
scor_allrank_mm <- readRDS(file.path(cmat_dir_mm, "aggregate_cormat_allrank_filter_scor_mm.RDS"))
# scor_fz_mm <- readRDS(file.path(cmat_dir_mm, "aggregate_cormat_FZ_filter_scor_mm.RDS"))
fz_mm <- readRDS("~/mcgdat/Cormats/Mm_pcor/aggregate_cormat_FZ_mm.RDS")
csc_mm <- readRDS("~/mcgdat/CSCORE_aggregate_microglia_mm.RDS")
grn_avg_mm <- readRDS("~/mcgdat/GRNBoost2_average_mat_mm.RDS")
grn_allrank_mm <- readRDS("~/mcgdat/GRNBoost2_allrank_mat_mm.RDS")

# List of gene count measurement summaries
count_summ <- readRDS(mcg_count_summ_list_path)


# Filter genes: require measured in at least a third of datasets
keep_mm <- count_summ$Mouse$Filter_genes 
keep_tfs_mm <- intersect(tfs_mm$Symbol, keep_mm)
msr_mm <- count_summ$Mouse$Summ_df


identical(rownames(csc_mm), keep_mm)
identical(rownames(allrank_sub_mm$Agg_mat), keep_mm)
identical(rownames(grn_avg_mm), keep_mm)
identical(rownames(grn_allrank_mm), keep_mm)
identical(colnames(grn_avg_mm), keep_tfs_mm)


# R = All rank
# FZ = Fisher's Z transformation
# CSC = CSCORE  ** Only has threshold genes
# GRN = GRNBoost2  ** Only has TFs
# Scor = Spearman's cor  ** Only has threshold genes
# *Sub = Using threshold genes


comparisons <- list(
  R_vs_RSub = list(allrank_mm$Agg_mat[keep_mm, keep_mm], allrank_sub_mm$Agg_mat[keep_mm, keep_mm]), 
  FZ_vs_R = list(fz_mm$Agg_mat, allrank_mm$Agg_mat),
  FZSub_vs_RSub = list(fz_mm$Agg_mat[keep_mm, keep_mm], allrank_sub_mm$Agg_mat[keep_mm, keep_mm]),
  FZSub_vs_CSC = list(fz_mm$Agg_mat[keep_mm, keep_mm], csc_mm),
  RSub_vs_CSC = list(allrank_sub_mm$Agg_mat, csc_mm),
  RSub_vs_GRN = list(allrank_sub_mm$Agg_mat[, keep_tfs_mm], grn_avg_mm),
  FZSub_vs_GRN = list(fz_mm$Agg_mat[keep_mm, keep_tfs_mm], grn_avg_mm),
  CSC_vs_GRN = list(csc_mm[, keep_tfs_mm], grn_avg_mm),
  ScorR_vs_Rsub = list(scor_allrank_mm$Agg_mat, allrank_sub_mm$Agg_mat)
  # ScorFZ_vs_FZ = list(scor_fz_mm$Agg_mat, fz_mm$Agg_mat[keep_mm, keep_mm])
)


# Perform all comparisons and store results
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
  ScorR_vs_Rsub = res$ScorR_vs_Rsub
  # ScorFZ_vs_FZ = tt10
)


compare_df <- compare_df %>% 
  mutate(
    Avg_overlap = rowMeans(select_if(., is.numeric), na.rm = TRUE),
    Is_TF = Symbol %in% keep_tfs_mm
    ) %>% 
  left_join(msr_mm, by = "Symbol")



# Which aggregate profiles are most consistent?
top_overlap <- compare_df %>% 
  arrange(desc(Avg_overlap)) %>% 
  pull(Symbol)

top_tfs_overlap <- top_overlap[top_overlap %in% keep_tfs_mm]



ggplot(compare_df, aes(x = R_vs_RSub)) +
  geom_histogram(bins = 30, fill = "#38761dff", alpha = 0.8) +
  xlab("Top200") +
  ylab("Count of genes") +
  ggtitle("Ranking with all genes vs. thresholded genes") +
  theme_classic() +
  theme(text = element_text(size = 20))


null_topk <- pair_shuffle_topk(
  comparisons$ScorR_vs_Rsub[[1]][, keep_tfs_mm],
  comparisons$ScorR_vs_Rsub[[2]][, keep_tfs_mm],
  k = k,
  ncores = ncore
)


ggplot(filter(compare_df, Is_TF), aes(x = ScorR_vs_Rsub)) +
  geom_histogram(bins = 30, fill = "#38761dff", alpha = 0.8) +
  geom_vline(xintercept = mean(null_topk), linetype = "dashed", colour = "darkgrey") +
  xlab("Top200") +
  ylab("Count of genes") +
  ggtitle("Ranked Pcor vs. Ranked Scor") +
  theme_classic() +
  theme(text = element_text(size = 20))



plot_mat <- as.matrix(select_if(compare_df, is.numeric))
# plot_mat <- plot_mat[names(avg_topk), ]
# plot_mat <- plot_mat[names(avg_tfs_topk), ]
plot_mat <- plot_mat[names(avg_tfs_topk)[1:100], ]


pheatmap(plot_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE)


plot(msr_mm[keep_mm, "QN_avg"], avg_topk[keep_mm])
plot(msr_mm[keep_tfs_mm, "QN_avg"], avg_topk[keep_tfs_mm])

plot(msr_mm[keep_mm, "N_msr"], avg_topk[keep_mm])
plot(msr_mm[keep_tfs_mm, "N_msr"], avg_topk[keep_tfs_mm])



check_gene <- "Mef2c"

df <- data.frame(
  Symbol = keep_mm,
  Allrank = allrank_mm$Agg_mat[keep_mm, check_gene],
  Allrank_sub = allrank_sub_mm$Agg_mat[keep_mm, check_gene],
  Scor_allrank = scor_allrank_mm$Agg_mat[keep_mm, check_gene],
  # Scor_FZ = scor_fz_mm$Agg_mat[keep_mm, check_gene],
  FZ = fz_mm$Agg_mat[keep_mm, check_gene],
  CSCORE = csc_mm[keep_mm, check_gene],
  GRN_avg = grn_avg_mm[, check_gene],
  GRN_allrank = grn_allrank_mm[, check_gene]
  ) %>% 
  left_join(msr_mm, by = "Symbol")
rownames(df) <- keep_mm
df <- filter(df, Symbol != check_gene)


cor(select_if(df, is.numeric), method = "spearman")
colwise_topk_intersect(as.matrix(select_if(df, is.numeric)), k = 200)
colwise_topk_intersect(-as.matrix(select_if(df, is.numeric)), k = 200)


plot(df$FZ, df$CSCORE)
plot(df$Allrank, df$CSCORE)


plot(df$FZ, df$Allrank)
plot(df$FZ, df$Allrank_sub)
plot(df$FZ, df$GRN_avg)
plot(df$FZ, df$GRN_allrank)

plot(df$CSCORE, df$Allrank_sub)
plot(df$CSCORE, df$GRN_avg)
plot(df$CSCORE, df$GRN_allrank)

plot(df$Allrank, df$Allrank_sub)

plot(df$FZ, df$N_msr)
plot(df$FZ, df$QN_avg)

plot(df$Allrank, df$N_msr)
plot(df$Allrank, df$QN_avg)

plot(df$CSCORE, df$N_msr)
plot(df$CSCORE, df$QN_avg)


plot(df$Allrank_sub, df$Scor_allrank)
plot(df$CSCORE, df$QN_avg)


# https://pmc.ncbi.nlm.nih.gov/articles/PMC6377292/
# https://www.frontiersin.org/journals/molecular-neuroscience/articles/10.3389/fnmol.2022.1072046/full
tfs <- c("Runx1", "Spi1", "Mef2a", "Mef2c", "Sall1", "Mafb", "Irf8")
pmat <- fz_mm$Agg_mat[keep_mm, tfs]
pmat[is.na(pmat)] <- 0
pmat[tfs, tfs] <- 0
pheatmap::pheatmap(pmat)


keep_tfs_mm <- intersect(keep_mm, tfs_mm$Symbol)
tt1 <- fz_mm$Agg_mat[keep_mm, keep_tfs_mm]
tt1[is.na(tt1) | is.infinite(tt1)] <- 0
tt2 <- colwise_topk_intersect(tt1, k = 200)
tt3 <- mat_to_df(tt2, value = "Topk")



# Load cormat

keep_tfs <- colnames(avg_grn_mat)
keep_genes <- rownames(avg_grn_mat)


cmat_l <- lapply(ids_mm, function(x) {
  fread_to_mat(file.path(cmat_dir_mm, paste0(x, "_cormat.tsv")), 
               genes = keep_genes, 
               sub_genes = keep_tfs)
})

msr_l <- lapply(cmat_l, function(x) !is.na(x))
msr_mat <- Reduce("+", msr_l)
msr_tfs <- sapply(keep_tfs, function(x) msr_mat[x, x])




fz_l <- lapply(cmat_l, function(x) {
  x[is.na(x)] <- 0
  aggtools::fisherz(x)
})
names(fz_l) <- ids_mm

fz_mat <- Reduce("+", fz_l) / msr_mat

# TODO: double check this
fz_mat[is.infinite(fz_mat)] <- 0
fz_mat[is.na(fz_mat)] <- 0

# which(is.infinite(fz_mat), arr.ind = TRUE)
# which(is.na(fz_mat), arr.ind = TRUE)


k <- 200
compare_topk <- pair_colwise_topk(fz_mat, avg_grn_mat, k = k, ncores = ncore)
compare_topk_abs <- pair_colwise_topk(abs(fz_mat), avg_grn_mat, k = k, ncores = ncore)
compare_scor <- pair_colwise_cor(fz_mat, avg_grn_mat, ncores = ncore)
compare_scor_abs <- pair_colwise_cor(abs(fz_mat), avg_grn_mat, ncores = ncore)


compare_df <- data.frame(
  Symbol = keep_tfs,
  N_msr = msr_tfs,
  Topk = compare_topk,
  Topk_abs = compare_topk_abs,
  Scor = compare_scor,
  Scor_abs = compare_scor_abs
)
compare_df <- filter(compare_df, N_msr >= 5)





# Hists of the similarity between the averaged network and pcor
ggplot(compare_df, aes(x = Topk)) +
  geom_histogram(bins = 100) +
  xlab(paste0("Top", k, "between averaged GRN and FZ")) +
  ylab("Count TFs") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))




check_tf <- "Spi1"


check_df <- data.frame(
  Symbol = rownames(fz_mat),
  N_msr = msr_mat[rownames(fz_mat), check_tf],
  FZ = fz_mat[, check_tf],
  GRN = avg_grn_mat[, check_tf]
)

check_df <- filter(check_df, N_msr >= 5)


ggplot(check_df, aes(x = FZ, y = GRN)) +
  geom_point(shape = 21, size = 2.2) +
  ggtitle(check_tf) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "firebrick") +
  theme_classic() +
  theme(text = element_text(size = 20))

cor(check_df$FZ, check_df$GRN, use = "pairwise.complete.obs")


slice_max(check_df, GRN, n = k) %>% filter(FZ < 0.01 & FZ > -0.01)
check_gene <- "Coro1a"


check_data <- lapply(ids_mm, function(x) {
  
  fz <- fz_l[[x]]
  grn <- avg_grn_l[[x]]
  
  if ((check_gene %!in% rownames(fz_mat)) || (check_tf %!in% colnames(fz_mat))) {
    return(NA)
  }
  
  if ((check_gene %!in% rownames(grn)) || (check_tf %!in% colnames(grn))) {
    return(NA)
  }
  
  
  fz_tf <- fz[, check_tf, drop = FALSE]
  fz_tf <- fz_tf[setdiff(rownames(fz_tf), check_tf), , drop = FALSE]
  rank_fz <- aggtools::colrank_mat(fz_tf)
  
  grn_tf <- grn[, check_tf, drop = FALSE]
  rank_grn <- aggtools::colrank_mat(grn_tf)
  
  data.frame(ID = x,
             FZ = fz_tf[check_gene, ], 
             GRN = grn_tf[check_gene, ],
             Rank_FZ = rank_fz[check_gene, ],
             GRN_rank = rank_grn[check_gene, ])
})

check_data <- check_data[!is.na(check_data)]
check_data <- do.call(rbind, check_data) %>% arrange(GRN_rank)


ggplot(check_data, aes(x = FZ)) +
  geom_histogram(bins = 10) +
  geom_vline(xintercept = 0, col = "firebrick") +
  xlim(c(-0.2, 0.2)) +
  ggtitle(paste(check_tf, check_gene)) +
  theme_classic() +
  theme(text = element_text(size = 20))


ggplot(check_data, aes(x = GRN)) +
  geom_histogram(bins = 10) +
  geom_vline(xintercept = 0, col = "firebrick") +
  # xlim(c(-0.2, 0.2)) +
  ggtitle(paste(check_tf, check_gene)) +
  theme_classic() +
  theme(text = element_text(size = 20))



# Compare to median GRN

# bind_cols(lapply(avg_grn_l, function(grn) {
#   
#   if ((check_gene %!in% rownames(grn)) || (check_tf %!in% colnames(grn))) {
#     return(NA)
#   }
#   
#   grn[, check_tf, drop = FALSE]
#   
# }))

