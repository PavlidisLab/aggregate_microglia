# Comparing aggregate ranks using HPA (single dataset with many cell types) as
# well as microglial aggregate (single cell type multiple datasets)
# ------------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Shuffle procedure for null comparison
set.seed(78)

#
mcg_p_path <- file.path(amat_dir, "Microglia", "Mm_pcor", "aggregate_matrix.tsv")
mcg_s_path <- file.path(amat_dir, "Microglia", "Mm_scor", "aggregate_matrix.tsv")
mcg_na_path <- file.path(amat_dir, "Microglia", "Mm_pcor", "NA_matrix.tsv")
mcg_pneg_path <- file.path(amat_dir, "Microglia", "Mm_pcor", "count_negcor_matrix.tsv")
mcg_sneg_path <- file.path(amat_dir, "Microglia", "Mm_scor", "count_negcor_matrix.tsv")

hpa_p_path <- file.path(amat_dir, "HPA", "HPA_RSR_allrank_CPM.tsv")
hpa_s_path <- file.path(amat_dir, "HPA", "HPA_RSR_allrank_scor.tsv")
hpa_na_path <- file.path(amat_dir, "HPA", "HPA_NA_mat_CPM.tsv")


mcg_p_ids <- list.files(file.path(amat_dir, "Microglia", "Mm_pcor"),
                        full.names = TRUE, pattern = "cormat")

mcg_s_ids <- list.files(file.path(amat_dir, "Microglia", "Mm_scor"),
                        full.names = TRUE, pattern = "cormat")


#
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
tfs_hg <- read.delim(tfs_hg_path, stringsAsFactors = FALSE)
tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)

#
mcg_p <- fread_to_mat(mcg_p_path, genes = pc_mm$Symbol)
mcg_s <- fread_to_mat(mcg_s_path, genes = pc_mm$Symbol)
# hpa_p <- fread_to_mat(hpa_p_path, genes = pc_hg$Symbol)
# hpa_s <- fread_to_mat(hpa_s_path, genes = pc_hg$Symbol)

#
mcg_na <- fread_to_mat(mcg_na_path, genes = pc_mm$Symbol)
# hpa_na <- fread_to_mat(hpa_na_path, genes = pc_hg$Symbol)

#
mcg_pneg <- fread_to_mat(mcg_pneg_path, genes = pc_mm$Symbol)
mcg_sneg <- fread_to_mat(mcg_sneg_path, genes = pc_mm$Symbol)



# sim with NA check
pair_colwise_cor_sub <- function(mat1,
                                 mat2,
                                 comsr_mat,
                                 cor_method = "spearman",
                                 min_comsr,
                                 ncores = 1) {
  
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  stopifnot(identical(colnames(mat1), colnames(comsr_mat)))

  cor_l <- mclapply(1:ncol(mat1), function(x) {
    
    keep <- names(which(comsr_mat[, x] >= min_comsr))
    if (length(keep) == 0) return(NA)
    cor(mat1[keep, x], mat2[keep, x], method = cor_method, use = "everything")
    
  }, mc.cores = ncores)
  
  names(cor_l) <- colnames(mat1)
  return(unlist(cor_l))
}



#
pair_colwise_topk_sub <- function(mat1,
                                  mat2,
                                  comsr_mat,
                                  k = 200,
                                  min_comsr,
                                  ncores = 1) {
  
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  stopifnot(identical(colnames(mat1), colnames(comsr_mat)))
  
  topk_l <- mclapply(1:ncol(mat1), function(x) {
    
    keep <- names(which(comsr_mat[, x] >= min_comsr))
    if (length(keep) == 0) return(NA)
    topk_intersect(topk_sort(vec = mat1[keep, x], k = k),
                   topk_sort(vec = mat2[keep, x], k = k))
    
  }, mc.cores = ncores)
  
  names(topk_l) <- colnames(mat1)
  return(unlist(topk_l))
}


# Construct a gene x dataset ID matrix of the raw cors for the given gene
load_or_generate_rawcor <- function(in_paths,
                                    out_path,
                                    genes,
                                    sub_gene) {
  
  if (!file.exists(out_path)) {
    cor_l <- lapply(in_paths, fread_to_mat, genes = genes, sub_genes = sub_gene)
    cmat <- do.call(cbind, cor_l)
    colnames(cmat) <- in_paths
    saveRDS(cmat, out_path)
  }
  
  cmat <- readRDS(out_path)
  
  return(cmat)
}



# Microglia comparison
# ------------------------------------------------------------------------------


# Count of times a gene was measured. Uses the count of NAs, and assumes that
# some genes were never measured (all NAs) and that this max equal n datasets
# TODO: should be tracked in microglia metadata

n_dat <- max(mcg_na)
mcg_n_msr <- n_dat - diag(mcg_na)
mcg_min_msr <- ceiling(n_dat * 0.1)

# Count of times genes were measured together
comsr_mat <- n_dat - mcg_na


# Checking k (are ties resulting in lower observed K? Answer: no)
# ks <- apply(mcg_p, 2, function(x) check_k(sort(-x), k = 200))


# Calculating similarity between P/Scor aggregates (using all genes)

mcg_cor <- pair_colwise_cor(mcg_p, mcg_s, ncores = ncore)
mcg_topk <- pair_colwise_topk(mcg_p, mcg_s, ncores = ncore)
mcg_btmk <- pair_colwise_topk(-mcg_p, -mcg_s, ncores = ncore)

# Requires minimum co-measurement between genes

mcg_cor_sub <- pair_colwise_cor_sub(mat1 = mcg_p, mat2 = mcg_s, comsr_mat = comsr_mat, min_comsr = mcg_min_msr, ncores = ncore)
mcg_topk_sub <- pair_colwise_topk_sub(mat1 = mcg_p, mat2 = mcg_s, comsr_mat = comsr_mat, min_comsr = mcg_min_msr, ncores = ncore)
mcg_btmk_sub <- pair_colwise_topk_sub(mat1 = -mcg_p, mat2 = -mcg_s, comsr_mat = comsr_mat, min_comsr = mcg_min_msr, ncores = ncore)

# Comparing to shuffled null: use only min msr to prevent inflation from non
# msrd genes having identical rankings

mcg_keep <- names(mcg_n_msr[mcg_n_msr >= mcg_min_msr])
mcg_p_sub <- mcg_p[mcg_keep, mcg_keep]
mcg_s_sub <- mcg_s[mcg_keep, mcg_keep]

mcg_null_cor <- pair_shuffle_cor(mcg_p_sub, mcg_s_sub, ncores = ncore)
mcg_null_topk <- pair_shuffle_topk(mcg_p_sub, mcg_s_sub, ncores = ncore)


# Diff in count of neg cors (does one form of correlation have more negatives?)
diff_negcor <- mcg_pneg - mcg_sneg

# Summarize difference in negcors using all genes
# Higher values == pcor more often had negative negatives than scor
avg_diff_negcor <- colMeans(mcg_pneg - mcg_sneg)

# Summarize using minimum co-measurement filter
avg_diff_negcor_sub <- vapply(pc_mm$Symbol, function(x) {
  keep <- names(which(comsr_mat[, x] >= mcg_min_msr))
  keep <- setdiff(keep, x)
  mean(diff_negcor[keep, x], na.rm = TRUE)
}, numeric(1))



stopifnot(identical(names(mcg_cor), pc_mm$Symbol))
stopifnot(identical(names(mcg_cor), rownames(comsr_mat)))


mcg_df <- data.frame(
  Symbol = pc_mm$Symbol,
  TF =  pc_mm$Symbol %in% tfs_mm$Symbol,
  Cor = mcg_cor,
  TopK = mcg_topk,
  BottomK = mcg_btmk,
  Cor_sub = mcg_cor_sub, 
  TopK_sub = mcg_topk_sub, 
  BottomK_sub = mcg_btmk_sub, 
  N_msr = mcg_n_msr,
  Avg_comsr = colMeans(comsr_mat, na.rm = TRUE),
  Avg_diff_negcor = avg_diff_negcor,
  Avg_diff_negcor_sub = avg_diff_negcor_sub
)


mcg_df_sub <- filter(mcg_df, N_msr >= mcg_min_msr)



summary(mcg_df)
summary(mcg_df_sub)
summary(mcg_null_cor)
summary(mcg_null_topk)


cor(select_if(mcg_df, is.numeric), method = "spearman", use = "pairwise.complete.obs")
cor(select_if(mcg_df_sub, is.numeric), method = "spearman")


# plot of gene measurement with filter overlaid
p1a <- plot_hist(mcg_df, stat_col = "N_msr") + 
  geom_vline(xintercept = mcg_min_msr - 1, col = "red") +
  ylab("N genes") +
  xlab("N datasets")


# Combine with scatters of similarity versus measurement
p1b <- plot_grid(
  p1a, 
  qplot(mcg_df_sub, xvar = "N_msr", yvar = "TopK") + geom_smooth(method = "loess", col = "red"),
  qplot(mcg_df_sub, xvar = "N_msr", yvar = "BottomK") + geom_smooth(method = "loess", col = "red"),
  qplot(mcg_df_sub, xvar = "N_msr", yvar = "Cor") + geom_smooth(method = "loess", col = "red"),
  ncol = 2)


# Hists of similarities
p2 <- plot_grid(
  plot_hist(mcg_df_sub, stat_col = "TopK"),
  plot_hist(mcg_df_sub, stat_col = "BottomK"),
  plot_hist(mcg_df_sub, stat_col = "Cor"), 
  ncol = 1)



# boxplot(mcg_df_sub$TopK ~ mcg_df_sub$TF)
# boxplot(mcg_df_sub$Bottom_k ~ mcg_df_sub$TF)
# boxplot(mcg_df_sub$Cor ~ mcg_df_sub$TF)



# Boxplots of null vs observed
p3a <- data.frame(
  Value = c(mcg_df_sub$TopK, mcg_df_sub$TopK_sub, mcg_null_topk),
  Group = c(rep("Min_msr", nrow(mcg_df_sub)),
            rep("Min_comsr", nrow(mcg_df_sub)),
            rep("Null", length(mcg_null_topk)))
) %>% 
  ggplot(., aes(x = Group, y = Value)) +
  geom_boxplot() +
  ylab("Top 200") +
  theme_classic() +
  theme(text = element_text(size = 25))



p3b <- data.frame(
  Value = c(mcg_df_sub$Cor, mcg_df_sub$Cor_sub, mcg_null_cor),
  Group = c(rep("Min_msr", nrow(mcg_df_sub)),
            rep("Min_comsr", nrow(mcg_df_sub)),
            rep("Null", length(mcg_null_cor)))
) %>% 
  ggplot(., aes(x = Group, y = Value)) +
  geom_boxplot() +
  ylab("Spearman's cor") +
  theme_classic() +
  theme(text = element_text(size = 25))
  



# Focus on a single gene to compare aggregate and raw values
# ------------------------------------------------------------------------------



# gene <- mcg_df %>%  filter(N_msr == 1) %>% slice_max(TopK)
# gene <- mcg_df %>%  filter(N_msr == 1) %>% slice_min(TopK)
# gene <- mcg_df %>%  filter(N_msr == 2) %>% slice_max(TopK)
# gene <- mcg_df %>%  filter(N_msr == 2) %>% slice_min(TopK)
# gene <- mcg_df %>%  filter(N_msr == max(N_msr)) %>% slice_max(TopK)
# gene <- mcg_df %>%  filter(N_msr == max(N_msr)) %>% slice_min(TopK)
gene <- mcg_df %>%  filter(Symbol == "Spi1")


# filter(mcg_df, Symbol == gene$Symbol)


# For loading/inspect individual gene raw cor vectors across datasets 
raw_pcor_path <- file.path(amat_dir, "Microglia", "Mm_pcor", paste0(gene$Symbol, "_rawcor.tsv"))
raw_scor_path <- file.path(amat_dir, "Microglia", "Mm_scor", paste0(gene$Symbol, "_rawcor.tsv"))


pmat <- load_or_generate_rawcor(
  in_paths = mcg_p_ids,
  out_path = raw_pcor_path,
  genes = pc_mm$Symbol,
  sub_gene = gene$Symbol
)


smat <- load_or_generate_rawcor(
  in_paths = mcg_s_ids,
  out_path = raw_scor_path,
  genes = pc_mm$Symbol,
  sub_gene = gene$Symbol
)


# Difference in raw cors: higher values mean Pcor was higher
diff_mat <- pmat - smat

# Co-msr (how many times were the genes co-measured with given gene across data)
comsr <- rowSums(!is.na(pmat))

# Count of measured coexpressed genes for each experiment (gene coverage)
cvg <- colSums(!is.na(pmat))


# Summary df for given gene


gene_rank <- data.frame(
  Symbol = pc_mm$Symbol,
  Comsr = comsr,
  Agg_pcor = mcg_p[, gene$Symbol],
  Agg_scor = mcg_s[, gene$Symbol],
  Avg_pcor = rowMeans(pmat, na.rm = TRUE),
  Avg_scor = rowMeans(smat, na.rm = TRUE),
  LT0_pcor = rowSums(pmat < 0, na.rm = TRUE),
  LT0_scor = rowSums(smat < 0, na.rm = TRUE),
  Agg_diff = mcg_p[, gene$Symbol] - mcg_s[, gene$Symbol]
)
gene_rank <- filter(gene_rank, Symbol != gene$Symbol)
gene_rank$Avg_diff <- gene_rank$Avg_pcor - gene_rank$Avg_scor
gene_rank$LT0_diff <- gene_rank$LT0_pcor - gene_rank$LT0_scor
gene_rank$Rank_agg_pcor <- rank(-gene_rank$Agg_pcor, ties.method = "min")
gene_rank$Rank_agg_scor <- rank(-gene_rank$Agg_scor, ties.method = "min")



# Checking that has same count of ties
stopifnot(identical(
  sort(table(gene_rank$Agg_pcor), decreasing = TRUE)[[1]],
  sort(table(gene_rank$Agg_scor), decreasing = TRUE)[[1]]
))



gene_rank_sub <- filter(gene_rank, Comsr >= mcg_min_msr)
gene_rank_sub$Rank_agg_pcor <- rank(-gene_rank_sub$Agg_pcor, ties.method = "min")
gene_rank_sub$Rank_agg_scor <- rank(-gene_rank_sub$Agg_scor, ties.method = "min")


summary(gene_rank)
summary(gene_rank_sub)

cor(select_if(gene_rank, is.numeric), use = "pairwise.complete.obs", method = "spearman")
cor(select_if(gene_rank_sub, is.numeric), use = "pairwise.complete.obs", method = "spearman")

# plot_hist(gene_rank, stat_col = "Comsr")
# plot_hist(gene_rank_sub, stat_col = "Comsr")



# Inspecting aggregate scores
p4a <- plot_grid(
  plot_hist(gene_rank, stat_col = "Agg_pcor", title = gene$Symbol),
  plot_hist(gene_rank, stat_col = "Agg_scor", title = gene$Symbol),
  ncol = 1)

p4b <- plot_grid(
  plot_hist(gene_rank_sub, stat_col = "Agg_pcor", title = gene$Symbol),
  plot_hist(gene_rank_sub, stat_col = "Agg_scor", title = gene$Symbol),
  ncol = 1)


p5a <- qplot(gene_rank, xvar = "Agg_pcor", yvar = "Agg_scor", title = gene$Symbol)
p5b <- qplot(gene_rank_sub, xvar = "Agg_pcor", yvar = "Agg_scor", title = gene$Symbol)



# Inspecting average raw cor and the difference (higher values means higher pcor)
p6a <- qplot(gene_rank, xvar = "Avg_pcor", yvar = "Avg_scor", title = gene$Symbol) + xlab("Avg. Pcor") + ylab("Avg. Scor")
p6b <- qplot(gene_rank, xvar = "Comsr", yvar = "Avg_diff", title = gene$Symbol) + ylab("Pcor - Scor") + geom_smooth(method = "loess", colour = "red")
p6c <- qplot(gene_rank_sub, xvar = "Avg_pcor", yvar = "Avg_scor", title = gene$Symbol) + xlab("Avg. Pcor") + ylab("Avg. Scor")
p6d <- qplot(gene_rank_sub, xvar = "Comsr", yvar = "Avg_diff", title = gene$Symbol) + ylab("Pcor - Scor") + geom_smooth(method = "loess", colour = "red")
p6 <- plot_grid(p6a, p6b, p6c, p6d, nrow = 2)


# qplot(gene_rank, xvar = "Comsr", yvar = "Agg_diff", title = gene$Symbol) + ylab("Pcor agg - Scor agg") + geom_smooth(method = "loess", colour = "red")
# qplot(gene_rank_sub, xvar = "Comsr", yvar = "Agg_diff", title = gene$Symbol) + ylab("Pcor agg - Scor agg") + geom_smooth(method = "loess", colour = "red")


# Aggregate versus average
p7a <- qplot(gene_rank_sub, xvar = "Avg_pcor", yvar = "Agg_pcor", title = gene$Symbol) + geom_vline(xintercept = 0, col = "red") + xlab("Avg. Pcor") + ylab("Agg. Pcor")
p7b <- qplot(gene_rank_sub, xvar = "Avg_scor", yvar = "Agg_scor", title = gene$Symbol) + geom_vline(xintercept = 0, col = "red") + xlab("Avg. Scor") + ylab("Agg. Scor")
p7 <- plot_grid(p7a, p7b, nrow = 1)

# Why are there genes with high agg Pcor but low avg Pcor, but not for Scor? 
filter(gene_rank, Agg_pcor > 0.999) %>% slice_min(Avg_pcor, n = 10)
filter(gene_rank, Agg_scor > 0.999) %>% slice_min(Avg_scor, n = 10)



# Inspecting genes with negative cors (tanking aggregate?)


negcor_scatter <- function(pmat, smat, gene1, gene2) {
  
  plot_df <- data.frame(
    Raw_Pcor = pmat[gene2, ],
    Raw_Scor = smat[gene2, ]) %>% 
    mutate(
      Group = factor(case_when(
        Raw_Pcor < 0 & Raw_Scor < 0 ~ 0,
        Raw_Pcor < 0 | Raw_Scor < 0 ~ 1,
        TRUE ~ 2)))
  
  ggplot(plot_df, aes(x = Raw_Pcor, y = Raw_Scor, fill = Group)) +
    geom_point(shape = 21, size = 2.6) +
    geom_hline(yintercept = 0, col = "red") +
    geom_vline(xintercept = 0, col = "red") +
    ggtitle(paste(gene1, gene2, sep = " - ")) +
    scale_fill_manual(values = c("#b2182b", "#878787", "white"), guide = "none") +
    theme_classic() +
    theme(text = element_text(size = 25),
          plot.margin = margin(c(10, 20, 10, 10)))
  
}


filter(gene_rank_sub, Avg_pcor > 0.1) %>% slice_min(Agg_pcor, n = 10)
filter(gene_rank_sub, Avg_scor > 0.1) %>% slice_min(Avg_scor, n = 10)


negcor_scatter(pmat, smat, gene1 = gene$Symbol, gene2 = "Ets2")


plot(density(gene_rank_sub$LT0_pcor), main = gene$Symbol, xlab = "Cor LT0", cex.lab = 1.5)
lines(density(gene_rank_sub$LT0_scor), col = "blue")
plot_hist(gene_rank_sub, stat_col = "LT0_diff", title = paste0(gene$Symbol), xlab = "Pcor LT0 - Scor LT0")

  

stopifnot(identical(
  as.numeric(mean(gene_rank_sub$LT0_diff)), 
  as.numeric(avg_diff_negcor_sub[gene$Symbol])))



# Dist of average difference of count of neg cors, overlaid with gene of interest
plot_hist(mcg_df_sub, 
          stat_col = "Avg_diff_negcor_sub", 
          title = "All measured genes",
          xlab = "Average Pcor LT0 - Scor LT0") +
  geom_vline(xintercept = avg_diff_negcor_sub[gene$Symbol], col = "red")


# Difference in count of negcor versus similarity

plot_grid(plotlist = list(
  
  qplot(mcg_df_sub, xvar = "Avg_diff_negcor_sub", yvar = "TopK_sub") +
  xlab("Average Pcor LT0 - Scor LT0") +
  ylab("Top 200"),
  
  qplot(mcg_df_sub, xvar = "Avg_diff_negcor_sub", yvar = "BottomK_sub") +
  xlab("Average Pcor LT0 - Scor LT0") +
  ylab("Bottom 200"),
  
  qplot(mcg_df_sub, xvar = "Avg_diff_negcor_sub", yvar = "Cor_sub") +
  xlab("Average Pcor LT0 - Scor LT0") +
  ylab("Cor")),
  
  nrow = 1)




plot_grid(plotlist = list(
  qplot(gene_rank_sub, xvar = "LT0_pcor", yvar = "Agg_pcor", title = gene$Symbol),
  qplot(gene_rank_sub, xvar = "LT0_pcor", yvar = "Avg_pcor", title = gene$Symbol),
  qplot(gene_rank_sub, xvar = "LT0_scor", yvar = "Agg_scor", title = gene$Symbol),
  qplot(gene_rank_sub, xvar = "LT0_scor", yvar = "Avg_scor", title = gene$Symbol)
), nrow = 2)
  


# Difference in aggregate versus difference in count of negcors
qplot(gene_rank_sub, xvar = "LT0_diff", yvar = "Agg_diff")

# Big change in aggregate but little change in negcor 
# (pay attention to ranks... delta may be driven by magnitiude)
filter(gene_rank_sub, abs(Agg_diff) > 0.3 & abs(LT0_diff) < 1) 

# Little change in aggregate despite big change in negcor
# Looks like mostly changes in negcor
filter(gene_rank_sub, abs(Agg_diff) < 0.1 & abs(LT0_diff) > 10) 



# Look at genes that were top 200 in each ranking
top_agg_pcor <- gene_rank %>% slice_min(Rank_agg_pcor, n = 200)
top_agg_scor <- gene_rank %>% slice_min(Rank_agg_scor, n = 200)
top_agg_pcor_sub <- gene_rank_sub %>% slice_min(Rank_agg_pcor, n = 200)
top_agg_scor_sub <- gene_rank_sub %>% slice_min(Rank_agg_scor, n = 200)


length(intersect(top_agg_pcor$Symbol, top_agg_scor$Symbol))
mcg_topk[gene$Symbol]
length(intersect(top_agg_pcor_sub$Symbol, top_agg_scor_sub$Symbol))


# No filtering
plot(top_agg_pcor$Rank_agg_scor, main = gene$Symbol, xlab = "Pcor aggregate top 200", ylab = "Scor aggregate rank", cex.lab = 1.5)
abline(h = 200, col = "red")
hist(top_agg_pcor$Comsr, breaks = 100, main = paste(gene$Symbol, "Pcor aggregate Top 200"), xlab = "Count co-measurement")

plot(top_agg_scor$Rank_agg_pcor, main = gene$Symbol, xlab = "Scor aggregate top 200", ylab = "Pcor aggregate rank", cex.lab = 1.5)
abline(h = 200, col = "red")
hist(top_agg_scor$Comsr, breaks = 100, main = paste(gene$Symbol, "Scor aggregate Top 200"), xlab = "Count co-measurement")

plot(density(top_agg_pcor$Comsr))
lines(density(top_agg_scor$Comsr), col = "red")


# Filtering
plot(top_agg_pcor_sub$Rank_agg_scor, main = gene$Symbol, xlab = "Pcor aggregate top 200", ylab = "Scor aggregate rank", cex.lab = 1.5)
abline(h = 200, col = "red")
hist(top_agg_pcor_sub$Comsr, breaks = 100, main = paste(gene$Symbol, "Pcor aggregate Top 200"), xlab = "Count co-measurement")

plot(top_agg_scor_sub$Rank_agg_pcor, main = gene$Symbol, xlab = "Scor aggregate top 200", ylab = "Pcor aggregate rank", cex.lab = 1.5)
abline(h = 200, col = "red")
hist(top_agg_scor_sub$Comsr, breaks = 100, main = paste(gene$Symbol, "Scor aggregate Top 200"), xlab = "Count co-measurement")

plot(density(top_agg_pcor_sub$Comsr))
lines(density(top_agg_scor_sub$Comsr), col = "red")




fit_p <- lm(Agg_pcor ~ Avg_pcor + LT0_pcor + Comsr, data = gene_rank)
fit_p_sub <- lm(Agg_pcor ~ Avg_pcor + LT0_pcor + Comsr, data = gene_rank_sub)

fit_s <- lm(Agg_scor ~ Avg_scor + LT0_scor + Comsr, data = gene_rank)
fit_s_sub <- lm(Agg_scor ~ Avg_scor + LT0_scor + Comsr, data = gene_rank_sub)


# Average by column (do some datasets skew towards higher values for P/S cor?)
avg_diff_exp <- colMeans(diff_mat, na.rm = TRUE)
hist(avg_diff_exp, breaks = 30)
plot(avg_diff_exp, cvg)

