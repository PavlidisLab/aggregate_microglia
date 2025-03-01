library(tidyverse)
library(data.table)
library(aggtools)
library(limma)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")

pc_mm <- read.delim(ref_mm_path)
tfs_mm <- read.delim(tfs_mm_path)

mcg_meta <- read.delim(mcg_meta_dedup_path)
macro_meta <- read.delim(macro_meta_dedup_path)

mcg_summ <- readRDS(mcg_count_summ_list_path)
macro_summ <- readRDS(macro_count_summ_list_path)

mcg_agg <- readRDS("/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/aggregate_cormat_FZ_mm.RDS")
macro_agg <- readRDS("/space/scratch/amorin/aggregate_microglia/Cormats/Macrophage_mm/aggregate_cormat_FZ_macrophage_mm.RDS")

mcg_cmat_paths <- list.files("/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/", pattern = "_cormat.tsv", full.names = TRUE)
macro_cmat_paths <- list.files("/space/scratch/amorin/aggregate_microglia/Cormats/Macrophage_mm/", pattern = "_cormat.tsv", full.names = TRUE)


keep_genes <- intersect(mcg_summ$Mouse$Filter_genes, macro_summ$Mouse$Filter_genes)
keep_tfs <- intersect(keep_genes, tfs_mm$Symbol)

mcg_mat <- mcg_agg$Agg_mat[keep_genes, keep_genes]
macro_mat <- macro_agg$Agg_mat[keep_genes, keep_genes]

mcg_macro_tf_cmat_l_path <- "/space/scratch/amorin/aggregate_microglia/mcg_macro_cmat_l.RDS"



# Functions
# ------------------------------------------------------------------------------


load_fz_mat <- function(path, keep_genes, sub_genes) {
  
  cmat <- fread_to_mat(path, genes = keep_genes, sub_genes = sub_genes)
  cmat[is.na(cmat)] <- 0
  fz <- fisherz(cmat)
  fz[is.infinite(fz)] <- 0 
  
  return(fz)
  
}



load_all_fz_mat <- function(paths, keep_genes, sub_genes) {
  mat_l <- lapply(paths, load_fz_mat, keep_genes, sub_genes)
  names(mat_l) <- paths
  return(mat_l)
}




# prepare_tf_mat <- function(mcg_l, macro_l, tf) {
#   
#   mcg_mat <- do.call(cbind, lapply(mcg_l, function(x) x[, tf, drop = FALSE]))
#   macro_mat <- do.call(cbind, lapply(macro_l, function(x) x[, tf, drop = FALSE]))
#   tf_mat <- cbind(mcg_mat, macro_mat)
#   colnames(tf_mat) <- c(rep("Microglia", ncol(mcg_mat)), rep("Macrophage", ncol(macro_mat)))
#   tf_mat <- tf_mat[, which(colSums(tf_mat) != 0)]
#   
#   return(tf_mat)
# }



prepare_tf_mat <- function(tf, mcg_l, macro_l) {
  
  mcg_mat <- do.call(cbind, lapply(mcg_l, function(x) x[, tf, drop = FALSE]))
  colnames(mcg_mat) <- paste0("Microglia_", names(mcg_l))
  
  macro_mat <- do.call(cbind, lapply(macro_l, function(x) x[, tf, drop = FALSE]))
  colnames(macro_mat) <- paste0("Macrophage_", names(macro_l))
  
  tf_mat <- cbind(mcg_mat, macro_mat)
  tf_mat <- tf_mat[, which(colSums(tf_mat) != 0)]
  
  return(tf_mat)
}



fit_model <- function(tf_mat) {
  
  cts <- str_extract(colnames(tf_mat), "Microglia|Macrophage")
  cts <- factor(cts, levels = rev(c("Microglia", "Macrophage")))
  design <- model.matrix(~ cts)
  fit <- eBayes(lmFit(tf_mat, design))
  
  res <- topTable(fit,
                  coef = "ctsMicroglia",
                  number = Inf,
                  adjust.method = "BH") %>%
    rownames_to_column(var = "Symbol")
  
  return(res)
}



fit_all_model <- function(mcg_l, macro_l, tfs, ncore = 1) {
  
  res_l <- mclapply(tfs, function(tf) {
    tf_mat <- prepare_tf_mat(tf, mcg_l, macro_l)
    res <- fit_model(tf_mat)
  }, mc.cores = ncore)
  
  names(res_l) <- tfs
  return(res_l)
}



# Load individual FZ mats for microglia and macrophage
# ------------------------------------------------------------------------------

if (!file.exists(mcg_macro_tf_cmat_l_path)) {
  
  mcg_l <- load_all_fz_mat(mcg_cmat_paths, keep_genes, keep_tfs)
  macro_l <- load_all_fz_mat(macro_cmat_paths, keep_genes, keep_tfs)
  saveRDS(list(mcg = mcg_l, macro = macro_l), file = mcg_macro_tf_cmat_l_path)
  
} else {
  
  dat <- readRDS(mcg_macro_tf_cmat_l_path)
  mcg_l <- dat$mcg
  macro_l <- dat$macro
  rm(dat)
}



# Inspecting rank differences of aggregate
# ------------------------------------------------------------------------------


mcg_rank_mat <- colrank_mat(mcg_mat)
macro_rank_mat <- colrank_mat(macro_mat)
diff_mat <- mcg_rank_mat - macro_rank_mat


check_tf <- "Mef2c"


# Scaling rank difference to -1 and 1
# https://stackoverflow.com/questions/66602071/r-scaling-a-vector-between-1-and-1
rescale_minMax <- function(x) {
  2 * (1 - (x - max(x)) / (min(x) - max(x))) - 1 
}


diffrank_df <- data.frame(
  Symbol = keep_genes,
  Agg_mcg = mcg_mat[, check_tf],
  Agg_macro = macro_mat[, check_tf],
  Rank_mcg = mcg_rank_mat[, check_tf],
  Rank_macro = macro_rank_mat[, check_tf],
  Diff = diff_mat[, check_tf]
)
diffrank_df$Scale_diff <- rescale_minMax(diffrank_df$Diff)




# ggplot(diffrank_df, aes(x = Diff)) +
ggplot(diffrank_df, aes(x = Scale_diff)) +
  geom_histogram(bins = 100) +
  ylab("Count of genes") +
  xlab("Diff. rank of aggregate coexpr") +
  ggtitle(check_tf) +
  theme_classic() +
  theme(text = element_text(size = 25))


# Null rank differences by permuting datasets


tf_mat <- prepare_tf_mat(check_tf, mcg_l, macro_l)

mcg_ix <- which(str_detect(colnames(tf_mat),  "Microglia"))
macro_ix <- which(str_detect(colnames(tf_mat),  "Macrophage"))


set.seed(15)


null_diff <- mclapply(1:1000, function(x) {
  
  sample_ix <-  sample(1:ncol(tf_mat), ncol(tf_mat), replace = FALSE)
  group1 <- sample_ix[1:length(mcg_ix)]
  group2 <- sample_ix[length(mcg_ix) + 1:length(macro_ix)]
  
  avg_group1 <- rowMeans(tf_mat[, group1])
  rank_group1 <- rank(-avg_group1, ties.method = "min")
  
  avg_group2 <- rowMeans(tf_mat[, group2])
  rank_group2 <- rank(-avg_group2, ties.method = "min")
  
  diff <- rank_group1 - rank_group2

  
}, mc.cores = ncore)


null_diff <- do.call(cbind, null_diff)
null_diff_scale <- apply(null_diff, 2, rescale_minMax)



# hist(null_diff_scale["Csmd3", ], breaks = 1000, xlim = c(-1, 1))
# abline(v = diffrank_df["Csmd3", "Scale_diff"], col = "red")
# hist(null_diff["Csmd3", ], breaks = 1000, xlim = c(min(diffrank_df$Diff), max(diffrank_df$Diff)))
# abline(v = diffrank_df["Csmd3", "Diff"], col = "red")
# 
# null_diff_mean <- rowMeans(null_diff)
# null_diff_scale_mean <- rowMeans(null_diff_scale)
# hist(null_diff_mean, breaks = 1000)
# hist(null_diff_scale_mean, breaks = 1000)
# plot(null_diff_mean, diffrank_df$Scale_diff)


summ_null <- lapply(keep_genes, function(gene) {
  
  obs <- diffrank_df[gene, "Diff"]
  null <- null_diff[gene, ]
  null_mean <- mean(null)
  
  pval <- sum(abs(obs) <= abs(null)) / 1000
  z <- (obs - null_mean) / sd(null)
  fc <- obs / null_mean
  
  data.frame(Symbol = gene, 
             Obs_diff = obs,
             Mean_null_diff = null_mean,
             Pval = pval, 
             Z = z, 
             FC = fc)
  
})


summ_null <- do.call(rbind, summ_null)



summ_null_scale <- lapply(keep_genes, function(gene) {
  
  obs <- diffrank_df[gene, "Scale_diff"]
  null <- null_diff_scale[gene, ]
  null_mean <- mean(null)
  
  pval <- sum(abs(obs) <= abs(null)) / 1000
  z <- (obs - null_mean) / sd(null)
  fc <- obs / null_mean
  
  data.frame(Symbol = gene, 
             Obs_diff = obs,
             Mean_null_diff = null_mean,
             Pval = pval, 
             Z = z, 
             FC = fc)
  
})

summ_null_scale <- do.call(rbind, summ_null_scale)



# hist(summ_null$Pval, breaks = 1000)
# hist(summ_null$Z, breaks = 1000)
# hist(summ_null$FC, breaks = 1000)
# plot(summ_null$Obs_diff, summ_null$Z)
# 
# hist(summ_null_scale$Pval, breaks = 1000)
# hist(summ_null_scale$Z, breaks = 1000)
# hist(summ_null_scale$FC, breaks = 1000)
# plot(summ_null_scale$Obs_diff, summ_null_scale$Z)


# TODO: figure out why this isn't identical
# plot(rowMeans(tf_mat[keep_genes, mcg_ix]), mcg_mat[keep_genes, check_tf])
# 
# 
# plot(
#   rank(-rowMeans(tf_mat[keep_genes, mcg_ix]), ties.method = "min"),
#   diffrank_df[keep_genes, "Agg_mcg"]
# )




# Running diff coexpr using limma -- here no accountng for expression
# ------------------------------------------------------------------------------


res_l <- fit_all_model(mcg_l, macro_l, keep_tfs, ncore)



# Count of genes differentially coexpressed
n_sig <- sapply(res_l, function(x) sum(x$adj.P.Val < 0.05, na.rm = TRUE))
sum(n_sig > 0) / length(n_sig)
summary(n_sig)
hist(n_sig, breaks = 100)




# Are there any significant TF-gene pairs that don't have sig expr changes
# ------------------------------------------------------------------------------


# First perform DE between microglia and macrophage pseudobulks
# No voom since data on hand is log2CPM

expr_mat <- cbind(mcg_summ$Mouse$QN_Avg[keep_genes, ], 
                  macro_summ$Mouse$QN_Avg[keep_genes, ])


colnames(expr_mat) <- c(rep("Microglia", ncol(mcg_summ$Mouse$QN_Avg[keep_genes, ])),
                        rep("Macrophage", ncol(macro_summ$Mouse$QN_Avg[keep_genes, ])))

cts <- factor(colnames(expr_mat))

expr_design <- model.matrix(~ cts)

expr_v <- voom(expr_mat, expr_design)

expr_fit <- eBayes(lmFit(expr_mat, expr_design))
# expr_fit <- eBayes(lmFit(expr_v, expr_design))


expr_res <- topTable(expr_fit,
                     coef = "ctsMicroglia",
                     number = Inf,
                     adjust.method = "BH") %>% 
  rownames_to_column(var = "Symbol")



ggplot(expr_res, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(shape = 21, size = 2.4, alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_classic() +
  theme(text = element_text(size = 25))



# Iterate through diff coexpr pairs, checking if TF and genes are DE

diff_tfs <- names(which(n_sig > 0))


diff_l <- lapply(diff_tfs, function(tf) {
  
  tf_de <- filter(expr_res, Symbol == tf)
  genes_dc <- filter(res_l[[tf]], adj.P.Val < 0.05)$Symbol
  genes_de <- filter(expr_res, Symbol %in% genes_dc)
  n_dc_and_de <- sum(genes_de$adj.P.Val < 0.05)
  
  data.frame(
    TF = tf,
    TF_FC = tf_de$logFC,
    TF_FDR = tf_de$adj.P.Val,
    TF_n_DC = n_sig[tf],
    n_DC_and_DE = n_dc_and_de
  )
  
})


diff_df <- do.call(rbind, diff_l)


n_tf_de <- sum(diff_df$TF_FDR < 0.05)


ggplot(diff_df, aes(x = TF_FC, y = -log10(TF_FDR))) +
  geom_point(shape = 21, size = 2.4) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_classic() +
  theme(text = element_text(size = 25))






# Manual inspection and plotting of individual TFs
# ------------------------------------------------------------------------------


# Pval hist

ggplot(res_l[[check_tf]], aes(x = P.Value)) +
  geom_histogram(bins = 100) +
  ggtitle(paste0(check_tf, " differential coexpression")) +
  theme_classic() +
  theme(text = element_text(size = 25))
  

# Volcano plot
ggplot(res_l[[check_tf]], aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(shape = 21, size = 2.4, alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05)) +
  ggtitle(paste0(check_tf, " differential coexpression")) +
  theme_classic() +
  theme(text = element_text(size = 25))


# Inpsecting values for an individual TF-gene pair
check_gene <- "Dab2"


check_fz <- prepare_tf_mat(mcg_l, macro_l, check_tf)[check_gene, ]
check_fz <- data.frame(FZ = check_fz, CT = names(check_fz))


expr_df <- data.frame(
  TF = expr_mat[check_tf, ],
  Gene = expr_mat[check_gene, ],
  CT = colnames(expr_mat)
)


tf_pval <- filter(expr_res, Symbol == check_tf)$adj.P.Val
# suppressWarnings(wilcox.test(expr_df$TF ~ expr_df$CT))

gene_pval <- filter(expr_res, Symbol == check_gene)$adj.P.Val
# suppressWarnings(wilcox.test(expr_df$Gene ~ expr_df$CT))



pxa <- ggplot(check_fz, aes(x = CT, y = FZ)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1.8, alpha = 0.6, fill = "slategrey", width = 0.05) +
  ggtitle(paste0(check_tf, " - ", check_gene, " coexpr.")) +
  xlab(NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))



pxb <- ggplot(expr_df, aes(x = CT, y = TF)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1.8, alpha = 0.6, fill = "slategrey", width = 0.05) +
  labs(
    y = "Average log2 CPM",
    x = NULL,
    title = paste0(check_tf, " expression"), 
    subtitle = paste0("Adj. pval=", signif(tf_pval, digits = 3))) +
  theme_classic() +
  theme(text = element_text(size = 20))



pxc <- ggplot(expr_df, aes(x = CT, y = Gene)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1.8, alpha = 0.6, fill = "slategrey", width = 0.05) +
  labs(
    y = "Average log2 CPM",
    x = NULL,
    title = paste0(check_gene, " expression"), 
    subtitle = paste0("Adj. pval=", signif(gene_pval, digits = 3))) +
  theme_classic() +
  theme(text = element_text(size = 20))


egg::ggarrange(pxa, pxb, pxc, nrow = 1)





# Looking for genes that are diff coexpressed but not diff. expressed


res_df <- left_join(res_l[[check_tf]], expr_res, 
                    by = "Symbol", suffix = c("_DC", "_DE"))


res_df %>% 
  filter(adj.P.Val_DC < 0.05 & adj.P.Val_DE > 0.05) %>% 
  slice_max(logFC_DC, n = 10)



# Checking adding difference in expression to diff coexpr model
# -----------------------------------------------------------------------------


names(mcg_l) <- str_extract(names(mcg_l), paste(mcg_meta$ID, collapse = "|"))
names(macro_l) <- str_extract(names(macro_l), paste(macro_meta$ID, collapse = "|"))



tf_mat <- prepare_tf_mat2(check_tf, mcg_l, macro_l)
# tf_mat[1:5, 1:5]


ids_mcg <- colnames(tf_mat) %>% 
  str_extract(., "Microglia_.*") %>%
  na.omit() %>% 
  str_replace(., "^.*_", "")


ids_macro <- colnames(tf_mat) %>% 
  str_extract(., "Macrophage_.*") %>%
  na.omit() %>% 
  str_replace(., "^.*_", "")



expr_diff_mcg <- sweep(mcg_summ$Mouse$QN_Avg[keep_genes, ids_mcg], 2, 
                       mcg_summ$Mouse$QN_Avg[check_tf, ids_mcg], "-")


expr_diff_macro <- sweep(macro_summ$Mouse$QN_Avg[keep_genes, ids_macro], 2, 
                         macro_summ$Mouse$QN_Avg[check_tf, ids_macro], "-")


abs_expr_diff_mcg <- abs(expr_diff_mcg)
abs_expr_diff_macro <- abs(expr_diff_macro)



expr_add_mcg <- sweep(mcg_summ$Mouse$QN_Avg[keep_genes, ids_mcg], 2, 
                       mcg_summ$Mouse$QN_Avg[check_tf, ids_mcg], "+")


expr_add_macro <- sweep(macro_summ$Mouse$QN_Avg[keep_genes, ids_macro], 2, 
                         macro_summ$Mouse$QN_Avg[check_tf, ids_macro], "+")



stopifnot(identical(
  expr_diff_mcg[, "GSE123335"],
  mcg_summ$Mouse$QN_Avg[keep_genes, "GSE123335"] - mcg_summ$Mouse$QN_Avg[check_tf, "GSE123335"]
))


stopifnot(identical(
  abs_expr_diff_macro[, "Kozareva2020"],
  abs(macro_summ$Mouse$QN_Avg[keep_genes, "Kozareva2020"] - macro_summ$Mouse$QN_Avg[check_tf, "Kozareva2020"])
))



stopifnot(identical(
  expr_add_mcg[, "GSE102827"],
  mcg_summ$Mouse$QN_Avg[keep_genes, "GSE102827"] + mcg_summ$Mouse$QN_Avg[check_tf, "GSE102827"]
))



tf_df <- data.frame(
  TF = rep(check_tf, nrow(tf_mat) * ncol(tf_mat)),
  Gene = rep(rownames(tf_mat), each = ncol(tf_mat)),
  Dataset = rep(colnames(tf_mat), by = ncol(tf_mat)),
  FZ = as.numeric(t(tf_mat)),
  Expr_diff = c(as.numeric(t(expr_diff_mcg)), as.numeric(t(expr_diff_macro))),
  Expr_absdiff = c(as.numeric(t(abs_expr_diff_mcg)), as.numeric(t(abs_expr_diff_macro))),
  Expr_add = c(as.numeric(t(expr_add_mcg)), as.numeric(t(expr_add_macro)))
)

tf_df$Cell_type <- str_replace(tf_df$Dataset, "_.*", "")



# TODO: can this work for limma?
# tf_mat_long <- matrix(tf_df$FZ, ncol = nrow(tf_df), nrow = 1)
# colnames(tf_mat_long) <- tf_df$Gene
# design <- model.matrix( ~ Cell_type + Abs_diff, data = tf_df)
# fit <- eBayes(lmFit(tf_mat_long, design))
# res <- topTable(fit, coef = "Cell_typeMicroglia", number = Inf)


fit <- lm(FZ ~ Cell_type + Expr_absdiff + Expr_add, data = tf_df, subset = Gene == "Lgals3")
summary(fit)



fit_l <- mclapply(keep_genes, function(gene) {
  
  fit <- lm(FZ ~ Cell_type + Expr_absdiff + Expr_add, 
            data = tf_df, subset = Gene == gene)
  
  summ <- summary(fit)
  
  data.frame(
    TF = check_tf,
    Gene = gene,
    CT_est = summ$coefficients["Cell_typeMicroglia", "Estimate"],
    CT_pval = summ$coefficients["Cell_typeMicroglia", "Pr(>|t|)"],
    Expr_absdiff_est = summ$coefficients["Expr_absdiff", "Estimate"],
    Expr_absdiff_pval = summ$coefficients["Expr_absdiff", "Pr(>|t|)"],
    Expr_add_est = summ$coefficients["Expr_add", "Estimate"],
    Expr_add_pval = summ$coefficients["Expr_add", "Pr(>|t|)"]
  )
  
}, mc.cores = ncore)



fit_l2 <- mclapply(keep_genes, function(gene) {
  
  fit <- lm(FZ ~ Cell_type, data = tf_df, subset = Gene == gene)
  summ <- summary(fit)
  
  
  data.frame(
    TF = check_tf,
    Gene = gene,
    CT_est = summ$coefficients["Cell_typeMicroglia", "Estimate"],
    CT_pval = summ$coefficients["Cell_typeMicroglia", "Pr(>|t|)"]
  )
  
}, mc.cores = ncore)




fit_df <- do.call(rbind, fit_l) %>% 
   mutate(CT_adj_pval = p.adjust(CT_pval, method = "BH"),
         Expr_absdiff_adj_pval = p.adjust(Expr_absdiff_pval, method = "BH"))



fit_df2 <- do.call(rbind, fit_l2) %>% 
  mutate(CT_adj_pval = p.adjust(CT_pval, method = "BH"))



compare_df <- 
  left_join(fit_df, fit_df2, by = "Gene", suffix = c("_w_expr", "_wo_expr")) %>% 
  left_join(res_l[[check_tf]], by = c("Gene" = "Symbol")) %>% 
  left_join(diffrank_df, by = c("Gene" = "Symbol")) %>% 
  left_join(summ_null, by = c("Gene" = "Symbol"))


cor(select_if(compare_df, is.numeric), use = "pairwise.complete.obs")


plot(compare_df$CT_est_wo_expr, compare_df$logFC, xlab = "lm w/o expr", ylab = "limma w/o expr", cex.lab = 1.6)
plot(compare_df$CT_est_w_expr, compare_df$CT_est_wo_expr, xlab = "lm w/ expr", ylab = "lm w/o expr", cex.lab = 1.6)
plot(-log10(compare_df$CT_adj_pval_w_expr), -log10(compare_df$adj.P.Val), xlab = "lm w/o expr", ylab = "limma w/o expr", cex.lab = 1.6)


ggplot(compare_df, aes(x = -log10(CT_adj_pval_w_expr), y = -log10(CT_adj_pval_wo_expr))) +
  geom_point(shape = 21, size = 2.4, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "firebrick") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "firebrick") +
  ylab("-log10 adj. pval for lm (-expr)") +
  xlab("-log10 adj. pval for lm (+expr)") +
  theme_classic() +
  theme(text = element_text(size = 20))



# compare_df %>% 
#   mutate(
#     Score_w_expr = -log10(CT_adj_pval_w_expr),
#     Score_wo_expr = -log10(CT_adj_pval_wo_expr),
#     Diff = Score_w_expr - Score_wo_expr) %>% 
#   dplyr::select(Gene, Score_w_expr, Score_wo_expr, Diff) %>% 
#   view()



compare_df %>% 
  slice_max(abs(Diff), n = 30) %>% 
  arrange(CT_adj_pval_w_expr)



tt1 <- fit_df %>% 
  left_join(diffrank_df, by = c("Gene" = "Symbol")) %>% 
  left_join(summ_null, by = c("Gene" = "Symbol")) %>% 
  left_join(expr_res, by = c("Gene" = "Symbol"))


plot(tt1$CT_est, tt1$logFC)
plot(tt1$Z, tt1$logFC)


ggplot(tt1, aes(x = logFC, y = Z)) +
  geom_point(shape = 21, size = 2.4, fill = "slategrey", colour = "slategrey", alpha = 0.8) +
  geom_hline(yintercept = -2, linetype = "dashed", colour = "firebrick") +
  geom_hline(yintercept = 2, linetype = "dashed", colour = "firebrick") +
  geom_vline(xintercept = -2, linetype = "dashed", colour = "firebrick") +
  geom_vline(xintercept = 2, linetype = "dashed", colour = "firebrick") +
  ylab("Differential coexpression Z-score") +
  xlab("Differential expression logFC") +
  ggtitle(check_tf) +
  theme_classic() +
  theme(text = element_text(size = 25))




# Model over all TFs



ready_tf_df <- function(tf, keep_genes, tf_mat, mcg_summ, macro_summ) {
  
  
  ids_mcg <- colnames(tf_mat) %>% 
    str_extract(., "Microglia_.*") %>%
    na.omit() %>% 
    str_replace(., "^.*_", "")
  
  
  ids_macro <- colnames(tf_mat) %>% 
    str_extract(., "Macrophage_.*") %>%
    na.omit() %>% 
    str_replace(., "^.*_", "")
  
  
  
  expr_diff_mcg <- sweep(mcg_summ$Mouse$QN_Avg[keep_genes, ids_mcg], 2, 
                         mcg_summ$Mouse$QN_Avg[tf, ids_mcg], "-")
  
  
  expr_diff_macro <- sweep(macro_summ$Mouse$QN_Avg[keep_genes, ids_macro], 2, 
                           macro_summ$Mouse$QN_Avg[tf, ids_macro], "-")
  
  
  abs_expr_diff_mcg <- abs(expr_diff_mcg)
  abs_expr_diff_macro <- abs(expr_diff_macro)
  
  
  
  expr_add_mcg <- sweep(mcg_summ$Mouse$QN_Avg[keep_genes, ids_mcg], 2, 
                        mcg_summ$Mouse$QN_Avg[tf, ids_mcg], "+")
  
  
  expr_add_macro <- sweep(macro_summ$Mouse$QN_Avg[keep_genes, ids_macro], 2, 
                          macro_summ$Mouse$QN_Avg[tf, ids_macro], "+")
  
  
  tf_df <- data.frame(
    TF = rep(tf, nrow(tf_mat) * ncol(tf_mat)),
    Gene = rep(rownames(tf_mat), each = ncol(tf_mat)),
    Dataset = rep(colnames(tf_mat), by = ncol(tf_mat)),
    FZ = as.numeric(t(tf_mat)),
    Expr_diff = c(as.numeric(t(expr_diff_mcg)), as.numeric(t(expr_diff_macro))),
    Expr_absdiff = c(as.numeric(t(abs_expr_diff_mcg)), as.numeric(t(abs_expr_diff_macro))),
    Expr_add = c(as.numeric(t(expr_add_mcg)), as.numeric(t(expr_add_macro)))
  )
  
  tf_df$Cell_type <- str_replace(tf_df$Dataset, "_.*", "")
  
  return(tf_df)
}





fit_diffcoexpr_all_tfs <- function(genes, tfs, mcg_l, macro_l, mcg_summ, macro_summ) {
  
  tf_l <- mclapply(tfs, function(tf) {
    
    message(paste(tf, Sys.time()))
    tf_mat <- prepare_tf_mat2(tf, mcg_l, macro_l)
    tf_df <- ready_tf_df(tf, genes, tf_mat, mcg_summ, macro_summ)
    
    fit_l <- lapply(keep_genes, function(gene) {
      
      fit <- lm(FZ ~ Cell_type + Expr_absdiff + Expr_add, 
                data = tf_df, subset = Gene == gene)
      
      summ <- summary(fit)
      
      data.frame(
        TF = tf,
        Gene = gene,
        CT_est = summ$coefficients["Cell_typeMicroglia", "Estimate"],
        CT_pval = summ$coefficients["Cell_typeMicroglia", "Pr(>|t|)"],
        Expr_absdiff_est = summ$coefficients["Expr_absdiff", "Estimate"],
        Expr_absdiff_pval = summ$coefficients["Expr_absdiff", "Pr(>|t|)"],
        Expr_add_est = summ$coefficients["Expr_add", "Estimate"],
        Expr_add_pval = summ$coefficients["Expr_add", "Pr(>|t|)"]
      )
      
    })
    
    fit_df <- do.call(rbind, fit_l) %>% 
      mutate(CT_adj_pval = p.adjust(CT_pval, method = "BH"))
    
  }, mc.cores = ncore)
  
  names(tf_l) <- tfs
  
  return(tf_l)
}


tf_fit_l <- fit_diffcoexpr_all_tfs(genes = keep_genes, 
                                   tfs = keep_tfs, 
                                   mcg_l = mcg_l,
                                   macro_l = macro_l, 
                                   mcg_summ = mcg_summ,
                                   macro_summ = macro_summ)



saveRDS(tf_fit_l, "/space/scratch/amorin/aggregate_microglia/mcg_macro_diffcoexpr_tf_lmfit.RDS")



n_sig2 <- sapply(tf_fit_l, function(x) sum(x$CT_adj_pval < 0.05, na.rm = TRUE))
sum(n_sig2 > 0) / length(n_sig2)
summary(n_sig2)
hist(n_sig2, breaks = 100)

plot(n_sig, n_sig2)


diff_tfs2 <- names(which(n_sig2 > 0))


diff_l2 <- lapply(diff_tfs2, function(tf) {
  
  tf_de <- filter(expr_res, Symbol == tf)
  genes_dc <- filter(tf_fit_l[[tf]], CT_adj_pval < 0.05)$Gene
  genes_de <- filter(expr_res, Symbol %in% genes_dc)
  n_dc_and_de <- sum(genes_de$adj.P.Val < 0.05)
  
  data.frame(
    TF = tf,
    TF_FC = tf_de$logFC,
    TF_FDR = tf_de$adj.P.Val,
    TF_n_DC = n_sig2[tf],
    n_DC_and_DE = n_dc_and_de
  )
  
})


diff_df2 <- do.call(rbind, diff_l2)


n_tf_de2 <- sum(diff_df2$TF_FDR < 0.05)


# Attempting to visualize FZ ~ Cell_type on residuals
tf <- "Mef2c"
gene <- "Dab2"
tf_mat <- prepare_tf_mat2(tf, mcg_l, macro_l)
tf_df <- ready_tf_df(tf, keep_genes, tf_mat, mcg_summ, macro_summ)
tf_df <- filter(tf_df, Gene == gene)


boxplot(tf_df$Expr_diff ~ tf_df$Cell_type)
boxplot(tf_df$Expr_absdiff ~ tf_df$Cell_type)
boxplot(tf_df$Expr_add ~ tf_df$Cell_type)
boxplot(tf_df$FZ ~ tf_df$Cell_type)


expr_df$Diff <- expr_df$TF - expr_df$Gene
expr_df$Abs_diff <- abs(expr_df$Diff)
expr_df$Add <- expr_df$TF + expr_df$Gene

boxplot(expr_df$Diff ~ expr_df$CT)
boxplot(expr_df$Abs_diff ~ expr_df$CT)
boxplot(expr_df$Add ~ expr_df$CT)



fit <- lm(FZ ~ Cell_type + Expr_absdiff + Expr_add, data = tf_df)

residuals_expr_corrected <- resid(fit) + coef(fit)["(Intercept)"] + (coef(fit)["Cell_typeMicroglia"] * as.numeric(tf_df$Cell_type == "Microglia"))

tf_df$FZ_corrected <- residuals_expr_corrected

boxplot(FZ ~ Cell_type, data = tf_df)
boxplot(FZ_corrected ~ Cell_type, data = tf_df)
