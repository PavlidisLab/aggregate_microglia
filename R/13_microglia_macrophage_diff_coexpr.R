## This script explores performing differential coexpression between microglia
## and macrophage coexpression in three different ways:
## 1) Rank differences of the respective aggregates. Permutation/shuffle for null.
## 2) Linear model using underlying FZ values from each dataset
## 3) Another linear model that attempts to control for differential expression.
## In addition, it perform diff *expression* between the cell types using limma.
## The focus here is on mouse data, given there is more microglia data
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
library(limma)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")

set.seed(15)

# Gene tables
pc_mm <- read.delim(ref_mm_path)
tfs_mm <- read.delim(tfs_mm_path)

# Dataset meta
mcg_meta <- read.delim(mcg_meta_dedup_path)
macro_meta <- read.delim(macro_meta_dedup_path)

# Lists of gene count measurement summaries
mcg_summ <- readRDS(mcg_count_summ_list_path)
macro_summ <- readRDS(macro_count_summ_list_path)

# Aggregates: focus on FZ
mcg_agg <- readRDS(mcg_fz_mm_path)
macro_agg <- readRDS(macro_fz_mm_path)

# Paths of the individual coexpr matrices
mcg_cmat_paths <- list.files(cmat_dir_mcg_mm, pattern = "_cormat.tsv", full.names = TRUE)
macro_cmat_paths <- list.files(cmat_dir_macro_mm, pattern = "_cormat.tsv", full.names = TRUE)

# Path for list of coexpr mats together (faster to load than individual)
mcg_macro_tf_cmat_l_path <- file.path(data_out_dir, "mcg_macro_cmat_l.RDS")

# Path of TF list of diff coexpr model controlling for expression
tf_fit_l_path <- file.path(data_out_dir, "mcg_macro_diffcoexpr_tf_lmfit.RDS")

# Considering only genes that met measurement filter in micro and macro
keep_genes <- intersect(mcg_summ$Mouse$Filter_genes, macro_summ$Mouse$Filter_genes)
keep_tfs <- intersect(keep_genes, tfs_mm$Symbol)

# Isolated agg matrices with kept genes
mcg_mat <- mcg_agg$Agg_mat[keep_genes, keep_genes]
macro_mat <- macro_agg$Agg_mat[keep_genes, keep_genes]


# This is the TF that will be inspected across diff coexpr approaches
check_tf <- "Maf"

# This gene gets plotted for its relationship with check_tf.
# Maf - Timp2 selected because this was ostensibly a case of diff coexpression
# that was not explained by differential expression
check_gene <- "Timp2"



# Functions
# ------------------------------------------------------------------------------


# Load raw correlation, FZ transform, and NA/inf to 0s

load_fz_mat <- function(path, keep_genes, sub_genes) {
  
  cmat <- fread_to_mat(path, genes = keep_genes, sub_genes = sub_genes)
  cmat[is.na(cmat)] <- 0
  fz <- fisherz(cmat)
  fz[is.infinite(fz)] <- 0 
  
  return(fz)
  
}


# Load all FZ mats into a list

load_all_fz_mat <- function(paths, keep_genes, sub_genes) {
  mat_l <- lapply(paths, load_fz_mat, keep_genes, sub_genes)
  names(mat_l) <- paths
  return(mat_l)
}


# For a given TF, load its profiles across datasets into one matrix

prepare_tf_mat <- function(tf, mcg_l, macro_l) {
  
  mcg_mat <- do.call(cbind, lapply(mcg_l, function(x) x[, tf, drop = FALSE]))
  colnames(mcg_mat) <- paste0("Microglia_", names(mcg_l))
  
  macro_mat <- do.call(cbind, lapply(macro_l, function(x) x[, tf, drop = FALSE]))
  colnames(macro_mat) <- paste0("Macrophage_", names(macro_l))
  
  tf_mat <- cbind(mcg_mat, macro_mat)
  tf_mat <- tf_mat[, which(colSums(tf_mat) != 0)]
  
  return(tf_mat)
}


# Null rank differences by permuting datasets

generate_null_diffrank <- function(check_tf, mcg_l, macro_l, iters = 1000) {
  
  # Ready TF matrix and get the indices of micro/macro profiles
  tf_mat <- prepare_tf_mat(check_tf, mcg_l, macro_l)
  mcg_ix <- which(str_detect(colnames(tf_mat),  "Microglia"))
  macro_ix <- which(str_detect(colnames(tf_mat),  "Macrophage"))
  
  # Iteratively sample datasets and then average, rank, and take difference
  null_diffrank <- mclapply(1:iters, function(x) {
    
    sample_ix <-  sample(1:ncol(tf_mat), ncol(tf_mat), replace = FALSE)
    group1 <- sample_ix[1:length(mcg_ix)]
    group2 <- sample_ix[length(mcg_ix) + 1:length(macro_ix)]
    
    avg_group1 <- rowMeans(tf_mat[, group1])
    rank_group1 <- rank(-avg_group1, ties.method = "min")
    
    avg_group2 <- rowMeans(tf_mat[, group2])
    rank_group2 <- rank(-avg_group2, ties.method = "min")
    
    diff <- rank_group1 - rank_group2
    
  }, mc.cores = ncore)
  
  
  null_diffrank <- do.call(cbind, null_diffrank)
  
  return(null_diffrank)
  
}


# Summarize null rank differences in three ways: 
# Z scores: (Observed difference - mean null difference) / s.d of null
# Fold change: Observed difference / mean null difference
# Empirical pval: Proportion of absolute observed differences less than all absolute nulls

summarize_diffrank <- function(diffrank_df, keep_genes, null_diffrank, iters = 1000) {
  
  summ_null_diffrank <- lapply(keep_genes, function(gene) {
    
    obs <- diffrank_df[gene, "Diff"]
    null <- null_diffrank[gene, ]
    null_mean <- mean(null)
    
    pval <- sum(abs(obs) <= abs(null)) / iters
    z <- (obs - null_mean) / sd(null)
    fc <- obs / null_mean
    
    data.frame(Symbol = gene, 
               Diff = obs,
               Mean_null_diff = null_mean,
               Pval_diff = pval, 
               Z_diff = z, 
               FC_diff = fc)
    
  }) %>% # Join comparison to null with diff rank df
    do.call(rbind, .) %>% 
    left_join(., diffrank_df, by = c("Symbol", "Diff"))
  
}


# Use limma to fit diff coexpr model of cell type status w/o expression covariates

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



# Fit limma diff coexpr model over all TFs

fit_all_model <- function(mcg_l, macro_l, tfs, ncore = 1) {
  
  res_l <- mclapply(tfs, function(tf) {
    tf_mat <- prepare_tf_mat(tf, mcg_l, macro_l)
    res <- fit_model(tf_mat)
  }, mc.cores = ncore)
  
  names(res_l) <- tfs
  return(res_l)
}


# For a given TF, ready a dataframe that includes expression info to perform
# diff coexpr while controlling for expression

prepare_gene_df <- function(tf, gene, mcg_l, macro_l, expr_mat) {
  
  
  fz <- prepare_tf_mat(tf, mcg_l, macro_l)[gene, ]
  cts <- str_extract(names(fz), "Microglia|Macrophage")
  ids <- names(fz)
  
  gene_df <- data.frame(
    FZ = fz, 
    CT = cts,
    ID = ids,
    Expr_TF = expr_mat[tf, ids],
    Expr_gene = expr_mat[gene, ids],
    Expr_diff = expr_mat[tf, ids] - expr_mat[gene, ids],
    Expr_absdiff = abs(expr_mat[tf, ids] - expr_mat[gene, ids]),
    Expr_sum = expr_mat[tf, ids] + expr_mat[gene, ids]
  )
  
  return(gene_df)
  
}


# For the given TF, loop over all genes, generating a dataframe of the TF-gene
# FZs across datasets as their expression levels, then use lm() to fit the 
# FZs as a function of cell type status and expression levels

fit_model_w_expr <- function(tf, genes, mcg_l, macro_l, expr_mat, ncores) {
  
  fit_expr_l <- mclapply(genes, function(gene) {
    
    if (gene == tf) return(NA) # Ignore self cor
    
    # data frame of FZ and expression values across datasets for TR and gene
    gene_df <- prepare_gene_df(tf, gene, mcg_l, macro_l, expr_mat)
    
    # summarize model fit
    fit <- lm(FZ ~ CT + Expr_absdiff + Expr_sum, data = gene_df)
    summ <- summary(fit)
    
    # Extracting coefficients and pvals of cell type and expression covariates
    data.frame(
      TF = tf,
      Gene = gene,
      CT_est = summ$coefficients["CTMicroglia", "Estimate"],
      CT_pval = summ$coefficients["CTMicroglia", "Pr(>|t|)"],
      Expr_absdiff_est = summ$coefficients["Expr_absdiff", "Estimate"],
      Expr_absdiff_pval = summ$coefficients["Expr_absdiff", "Pr(>|t|)"],
      Expr_sum_est = summ$coefficients["Expr_sum", "Estimate"],
      Expr_sum_pval = summ$coefficients["Expr_sum", "Pr(>|t|)"]
    )
    
  }, mc.cores = ncore)
  
  # Bind all gene fits into one df and adjust pvals
  fit_expr_df <- do.call(rbind, fit_expr_l) %>%
    mutate(CT_adj_pval = p.adjust(CT_pval, method = "BH"),
           Expr_absdiff_adj_pval = p.adjust(Expr_absdiff_pval, method = "BH"),
           Expr_sum_adj_pval = p.adjust(Expr_sum_pval, method = "BH"))
  
  return(fit_expr_df)
}



# This is the same approach but does not control for expression levels. This
# was used to ensure that it was giving equivalent answers to the limma model 
# with no expression covariates

fit_model_wo_expr <- function(tf, genes, mcg_l, macro_l, expr_mat, ncores) {
  
  fit_noexpr_l <- mclapply(genes, function(gene) {
    
    if (gene == tf) return(NA) # Ignore self cor
    
    # data frame of FZ and expression values across datasets for TR and gene
    gene_df <- prepare_gene_df(tf, gene, mcg_l, macro_l, expr_mat)
    
    # summarize model fit
    fit <- lm(FZ ~ CT, data = gene_df)
    summ <- summary(fit)
    
    # Extracting coefficients and pvals of cell type
    data.frame(
      TF = check_tf,
      Gene = gene,
      CT_est = summ$coefficients["CTMicroglia", "Estimate"],
      CT_pval = summ$coefficients["CTMicroglia", "Pr(>|t|)"]
    )
    
  }, mc.cores = ncore)
  
  # Bind all gene fits into one df and adjust pvals
  fit_noexpr_df <- do.call(rbind, fit_noexpr_l) %>% 
    mutate(CT_adj_pval = p.adjust(CT_pval, method = "BH"))
  
  return(fit_noexpr_df)
}



# Fit differential coexpression model controlling for expression over all TFs.
# This is very slow, would be good to figure out how to include gene-level 
# covariates into limma...

fit_model_w_expr <- function(tfs, genes, mcg_l, macro_l, expr_mat, ncores) {
  
  
  fit_l <- mclapply(tfs, function(tf) {
    fit_model_w_expr(tf, genes, mcg_l, macro_l, expr_mat, ncores = 1)
  }, mc.cores = ncores)
  names(fit_l) <- tfs
  
  return(fit_l)
}


# Load individual FZ mats for microglia and macrophage into a list
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



# Differential coexpression via rank differences between aggregates. Note that
# here, I am taking the FZ aggregates, then *column* ranking each (i.e., each
# TR profile is ranked individually), and the taking the difference.
# ------------------------------------------------------------------------------


mcg_rank_mat <- colrank_mat(mcg_mat)
macro_rank_mat <- colrank_mat(macro_mat)
diff_mat <- mcg_rank_mat - macro_rank_mat



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



null_diffrank <- generate_null_diffrank(check_tf, mcg_l, macro_l)
diffrank_df <- summarize_diffrank(diffrank_df, keep_genes, null_diffrank)



# Running diff coexpr using limma -- here not accounting for expression
# ------------------------------------------------------------------------------


res_l <- fit_all_model(mcg_l, macro_l, keep_tfs, ncore)


# Count of genes differentially coexpressed w/o controlling for expression
n_sig_wo_expr <- sapply(res_l, function(x) sum(x$adj.P.Val < 0.05, na.rm = TRUE))
prop_sig_wo_expr <- sum(n_sig_wo_expr > 0) / length(n_sig_wo_expr)



# Perform differential *expression* between microglia and macrophages, using
# the averaged/pseudobulked log2 cpm quantile norm'd data
# ------------------------------------------------------------------------------


# No voom since data on hand is log2CPM

expr_mat <- cbind(mcg_summ$Mouse$QN_Avg[keep_genes, ], 
                  macro_summ$Mouse$QN_Avg[keep_genes, ])


expr_cts <- c(rep("Microglia", ncol(mcg_summ$Mouse$QN_Avg)),
              rep("Macrophage", ncol(macro_summ$Mouse$QN_Avg)))


colnames(expr_mat) <- paste0(expr_cts, "_", colnames(expr_mat))

expr_cts <- factor(expr_cts)
expr_design <- model.matrix(~ expr_cts)

expr_fit <- eBayes(lmFit(expr_mat, expr_design))

expr_res <- topTable(expr_fit,
                     coef = "expr_ctsMicroglia",
                     number = Inf,
                     adjust.method = "BH") %>% 
  rownames_to_column(var = "Symbol")




# Are there any significant TF-gene pairs that don't have sig expr changes?


# TFs that have at least 1 sig diff coexpr. gene (not controlling for expression)
diff_tfs_wo_expr <- names(which(n_sig_wo_expr > 0))


# Iterate through diff coexpr pairs, checking whether the TF itself is DE and 
# counting how many of its diff coexpr. genes are also DE

diff_df_wo_expr <- lapply(diff_tfs_wo_expr, function(tf) {
  
  tf_de <- filter(expr_res, Symbol == tf)
  genes_dc <- filter(res_l[[tf]], adj.P.Val < 0.05)$Symbol
  genes_de <- filter(expr_res, Symbol %in% genes_dc)
  n_dc_and_de <- sum(genes_de$adj.P.Val < 0.05)
  
  data.frame(
    TF = tf,
    TF_FC = tf_de$logFC,
    TF_FDR = tf_de$adj.P.Val,
    TF_n_DC = n_sig_wo_expr[tf],
    n_DC_and_DE = n_dc_and_de
  )
  
}) %>% 
  bind_rows()



# Count of TFs that were themselves DE
n_tf_de <- sum(diff_df_wo_expr$TF_FDR < 0.05)


# Looking for genes that are diff coexpressed but not diff. expressed for check tf

dc_not_de <- res_l[[check_tf]] %>% 
  left_join(., expr_res, by = "Symbol", suffix = c("_DC", "_DE")) %>% 
  filter(adj.P.Val_DC < 0.05 & adj.P.Val_DE > 0.05) %>% 
  slice_max(abs(logFC_DC), n = 10)
  


# Adding expression covariates into a differential coexpression model
# NOTE: I couldn't figure out how to build the design that includes feature/gene
# level covariates for limma, so this is done iteratively per TF using lm().
# -----------------------------------------------------------------------------


# Fit over all gene pairs with check tf, controlling for joint expression
fit_expr_df <- fit_model_w_expr(check_tf, keep_genes, mcg_l, macro_l, expr_mat, ncores)

# Fit over all gene pairs with check tf, without controlling for joint expression
fit_noexpr_df <- fit_model_wo_expr(check_tf, keep_genes, mcg_l, macro_l, expr_mat, ncores)


# Checking that limma -expr and lm -expr models are equivalent
# check_noexpr <- left_join(fit_noexpr_df, res_l[[check_tf]], by = c("Gene" = "Symbol"))
# plot(check_noexpr$CT_est, check_noexpr$logFC, xlab = "lm w/o expr", ylab = "limma w/o expr", cex.lab = 1.6)
# plot(check_noexpr$CT_pval, check_noexpr$P.Value, xlab = "lm w/o expr", ylab = "limma w/o expr", cex.lab = 1.6)


# Here I am joining the differential rank, differential expression, and diff
# coexpr models (no expr is using the lm, not limma fit) for check TF.

compare_df <- 
  left_join(fit_expr_df, fit_noexpr_df, by = "Gene", suffix = c("_w_expr", "_wo_expr")) %>% 
  left_join(expr_res, by = c("Gene" = "Symbol")) %>% 
  left_join(diffrank_df, by = c("Gene" = "Symbol")) %>% 
  dplyr::select(-c(TF_w_expr, TF_wo_expr))



cor_compre_df <- cor(select_if(compare_df, is.numeric), use = "pairwise.complete.obs")


# Fit diff coexpr model controlling for expression over all TFs (slow!!)

if (!file.exists(tf_fit_l_path)) {
  
  fit_model_w_expr(tfs = keep_tfs, 
                   genes = keep_genes, 
                   mcg_l = mcg_l, 
                   macro_l = macro_l, 
                   expr_mat = expr_mat, 
                   ncores = ncore)
  
  saveRDS(tf_fit_l, tf_fit_l_path)
  
} else {
  
  tf_fit_l <- readRDS(tf_fit_l_path)

}



# Count of sig genes across all TFs after controlling for expression

n_sig_w_expr <- unlist(
  lapply(tf_fit_l, function(x) sum(x$CT_adj_pval < 0.05, na.rm = TRUE))
)


prop_sig_w_expr <- sum(n_sig_w_expr > 0) / length(n_sig_w_expr)



# Same as before, iterate through TFs that have sig diff coexpr genes (this time
# after controlling for expression), and tally DE/diff coexpr status
diff_tfs_w_expr <- names(which(n_sig_w_expr > 0))


diff_df_w_expr <- lapply(diff_tfs_w_expr, function(tf) {
  
  tf_de <- filter(expr_res, Symbol == tf)
  genes_dc <- filter(tf_fit_l[[tf]], CT_adj_pval < 0.05)$Gene
  genes_de <- filter(expr_res, Symbol %in% genes_dc)
  n_dc_and_de <- sum(genes_de$adj.P.Val < 0.05)
  
  data.frame(
    TF = tf,
    TF_FC = tf_de$logFC,
    TF_FDR = tf_de$adj.P.Val,
    TF_n_DC = n_sig_w_expr[tf],
    n_DC_and_DE = n_dc_and_de
  )
  
}) %>% 
  bind_rows()




# Plots 
# ------------------------------------------------------------------------------


# Histogram of rank differences for check tf
p1 <- ggplot(diffrank_df, aes(x = Scale_diff)) +
  geom_histogram(bins = 100) +
  ylab("Count of genes") +
  xlab("Diff. rank of aggregate coexpr") +
  ggtitle(check_tf) +
  theme_classic() +
  theme(text = element_text(size = 25))


# Volcano plot of DE between micro and macro
p2 <- ggplot(expr_res, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(shape = 21, size = 2.4, alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_classic() +
  theme(text = element_text(size = 25))


# Volcano plot of DE focused on TFs
p3 <- ggplot(diff_df, aes(x = TF_FC, y = -log10(TF_FDR))) +
  geom_point(shape = 21, size = 2.4) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_classic() +
  theme(text = element_text(size = 25))


# Pval hist of diff coexpr (-expr covariates for check tf)
p4a <- ggplot(res_l[[check_tf]], aes(x = P.Value)) +
  geom_histogram(bins = 100) +
  ggtitle(paste0(check_tf, " differential coexpression")) +
  theme_classic() +
  theme(text = element_text(size = 25))
  

# Volcano plot of diff coexpr (-expr covariates for check tf)
p4b <- ggplot(res_l[[check_tf]], aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(shape = 21, size = 2.4, alpha = 0.4) +
  geom_hline(yintercept = -log10(0.05)) +
  ggtitle(paste0(check_tf, " differential coexpression")) +
  theme_classic() +
  theme(text = element_text(size = 25))


# Inspecting coexpr and expression values for an individual TF-gene pair

check_df <- prepare_gene_df(check_tf, check_gene, mcg_l, macro_l, expr_mat)

tf_pval <- filter(expr_res, Symbol == check_tf)$adj.P.Val
gene_pval <- filter(expr_res, Symbol == check_gene)$adj.P.Val


p5a <- ggplot(check_df, aes(x = CT, y = FZ)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1.8, alpha = 0.6, fill = "slategrey", width = 0.05) +
  ggtitle(paste0(check_tf, " - ", check_gene, " coexpr.")) +
  xlab(NULL) +
  theme_classic() +
  theme(text = element_text(size = 20))


p5b <- ggplot(check_df, aes(x = CT, y = Expr_TF)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1.8, alpha = 0.6, fill = "slategrey", width = 0.05) +
  labs(
    y = "Average log2 CPM",
    x = NULL,
    title = paste0(check_tf, " expression"), 
    subtitle = paste0("Adj. pval=", signif(tf_pval, digits = 3))) +
  theme_classic() +
  theme(text = element_text(size = 20))


p5c <- ggplot(check_df, aes(x = CT, y = Expr_gene)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 1.8, alpha = 0.6, fill = "slategrey", width = 0.05) +
  labs(
    y = "Average log2 CPM",
    x = NULL,
    title = paste0(check_gene, " expression"), 
    subtitle = paste0("Adj. pval=", signif(gene_pval, digits = 3))) +
  theme_classic() +
  theme(text = element_text(size = 20))


p5 <- egg::ggarrange(p5a, p5b, p5c, nrow = 1)


# Scatter of check tf-gene significance for fit with and without controlling for expression

p6 <- ggplot(compare_df, aes(x = -log10(CT_adj_pval_w_expr), y = -log10(CT_adj_pval_wo_expr))) +
  geom_point(shape = 21, size = 2.4, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "firebrick") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "firebrick") +
  ylab("-log10 adj. pval for lm (-expr)") +
  xlab("-log10 adj. pval for lm (+expr)") +
  ggtitle(check_tf) +
  theme_classic() +
  theme(text = element_text(size = 20))
