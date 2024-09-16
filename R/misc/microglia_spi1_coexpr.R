library(tidyverse)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

in_dir <- file.path(amat_dir, "Microglia", "Mm_pcor")
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
amat <- fread_to_mat(file.path(in_dir, "aggregate_matrix.tsv"), genes = pc_mm$Symbol)
rank_tf_mm <- readRDS(rank_tf_mm_path)

tf <- "Spi1"

# TODO: Should be handled with meta table of IDs
ids <- list.files(in_dir, full.names = TRUE, pattern = "cormat")


# Construct a matrix of the raw Pcor of Spi1 profiles for each microglia dataset
cmat <- 
  do.call(cbind,
          lapply(ids, fread_to_mat, genes = pc_mm$Symbol, sub_genes = tf))

colnames(cmat) <- ids


# cmat_cp <- cmat


# Pulling individual ranks (all + column)

# rank_list <- lapply(ids, function(x) {
#   
#   cmat <- fread_to_mat(x, genes = pc_mm$Symbol)
#   
#   cmat <- cmat %>%
#     na_to_zero() %>%
#     diag_to_one() %>%
#     upper_to_na()
#   
#   allrank <- allrank_mat(cmat, ties_arg = "min")
#   
#   colrank <- colrank_mat(cmat, ties_arg = "min")
#   
#   
# })




# Count NAs by experiment
exp_na <- sort(apply(cmat, 2, function(x) sum(is.na(x))))
hist(exp_na, breaks = 20)
hist(length(pc_mm$Symbol) - exp_na, breaks = 20, main = "Count co-measured")


# Remove exps that don't measure SPI1
keep <- exp_na != nrow(cmat)
keep <- names(keep[keep])
cmat <- cmat[, keep]


# Count NAs by genes
gene_na <- apply(cmat, 1, function(x) sum(is.na(x)))
comsr <- length(ids) - gene_na
summary(gene_na)
no_comsr <- sum(gene_na == ncol(cmat))
all_comsr <- sum(gene_na == 0)


# Demonstrating coverage

msr_mat <- (!is.na(cmat)) * 1

pheatmap::pheatmap(msr_mat,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color = c("black", "royalblue"),
                   cluster_rows = TRUE, 
                   cluster_cols = FALSE)


# Impute NA cors to 0 and rank within each dataset
cmat_nato0 <- na_to_zero(cmat)
cmat_rank <- colrank_mat(cmat, ties_arg = "min")
cmat_rank_nato0 <- colrank_mat(cmat_nato0, ties_arg = "min")


# Topk of raw pcors across experiments 

topk_mat <- cmat_nato0[setdiff(rownames(cmat_nato0), tf), ]  # prevent self cor
topk_mat <- colwise_topk_intersect(mat = cmat_nato0, k = 200)
topk_df <- mat_to_df(topk_mat)
summary(topk_df$Value)
exp_avg_topk <- rowMeans(topk_mat)


# Plotting topk versus global versus null (requires 04b...)

# null <- unlist(lapply(sim_null_mm, `[[`, "Topk"))
null <- sim_null_mm[[247]]$Topk
global <- sim_tf_mm$Spi1$Topk
mcg <- topk_df$Value

plot_df <- data.frame(
  Group = c(rep("Null", length(null)),
            rep("Global", length(global)),
            rep("Microglia", length(mcg))),
  Topk = c(null, global, mcg)
)
plot_df$Group <- factor(plot_df$Group, levels = unique(plot_df$Group))


ggplot(plot_df, aes(x = Group, y = Topk)) +
  geom_boxplot(width = 0.6, fill = "slategrey") +
  theme_classic() +
  theme(text = element_text(size = 25))


# Aggregate
avg_cor <- round(rowMeans(cmat, na.rm = TRUE), 5)
avg_cor_nato0 <- round(rowMeans(cmat_nato0, na.rm = TRUE), 5)
avg_rank <- rank(rowMeans(cmat_rank, na.rm = TRUE), ties.method = "min")
avg_rank_nato0 <- rank(rowMeans(cmat_rank_nato0, na.rm = TRUE), ties.method = "min")
max_cor <- apply(cmat_nato0, 1, max)
min_cor <- apply(cmat_nato0, 1, min)


# rp_cor <- rank(-rowSums(log(cmat_rank), na.rm = TRUE), ties.method = "min")
# view(data.frame(avg_rank_nato0))


stopifnot(identical(rownames(amat), names(avg_cor)))


all_rank_df <- rank_tf_mm$Spi1 %>%
  dplyr::rename(All_avg_aggr = Avg_aggr_coexpr,
                All_rank_aggr = Rank_aggr_coexpr,
                All_n_comsr = N_comeasured) %>%
  dplyr::select(Symbol, All_avg_aggr, All_rank_aggr, All_n_comsr)


# TODO: max cor, remove redundant, rank avg cor versus avg rank

rank_df <- data.frame(
  Symbol = rownames(amat),
  Microglia_n_comsr = comsr,
  Microglia_aggr = amat[, tf],
  Microglia_aggr_rank = rank(-amat[, tf], ties.method = "min"),
  Microglia_avg_pcor = avg_cor,
  Microglia_avg_pcor_nato0 = avg_cor_nato0,
  Microglia_avg_rank = avg_rank,
  Microglia_avg_rank_nato0 = avg_rank_nato0,
  Microglia_max_cor = max_cor,
  Microglia_min_cor = min_cor) %>% 
  left_join(all_rank_df, by = "Symbol") %>% 
  mutate(Diff = All_rank_aggr - Microglia_aggr_rank,
         Diff2 = All_rank_aggr - Microglia_avg_rank_nato0)



cor(select_if(rank_df, is.numeric), use = "pairwise.complete.obs", method = "spearman")


rank_df %>% 
  mutate(Low_n = Microglia_n_comsr < 5) %>% 
  ggplot(., aes(x = Microglia_aggr, y = All_avg_aggr, fill = Low_n)) +
  geom_point(shape = 21, alpha = 0.4, colour = "white") +
  ylab("All aggregate coexpr") +
  xlab("Microglia aggregate coexpr") +
  scale_fill_manual(values = c("black", "red")) +
  theme_classic() +
  theme(text = element_text(size = 20))




rank_df %>% 
  mutate(Low_n = Microglia_n_comsr < 5) %>% 
  ggplot(., aes(x = Microglia_n_comsr, y = Microglia_aggr, fill = Low_n)) +
  geom_point(shape = 21, alpha = 0.4, colour = "white") +
  xlab("Count of datasets") +
  ylab("Microglia aggregate coexpr") +
  scale_fill_manual(values = c("black", "red")) +
  theme_classic() +
  theme(text = element_text(size = 20))




# Inspect 
gene <- "Limd2"
amat[gene, tf]
cmat[gene, ]
cmat_rank[gene, ]
cmat_rank_nato0[gene, ]


summary(cmat[gene, ])
hist(cmat[gene, ], breaks = 30, xlim = c(-0.5, 0.5), main = paste(gene, "-", tf))



# Human SPI1 with bulk
rank_hg_mm <- readRDS(rank_tf_hg_path)
bulk_cor <- readRDS("~/Robj/HPA_cor_mat_list.RDS")

pcor_gtex <- data.frame(Pcor_GTEX = bulk_cor$GTEX[, "SPI1"]) %>% 
  rownames_to_column("Symbol") %>% 
  mutate(Rank_GTEX = rank(-Pcor_GTEX))

pcor_hpa <- data.frame(Pcor_HPA = bulk_cor$HPA[, "SPI1"]) %>% 
  rownames_to_column("Symbol") %>% 
  mutate(Rank_HPA = rank(-Pcor_HPA))

spi1_rank <- rank_tf_hg$SPI1 %>% 
  left_join(pcor_gtex, by = "Symbol") %>% 
  left_join(pcor_hpa, by = "Symbol")




# Pcor vs Scor
pamat <- fread_to_mat(file.path(file.path(amat_dir, "Microglia", "Mm_pcor"), "aggregate_matrix.tsv"), genes = pc_mm$Symbol)
samat <- fread_to_mat(file.path(file.path(amat_dir, "Microglia", "Mm_scor"), "aggregate_matrix.tsv"), genes = pc_mm$Symbol)
cor_df <- data.frame(Symbol = pc_mm$Symbol, Pcor = pamat[, "Spi1"], Scor = samat[, "Spi1"])
plot(cor_df$Pcor, cor_df$Scor)
cor(cor_df$Pcor, cor_df$Scor, method = "spearman")


pids <- list.files(file.path(amat_dir, "Microglia", "Mm_pcor"), full.names = TRUE, pattern = "cormat")
sids <- list.files(file.path(amat_dir, "Microglia", "Mm_scor"), full.names = TRUE, pattern = "cormat")


pmat <- 
  do.call(cbind,
          lapply(pids, fread_to_mat, genes = pc_mm$Symbol, sub_genes = tf))


smat <- 
  do.call(cbind,
          lapply(sids, fread_to_mat, genes = pc_mm$Symbol, sub_genes = tf))


colnames(pmat) <- colnames(smat) <- ids
