library(tidyverse)
library(pheatmap)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")


pc_df <- read.delim(ref_mm_path, stringsAsFactors = FALSE)


# Paths for the aggregate matrix lists
fz_path <- file.path(cmat_dir_mm, "aggregate_cormat_FZ_mm.RDS")
allrank_path <- file.path(cmat_dir_mm, "aggregate_cormat_allrank_mm.RDS")
allrank_len_path <- file.path(cmat_dir_mm, "aggregate_cormat_allrank_filter_lenient_mm.RDS")
allrank_str_path <- file.path(cmat_dir_mm, "aggregate_cormat_allrank_filter_stringent_mm.RDS")
colrank_path <- file.path(cmat_dir_mm, "aggregate_cormat_colrank_mm.RDS")
colrank_len_path <- file.path(cmat_dir_mm, "aggregate_cormat_colrank_filter_lenient_mm.RDS")
colrank_str_path <- file.path(cmat_dir_mm, "aggregate_cormat_colrank_filter_stringent_mm.RDS")

fz <- readRDS(fz_path)
na <- fz$NA_mat
fz <- fz$Agg_mat 

allrank <- readRDS(allrank_path)$Agg_mat
allrank_len <- readRDS(allrank_len_path)$Agg_mat
allrank_str <- readRDS(allrank_str_path)$Agg_mat

colrank <- readRDS(colrank_path)$Agg_mat
colrank_len <- readRDS(colrank_len_path)$Agg_mat
colrank_str <- readRDS(colrank_str_path)$Agg_mat


genes_all <- rownames(fz)
genes_len <- rownames(allrank_len)
genes_str <- rownames(allrank_str)





fz_vs_allrank <- pair_colwise_topk(fz, allrank, ncores = ncore)
fz_vs_colrank <- pair_colwise_topk(fz, colrank, ncores = ncore)

fz_vs_allrank_len <- pair_colwise_topk(fz[genes_len, genes_len], allrank_len[genes_len, genes_len], ncores = ncore)
fz_vs_colrank_len <- pair_colwise_topk(fz[genes_len, genes_len], colrank_len[genes_len, genes_len], ncores = ncore)

fz_vs_allrank_str <- pair_colwise_topk(fz[genes_str, genes_str], allrank_str[genes_str, genes_str], ncores = ncore)
fz_vs_colrank_str <- pair_colwise_topk(fz[genes_str, genes_str], colrank_str[genes_str, genes_str], ncores = ncore)


plot(density(fz_vs_allrank), col = "red", main = "All genes")
lines(density(fz_vs_colrank), col = "blue")

plot(density(fz_vs_allrank_len), col = "red", "Lenient cut-off")
lines(density(fz_vs_colrank_len), col = "blue")

plot(density(fz_vs_allrank_str), col = "red", "Stringent cut-off")
lines(density(fz_vs_colrank_str), col = "blue")
view(data.frame(fz_vs_allrank_str))



check <- "C1qc"


df_all <- data.frame(
  Symbol = genes_all, 
  N_msr = max(na) - na[genes_all, check],
  FZ = fz[genes_all, check],
  Allrank = allrank[genes_all, check],
  Colrank = colrank[genes_all, check]
)
df_all <- filter(df_all, Symbol != check)




df_len <- data.frame(
  Symbol = genes_len, 
  N_msr = max(na) - na[genes_len, check],
  FZ = fz[genes_len, check],
  Allrank = allrank[genes_len, check],
  Allrank_len = allrank_len[genes_len, check],
  Colrank = colrank[genes_len, check],
  Colrank_len = colrank_len[genes_len, check]
)
df_len <- filter(df_len, Symbol != check)




df_str <- data.frame(
  Symbol = genes_str, 
  N_msr = max(na) - na[genes_str, check],
  FZ = fz[genes_str, check],
  Allrank = allrank[genes_str, check],
  Allrank_str = allrank_str[genes_str, check],
  Colrank = colrank[genes_str, check],
  Colrank_str = colrank_str[genes_str, check]
)
df_str <- filter(df_str, Symbol != check)



cor(select_if(df_all, is.numeric), method = "spearman", use = "pairwise.complete.obs")
plot(df_all$FZ, df_all$Allrank)
plot(df_all$FZ, df_all$Colrank)
plot(df_all$Colrank, df_all$Allrank)


cor(select_if(df_len, is.numeric), method = "spearman", use = "pairwise.complete.obs")
plot(df_len$FZ, df_len$Allrank)
plot(df_len$FZ, df_len$Allrank_len)
plot(df_len$FZ, df_len$Colrank)
plot(df_len$FZ, df_len$Colrank_len)


cor(select_if(df_str, is.numeric), method = "spearman", use = "pairwise.complete.obs")
plot(df_str$FZ, df_str$Allrank)
plot(df_str$FZ, df_str$Allrank_str)
plot(df_str$FZ, df_str$Colrank)
plot(df_str$FZ, df_str$Colrank_str)



# Fbxl6 had the lowest Top K in the stringent gene set. Golgb1 had the highest
# FZ but weak/middling allrank. Ifngr1 had the highest allrank

cor_paths_mm <- list.files(cmat_dir_mm, full.names = TRUE, pattern = "cormat.tsv")

check <- "Fbxl6"
gene1 <- "Golgb1"
gene2 <- "Ifngr1"

cor_l <- lapply(cor_paths_mm, function(x) {
  fread_to_mat(x, genes = pc_df$Symbol, sub_genes = c(gene1, gene2))
})


vec1 <- sapply(cor_l, function(x) x[check, gene1])
vec1 <- vec1[!is.na(vec1)]

vec2 <- sapply(cor_l, function(x) x[check, gene2])
vec2 <- vec2[!is.na(vec2)]

summary(vec1)
summary(vec2)

sum(vec1 < 0)
sum(vec2 < 0)


hist(vec1, breaks = 100, main = paste(check, gene1, sep = "-"), xlab = "Pcor")
hist(vec2, breaks = 100, main = paste(check, gene2, sep = "-"), xlab = "Pcor")
