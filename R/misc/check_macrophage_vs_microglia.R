library(tidyverse)
library(pheatmap)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")

mcg <- readRDS("~/mcgdat/Cormats/Hg_pcor/aggregate_cormat_FZ_hg.RDS")
mac <- readRDS("~/mcgdat/Cormats/Macrophage_hg/aggregate_cormat_FZ_macrophage_hg.RDS")

meta_mcg <- read.delim(mcg_meta_path) %>% filter(Species == "Human")
meta_mac <- read.delim(macro_meta_path) %>% filter(Species == "Human")

n_mcg <- n_distinct(meta_mcg$ID)
n_mac <- n_distinct(meta_mac$ID)


msr_df <- data.frame(
  Symbol = rownames(mcg$NA_mat),
  N_msr_mcg = n_mcg - diag(mcg$NA_mat),
  N_msr_mac = n_mac - diag(mac$NA_mat),
  Prop_msr_mcg = 1 - (diag(mcg$NA_mat) / n_mcg),
  Prop_msr_mac = 1 - (diag(mac$NA_mat) / n_mac)
)
msr_df$Prop_msr_avg <- rowMeans(msr_df[, c("Prop_msr_mcg", "Prop_msr_mac")])


min_mcg <- filter(msr_df, Prop_msr_mcg >= 1/3)
min_mac <- filter(msr_df, Prop_msr_mac >= 1/3)
never <- filter(msr_df, Prop_msr_avg == 0)
always <- filter(msr_df, Prop_msr_avg == 1)
mostly <- filter(msr_df, Prop_msr_avg >= 0.9)


keep_genes <- mostly$Symbol

mcg_mat <- mcg$Agg_mat[keep_genes, keep_genes]
mac_mat <- mac$Agg_mat[keep_genes, keep_genes]

diff_mat <- mcg_mat - mac_mat
diff_df <- mat_to_df(diff_mat, value_name = "Diff_FZ")


shuffle_scor <- lapply(1:100, function(x) {
  pair_shuffle_cor(mcg_mat, mac_mat, ncores = ncore)
})


shuffle_topk <- lapply(1:100, function(x) {
  pair_shuffle_topk(mcg_mat, mac_mat, k = 200, ncores = ncore)
})



sim_df <- data.frame(
  Symbol = keep_genes,
  Scor = pair_colwise_cor(mcg_mat, mac_mat, ncores = ncore),
  Topk = pair_colwise_topk(mcg_mat, mac_mat, k = 200, ncores = ncore),
  Bottomk = pair_colwise_topk(-mcg_mat, -mac_mat, k = 200, ncores = ncore)
)





plot(density(sim_df$Scor))
abline(v = median(sim_df$Scor), col = "red")
abline(v = min(sapply(shuffle_scor, median)), col = "darkgrey")
abline(v = max(sapply(shuffle_scor, median)), col = "darkgrey")


plot(density(sim_df$Topk))
abline(v = median(sim_df$Topk), col = "red")
abline(v = min(sapply(shuffle_topk, median)), col = "darkgrey")
abline(v = max(sapply(shuffle_topk, median)), col = "darkgrey")



plot(density(diff_df$Diff_FZ))


# pheatmap(diff_mat, 
#          cluster_rows = TRUE,
#          cluster_cols = TRUE)



check <- "SPI1"


check_df <- data.frame(
  Symbol = keep_genes, 
  Prop_msr_mcg = 1 - (mcg$NA_mat[keep_genes, check] / n_mcg),
  Prop_msr_mac = 1 - (mac$NA_mat[keep_genes, check] / n_mac),
  FZ_mcg = mcg_mat[, check],
  FZ_mac = mac_mat[, check]
)


cor(select_if(check_df, is.numeric), method = "spearman")
plot(check_df$FZ_mcg, check_df$FZ_mac)
