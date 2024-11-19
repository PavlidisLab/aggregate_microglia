## 

library(tidyverse)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")

# Load genes/TFs
tf_hg <- read.delim("/home/amorin/Data/Metadata/TFs_human.tsv")
pc_hg <- read.delim(ref_hg_path)

# Load a single count matrix
id <- "GSE180928"
dat_l <- readRDS(mcg_dat_path)
mat <- dat_l[[id]]$Mat
rm(dat_l)
gc()

# Load corresponding raw pcor mat
cormat <- fread_to_mat("/space/scratch/amorin/aggregate_microglia/Cormats/Hg_pcor/GSE180928_cormat.tsv", genes = pc_hg$Symbol)

# Only keeping measured genes/TFs
keep_genes <- names(which(rowSums(mat > 0) >= 20))
keep_tfs <- intersect(tf_hg$Symbol, keep_genes)
mat <- as.matrix(mat[keep_genes, ])
cormat <- cormat[keep_genes, keep_tfs]


# Save filtered count matrix
mat_path <- "/space/scratch/amorin/R_objects/GSE180928_mcg_filt.tsv"
if (!file.exists(mat_path)) {
  fwrite_mat(t(mat), mat_path)
}


# Output folder from iteratively performing GRNBoost2
grn_dir <- "/space/scratch/amorin/R_objects/arboreto_test/"
grn_files <- list.files(grn_dir, full.names = TRUE)



# GRNBoost2 gives a list of TF-target-importance scores. Convert to a gene x TF
# matrix of importance scores

network_list_to_mat <- function(network_list, tfs, genes, ncore) {
  
  imp_vec <- setNames(rep(0, length(genes)), genes) # init importance vector
  
  imp_l <- mclapply(tfs, function(tf) { # extract gene scores for each TF
    
    tf_df <- filter(network_list, V1 == tf)
    gene_match <- match(tf_df$V2, genes)
    imp_vec[gene_match] <- tf_df$V3
    imp_vec
    
  }, mc.cores = ncore)
  
  # Bind importance vectors into a matrix
  imp_mat <- do.call(cbind, imp_l)
  colnames(imp_mat) <- tfs
  
  return(imp_mat)
}



# Generate a network/importance score matrix for each output network list

ready_networks <- function(paths, tfs, genes, ncore) {
  
  mat_l <- lapply(paths, function(x) {
    
    network_list <- read.delim(x, header = FALSE)
    if (!all(network_list$V1 %in% tfs)) stop(paste(x, "doesn't have all TFs"))
    network_mat <- network_list_to_mat(network_list, tfs, genes, ncore)

  })
  
  return(mat_l)
}




mat_l <- ready_networks(grn_files, keep_tfs, keep_genes, ncore)

n_per_tf <- bind_rows(lapply(mat_l, function(x) colSums(x > 0)))

summ_n <- bind_rows(lapply(keep_tfs, function(x) {
  n <- n_per_tf[, x, drop = TRUE]
  data.frame(TF = x, min = min(n), median = median(n), max = max(n))
}))



ggplot(summ_n, aes(x = reorder(TF, median))) +
  geom_linerange(aes(ymin = min, ymax = max)) +
  geom_point(aes(x = TF, y = median, colour = "red"), shape = 21) +
  xlab("TF") +
  ylab("Count of associated genes") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 25),
        legend.position = "none")



summ_n %>% 
  mutate(Diff = max - min) %>% 
  ggplot(., aes(x = Diff)) +
  geom_histogram(bins = 100) +
  xlab("Max-min of associated genes") +
  ylab("Count of TFs") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))




# Topk overlap for each TF over all iterations

k <- 30

iter_topk <- mclapply(keep_tfs, function(tf) {
  
  tf_mat <- as.matrix(do.call(cbind, lapply(mat_l, function(x) x[, tf])))
  topk_mat <- colwise_topk_intersect(tf_mat, k = k)
  rownames(topk_mat) <- colnames(topk_mat) <- paste0("rep", 1:ncol(topk_mat))
  topk_df <- mat_to_df(topk_mat, value_name = "Topk")
  mean(topk_df$Topk)
}, mc.cores = ncore)



data.frame(Topk = unlist(iter_topk)) %>% 
  ggplot(., aes(x = Topk)) +
  geom_histogram(bins = 100) +
  xlab("Mean Top30 between iterations") +
  ylab("Count TFs") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 25),
        legend.position = "none")



# Average network

avg_network <- Reduce("+", mat_l) / length(mat_l)

# Prepare comparison with raw pearson's cor mat

diag(cormat[keep_tfs, keep_tfs]) <- NA


stopifnot(identical(rownames(avg_network), rownames(cormat)),
          identical(colnames(avg_network), colnames(cormat)))


compare_topk <- pair_colwise_topk(cormat, avg_network, k = k, ncores = ncore)
compare_topk_abs <- pair_colwise_topk(abs(cormat), avg_network, k = k, ncores = ncore)
compare_scor <- pair_colwise_cor(cormat, avg_network, ncores = ncore)
compare_scor_abs <- pair_colwise_cor(abs(cormat), avg_network, ncores = ncore)


compare_df <- data.frame(
  Symbol = keep_tfs, 
  Top30 = compare_topk,
  Top30_abs = compare_topk_abs,
  Scor = compare_scor,
  Scor_abs = compare_scor_abs
)


check_tf <- "RUNX1"

check_df <- data.frame(
  Symbol = keep_genes,
  Pcor = cormat[, check_tf], 
  GRN = avg_network[, check_tf])


cor(check_df$Pcor, check_df$GRN, method = "spearman", use = "pairwise.complete.obs")
cor(abs(check_df$Pcor), check_df$GRN, method = "spearman", use = "pairwise.complete.obs")


topk_intersect(
  slice_max(check_df, Pcor, n = k)$Symbol,
  slice_max(check_df, GRN, n = k)$Symbol
)


topk_intersect(
  slice_max(check_df, abs(Pcor), n = k)$Symbol,
  slice_max(check_df, GRN, n = k)$Symbol
)


plot(check_df$Pcor, check_df$GRN)


check_gene <- "GMPS"
plot(mat[check_tf, ], mat[check_gene, ])
cor(mat[check_tf, ], mat[check_gene, ])
cor(mat[check_tf, ], mat[check_gene, ], method = "spearman")


# check warning produced for "GTCTCGTTCGTATCAG-1-5382" 
