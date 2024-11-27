## GRNBoost2 was ran 100 times on each microglia dataset. For each dataset ID,
## load each iterations and average into one matrix. Then, generate one global
## matrix that averages all dataset matrices.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/grnboost2_utils.R")

grn_dir <- file.path(data_out_dir, "GRNBoost2")

mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_hg <- filter(mcg_meta, Species == "Human")$ID
ids_mm <- filter(mcg_meta, Species == "Mouse")$ID

# Load genes/TFs
tf_hg <- read.delim("/home/amorin/Data/Metadata/TFs_human.tsv")
pc_hg <- read.delim(ref_hg_path)


# TODO: remove/reduce when everything is ran
# GRN output for each dataset
grn_files_list <- lapply(ids_hg, function(x) {
  list.files(file.path(grn_dir, x), full.names = TRUE, pattern = ".*_iter.*")
})

names(grn_files_list) <- ids_hg
grn_files_list <- grn_files_list[sapply(grn_files_list, length) > 0]
ids_hg <- names(grn_files_list)


# TODO: coerce common dimensions earlier in list, for ease of downstream comp.

avg_grn_l <- average_and_save_each_network(ids = ids_hg,
                                           grn_dir = grn_dir,
                                           ncore = ncore)


avg_grn_mat <- average_all_networks(avg_grn_l)



# Load cormat

keep_tfs <- colnames(avg_grn_mat)
keep_genes <- rownames(avg_grn_mat)


cmat_l <- lapply(ids_hg, function(x) {
  fread_to_mat(file.path(cmat_dir_hg, paste0(x, "_cormat.tsv")), 
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
names(fz_l) <- ids_hg

fz_mat <- Reduce("+", fz_l) / msr_mat
fz_mat[is.infinite(fz_mat)] <- NA



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




check_tf <- "RUNX1"


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
check_gene <- "RCN2"


check_data <- lapply(ids_hg, function(x) {
  
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




# Requires loading of demo_Rcistarget.R (I know, I know)

coexpr_genes <- slice_max(check_df, GRN, n = n)$Symbol

res_coexpr <- cisTarget(coexpr_genes, 
                        motifRankings = motif_rank, 
                        motifAnnot = motif_anno)

check_tf_motif_score <- lapply(check_df$Symbol, function(check_gene) {
  
  if (check_gene %!in% colnames(motif_rank_hg@rankings)) {
    return(NA)
  }
  
  rank <- left_join(
    motif_rank_hg@rankings[, c("motifs", check_gene)],
    motif_score_hg@rankings[, c("motifs", check_gene)],
    by = "motifs",
    suffix = c("_rank", "_score")
  ) %>%
    left_join(motif_anno_hg, by = c("motifs" = "motif")) %>%
    filter(TF == check_tf) %>% 
    slice_min(!!sym(paste0(check_gene, "_rank")), n = 1) %>% 
    pull(!!sym(paste0(check_gene, "_rank"))) %>% 
    unique()
  
  data.frame(Symbol = check_gene, Rank = rank)
  
})
names(check_tf_motif_score) <- coexpr_genes

check_tf_motif_score <- check_tf_motif_score[!is.na(check_tf_motif_score)]


check_tf_motif_score <- bind_rows(check_tf_motif_score)

check_df <- left_join(check_df, check_tf_motif_score, by = "Symbol")
check_df$Top_motif <- check_df$Rank <= 200
check_df$Top_motif[is.na(check_df$Top_motif)] <- FALSE



check_df <- left_join(check_df, tf_rank_hg[[check_tf]], by = "Symbol")



ggplot() +
  geom_point(data = filter(check_df, !Top_motif), aes(x = FZ, y = GRN), shape = 21, size = 2.2, fill = "white", alpha = 0.4) +
  geom_point(data = filter(check_df, Top_motif), aes(x = FZ, y = GRN), shape = 21, size = 3.2, fill = "goldenrod") +
  ggtitle(check_tf) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "firebrick") +
  scale_fill_manual(values = c("white", "goldenrod")) +
  theme_classic() +
  theme(text = element_text(size = 20))
