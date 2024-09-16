## TODO
## TODO: ccre -> GR -> intersect with SPI1 peaks
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(GenomicRanges)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Shuffle procedure for null comparison
set.seed(78)

# Pathing
mcg_root_hg <- file.path(amat_dir, "Microglia", "Hg_pcor")
mcg_root_mm <- file.path(amat_dir, "Microglia", "Mm_pcor")

# mcg_path_hg <- file.path(mcg_root_hg, "aggregate_matrix.tsv")
mcg_path_hg <- file.path(mcg_root_hg, "aggregate_matrix_fishersz.tsv")
# mcg_path_mm <- file.path(mcg_root_mm, "aggregate_matrix.tsv")
mcg_path_mm <- file.path(mcg_root_mm, "aggregate_matrix_fishersz.tsv")
mcg_na_path_hg <- file.path(mcg_root_hg, "NA_matrix.tsv")
mcg_na_path_mm <- file.path(mcg_root_mm, "NA_matrix.tsv")

# Individual dataset cor matrice paths 
mcg_ids_hg <- list.files(mcg_root_hg, full.names = TRUE, pattern = "cormat")
mcg_ids_mm <- list.files(mcg_root_mm, full.names = TRUE, pattern = "cormat")

# TODO
ccre_path_hg <- "/space/grp/amorin/human_microglia_ccre_gene_cors.txt"
ccre_path_mm <- "/space/grp/amorin/mouse_microglia_ccre_gene_cors.txt"
ccre_hg <- read.delim(ccre_path_hg, stringsAsFactors = FALSE)
ccre_mm <- read.delim(ccre_path_mm, stringsAsFactors = FALSE)

chip_path_hg <- "/space/scratch/amorin/R_objects/unibind_grlist_perm_human.RDS"
chip_path_mm <- "/space/scratch/amorin/R_objects/unibind_grlist_perm_mouse.RDS"
chip_hg <- readRDS(chip_path_hg)
chip_mm <- readRDS(chip_path_mm)

# Average bind scores and bind matrix
bind_summary <- readRDS(bind_summary_path)
bind_dat <- readRDS(bind_dat_path)

# Protein coding genes
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
tfs_hg <- read.delim(tfs_hg_path, stringsAsFactors = FALSE)
tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)

# Aggregate matrices
mcg_hg <- fread_to_mat(mcg_path_hg, genes = pc_hg$Symbol)
mcg_mm <- fread_to_mat(mcg_path_mm, genes = pc_mm$Symbol)
coexpr_hg <- readRDS(rank_tf_hg_path)
coexpr_mm <- readRDS(rank_tf_mm_path)

# Matrices tracking count of NAs
mcg_na_hg <- fread_to_mat(mcg_na_path_hg, genes = pc_hg$Symbol)
mcg_na_mm <- fread_to_mat(mcg_na_path_mm, genes = pc_mm$Symbol)


# FZ: Inf diag to NA and all NA to 0
# TODO: upstream?

mcg_hg <- diag_to_na(mcg_hg)
mcg_hg[is.na(mcg_hg)] <- 0

mcg_mm <- diag_to_na(mcg_mm)
mcg_mm[is.na(mcg_mm)] <- 0


# Remove low/unmeasured genes
keep_hg <- names(which(diag(mcg_na_hg) < floor(max(mcg_na_hg) * 0.9)))
keep_mm <- names(which(diag(mcg_na_mm) < floor(max(mcg_na_mm) * 0.9)))
mcg_sub_hg <- mcg_hg[keep_hg, keep_hg]
mcg_sub_mm <- mcg_mm[keep_mm, keep_mm]


gene_hg <- "SPI1"
gene_mm <- "Spi1"



# Sim of matched aggregate relative to all bind
# ------------------------------------------------------------------------------



# TODO: into function, if keeping
# Looking at how well matched aggregate vectors agree relative to rest



# All binding versus one TF microglia coexpr

sim_mcg_hg <- mclapply(colnames(bind_summary$Human_TF), function(x) {
  
  df <- data.frame(
    Coexpr = mcg_hg[keep_hg, gene_hg],
    Bind = bind_summary$Human_TF[keep_hg, x]
  )
  
  topk <- topk_intersect(
    topk_sort(mcg_hg[keep_hg, gene_hg], k = 200),
    topk_sort(bind_summary$Human_TF[keep_hg, x], k = 200, check_k_arg = FALSE)
  )
  
  scor <-  cor(df$Coexpr, df$Bind, method = "spearman", use = "pairwise.complete.obs")
  
  return(data.frame(Symbol = x, Top200 = topk, Scor = scor))
  
}, mc.cores = ncore)

sim_mcg_hg <- do.call(rbind, sim_mcg_hg)




sim_mcg_mm <- mclapply(colnames(bind_summary$Mouse_TF), function(x) {
  
  df <- data.frame(
    Coexpr = mcg_mm[keep_mm, gene_mm],
    Bind = bind_summary$Mouse_TF[keep_mm, x]
  )
  
  topk <- topk_intersect(
    topk_sort(mcg_mm[keep_mm, gene_mm], k = 200),
    topk_sort(bind_summary$Mouse_TF[keep_mm, x], k = 200, check_k_arg = FALSE)
  )
  
  scor <-  cor(df$Coexpr, df$Bind, method = "spearman", use = "pairwise.complete.obs")
  
  return(data.frame(Symbol = x, Top200 = topk, Scor = scor))
  
}, mc.cores = ncore)

sim_mcg_mm <- do.call(rbind, sim_mcg_mm)




ggplot(sim_mcg_hg, aes(x = Scor)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = filter(sim_mcg_hg, Symbol == gene_hg)[["Scor"]], col = "red") +
  ggtitle("Human SPI1") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(sim_mcg_hg, aes(x = Top200)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = filter(sim_mcg_hg, Symbol == gene_hg)[["Top200"]], col = "red") +
  ggtitle("Human SPI1") +
  theme_classic() +
  theme(text = element_text(size = 25))



ggplot(sim_mcg_mm, aes(x = Scor)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = filter(sim_mcg_mm, Symbol == gene_mm)[["Scor"]], col = "red") +
  ggtitle("Mouse SPI1") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(sim_mcg_mm, aes(x = Top200)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = filter(sim_mcg_mm, Symbol == gene_mm)[["Top200"]], col = "red") +
  ggtitle("Mouse SPI1") +
  theme_classic() +
  theme(text = element_text(size = 25))






# All binding versus one TF all coexpr


sim_all_hg <- mclapply(colnames(bind_summary$Human_TF), function(x) {
  
  df <- coexpr_hg[[gene_hg]]
  df$Bind <- bind_summary$Human_TF[df$Symbol, x]
  
  topk <- topk_intersect(
    topk_sort(setNames(df$Avg_aggr_coexpr, df$Symbol), k = 200),
    topk_sort(setNames(df$Bind, df$Symbol), k = 200, check_k_arg = FALSE)
  )
  
  scor <-  cor(df$Avg_aggr_coexpr, df$Bind, method = "spearman", use = "pairwise.complete.obs")
  
  return(data.frame(Symbol = x, Top200 = topk, Scor = scor))
  
}, mc.cores = ncore)

sim_all_hg <- do.call(rbind, sim_all_hg)



sim_all_mm <- mclapply(colnames(bind_summary$Mouse_TF), function(x) {
  
  df <- coexpr_mm[[gene_mm]]
  df$Bind <- bind_summary$Mouse_TF[df$Symbol, x]
  
  topk <- topk_intersect(
    topk_sort(setNames(df$Avg_aggr_coexpr, df$Symbol), k = 200),
    topk_sort(setNames(df$Bind, df$Symbol), k = 200, check_k_arg = FALSE)
  )
  
  scor <-  cor(df$Avg_aggr_coexpr, df$Bind, method = "spearman", use = "pairwise.complete.obs")
  
  return(data.frame(Symbol = x, Top200 = topk, Scor = scor))
  
}, mc.cores = ncore)

sim_all_mm <- do.call(rbind, sim_all_mm)



# Microglia SPI1 ChIP-seq
# ------------------------------------------------------------------------------


spi_cmat_hg <- readRDS(file.path(mcg_root_hg, "SPI1_rawcor.tsv"))
spi_cmat_mm <- readRDS(file.path(mcg_root_mm, "Spi1_rawcor.tsv"))

comsr_hg <- rowSums(!is.na(spi_cmat_hg))
comsr_mm <- rowSums(!is.na(spi_cmat_mm))

min_keep_hg <- ceiling(max(comsr_hg) * 0.1)
min_keep_mm <- ceiling(max(comsr_mm) * 0.1)


# Human SPI1 no microglia data
all_bind_ids_hg <- bind_dat$Permissive_hg$Meta %>% 
  filter(Symbol == gene_hg) %>% 
  pull(ID)


mcg_bind_ids_mm <- bind_dat$Permissive_mm$Meta %>% 
  filter(str_to_title(Symbol) == gene_mm & str_detect(ID, ".*microglia.*")) %>% 
  pull(ID)


all_bind_ids_mm <- bind_dat$Permissive_mm$Meta %>% 
  filter(str_to_title(Symbol) == gene_mm) %>% 
  pull(ID)



# TODO: show agreement of SPI1 coexpr aggregation between:
# 1) SPI1 microglia bind profiles
# 2) SPI1 all bind profiles (all + shuffle same n)


all_bind_hg <- bind_dat$Permissive_hg$Mat_qnl[, all_bind_ids_hg]
all_bind_mm <- bind_dat$Permissive_mm$Mat_qnl[, all_bind_ids_mm]
mcg_bind_mm <- bind_dat$Permissive_mm$Mat_qnl[, mcg_bind_ids_mm]


bind_df_hg <- data.frame(
  Symbol = pc_hg$Symbol,
  Coexpr_all = coexpr_hg[[gene_hg]][pc_hg$Symbol, "Avg_aggr_coexpr"],
  Coexpr_mcg = mcg_hg[pc_hg$Symbol, gene_hg],
  Bind_all = bind_summary$Human_TF[pc_hg$Symbol, gene_hg],
  N_na = mcg_na_hg[pc_hg$Symbol, gene_hg],
  Comsr = comsr_hg
)
bind_df_hg <- filter(bind_df_hg, Symbol != gene_hg)
bind_df_hg <- filter(bind_df_hg, Comsr >= min_keep_hg)




bind_df_mm <- data.frame(
  Symbol = pc_mm$Symbol,
  Coexpr_all = coexpr_mm[[gene_mm]][pc_mm$Symbol, "Avg_aggr_coexpr"],
  Coexpr_mcg = mcg_mm[pc_mm$Symbol, gene_mm],
  Bind_all = bind_summary$Mouse_TF[pc_mm$Symbol, gene_mm],
  Bind_wo_mcg = rowMeans(all_bind_mm[pc_mm$Symbol, setdiff(all_bind_ids_mm, mcg_bind_ids_mm)]),
  Bind_mcg = rowMeans(mcg_bind_mm[pc_mm$Symbol, ]),
  N_na = mcg_na_mm[pc_mm$Symbol, gene_mm],
  Comsr = comsr_mm
)
bind_df_mm <- filter(bind_df_mm, Symbol != gene_mm)
bind_df_mm <- filter(bind_df_mm, Comsr >= min_keep_mm)


plot(bind_df_hg$Coexpr_mcg, bind_df_hg$Bind_all)
plot(bind_df_hg$Coexpr_all, bind_df_hg$Bind_all)
plot(bind_df_mm$Coexpr_mcg, bind_df_mm$Bind_all)
plot(bind_df_mm$Coexpr_all, bind_df_mm$Bind_all)


cor(select_if(bind_df_hg, is.numeric), method = "spearman", use = "pairwise.complete.obs")
cor(select_if(bind_df_mm, is.numeric), method = "spearman", use = "pairwise.complete.obs")

colwise_topk_intersect(select_if(bind_df_hg, is.numeric), k = 200, check_k_arg = FALSE)
colwise_topk_intersect(select_if(bind_df_mm, is.numeric), k = 200, check_k_arg = FALSE)




bind_df_hg$Rank_coexpr_all <- rank(-bind_df_hg$Coexpr_all)
bind_df_hg$Rank_coexpr_mcg <- rank(-bind_df_hg$Coexpr_mcg)
bind_df_hg$Rank_bind_all <- rank(-bind_df_hg$Bind_all)


bind_df_mm$Rank_coexpr_all <- rank(-bind_df_mm$Coexpr_all)
bind_df_mm$Rank_coexpr_mcg <- rank(-bind_df_mm$Coexpr_mcg)
bind_df_mm$Rank_bind_all <- rank(-bind_df_mm$Bind_all)
bind_df_mm$Rank_bind_wo_mcg <- rank(-bind_df_mm$Bind_wo_mcg)
bind_df_mm$Rank_bind_mcg <- rank(-bind_df_mm$Bind_mcg)
bind_df_mm$Rank_bind_diff <- bind_df_mm$Rank_bind_wo_mcg - bind_df_mm$Rank_bind_mcg




bind_df_mm %>%
  filter(Rank_coexpr_mcg <= 200 &
        (Rank_bind_mcg <= 200 | Rank_bind_wo_mcg <= 200 | Rank_bind_all <= 200)) %>%
  view()


# Repression?
bind_df_mm %>%
  filter(Rank_bind_all <= 2000) %>%
  slice_max(Rank_coexpr, n = 200) %>%
  view()



# Bound regions as GRanges
all_gr_hg <- GRangesList(chip_hg[all_bind_ids_hg])
all_gr_mm <- GRangesList(chip_mm[all_bind_ids_mm])
mcg_gr_mm <- GRangesList(chip_mm[mcg_bind_ids_mm])

n_peaks_all_hg <- unlist(lapply(all_gr_hg, length))
n_peaks_all_mm <- unlist(lapply(all_gr_mm, length))
n_peaks_mcg_mm <- unlist(lapply(mcg_gr_mm, length))


# cCRE as Granges. Proximal == gene promoter, distal == assumed enhancer.
# Want to intersect distal regions with ChIP-seq to nominate binding

ccre_gr_hg <- ccre_hg %>% 
  dplyr::rename(
    Chromosome = Chromosome_distal,
    Start = Start_distal,
    End = End_distal) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


ccre_gr_mm <- ccre_mm %>% 
  dplyr::rename(
    Chromosome = Chromosome_distal,
    Start = Start_distal,
    End = End_distal) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


# TODO: filter
setdiff(ccre_gr_hg$Symbol, pc_hg$Symbol)
setdiff(ccre_gr_mm$Symbol, pc_mm$Symbol)


n_ccre_genes_hg <- n_distinct(ccre_gr_hg$Symbol)
n_ccre_genes_mm <- n_distinct(ccre_gr_mm$Symbol)

n_ccre_per_gene_hg <- sort(table(ccre_gr_hg$Symbol), decreasing = TRUE)
n_ccre_per_gene_mm <- sort(table(ccre_gr_mm$Symbol), decreasing = TRUE)




# Overlap between ChIP and cCREs
ol_all_hg <- findOverlaps(all_gr_hg, ccre_gr_hg)
ol_all_mm <- findOverlaps(all_gr_mm, ccre_gr_mm)
ol_mcg_mm <- findOverlaps(mcg_gr_mm, ccre_gr_mm)


n_hits_per_all_hg <- table(ol_all_hg@from)
n_hits_per_all_mm <- table(ol_all_mm@from)
n_hits_per_mcg_mm <- table(ol_mcg_mm@from)

# As expected, more peaks more cCRE hits
plot(as.integer(n_hits_per_all_hg), n_peaks_all_hg)
plot(as.integer(n_hits_per_mcg_mm), n_peaks_mcg_mm)


# 2 mouse experiments no hits ==> small count of peaks
all_bind_ids_mm[c(165, 170)]
c(165, 170) %in% ol_all_mm@from
n_peaks_all_mm[c(165, 170)]


# Genes that have the most hits
tt1 <- sort(table(ol_mcg_mm@to), decreasing = TRUE)
# tt1 <- sort(table(ol_all_mm@to), decreasing = TRUE)

gene_hits_mcg_mm <- data.frame(
  Symbol = ccre_gr_mm$Symbol[as.integer(names(tt1))],
  ID = ccre_gr_mm$cCRE_ID[as.integer(names(tt1))],
  N_hits = as.integer(tt1)
)


head(sort(table(gene_hits_mcg_mm$Symbol), decreasing = TRUE), 30)
head(sort(table(gene_hits_mcg_mm$ID), decreasing = TRUE), 30)
n_distinct(gene_hits_mcg_mm$Symbol)
n_distinct(gene_hits_mcg_mm$ID)


bind_df_mm$cCRE_hit <- bind_df_mm$Symbol %in% gene_hits_mcg_mm$Symbol


ggplot(bind_df_mm, aes(x = cCRE_hit, y = Coexpr_mcg)) + 
  geom_boxplot(width = 0.3, fill = "slategrey", alpha = 0.3) +
  ggtitle("Mouse Spi1") +
  ylab("Average Fisher's Z") +
  xlab("Has >= 1 cCRE-SPI1 hit") +
  theme_classic() +
  theme(text = element_text(size = 25))

wilcox.test(bind_df_mm$Coexpr_mcg ~ bind_df_mm$cCRE_hit)




ol_mcg_mm[ol_mcg_mm@to == 11]
ccre_gr_mm[11]                  

# TODO: confirm order of overlap and get count method. count in subject or y/n hit?
# TODO: also I'm sure you've already written similar code, find it




# TopK between all genes
# ------------------------------------------------------------------------------


topk_path <- "~/scratch/R_objects/microglia_aggr_topk_FZ.RDS"

if (!file.exists(topk_path)) {
  topk_hg <- colwise_topk_intersect(mcg_sub_hg, k = 200)
  topk_mm <- colwise_topk_intersect(mcg_sub_mm, k = 200)
  saveRDS(list(Human = topk_hg, Mouse = topk_mm), file = topk_path)
} else {
  topk_l <- readRDS(topk_path)
}



# Assume that if two genes have a high Top K (share a lot of coexpressed genes),
# they themselves are also coexpressed.
# TODO: double check co-msr


topk_df_hg <- data.frame(
  Symbol = keep_hg,
  Agg = mcg_hg[keep_hg, gene_hg],
  Topk = topk_l$Human[, gene_hg],
  N_na = mcg_na_hg[keep_hg, gene_hg])

topk_df_hg <- filter(topk_df_hg, Symbol != gene_hg)

plot(topk_df_hg$Agg, topk_df_hg$Topk)
cor(topk_df_hg$Agg, topk_df_hg$Topk, method = "spearman", use = "pairwise.complete.obs")



topk_df_mm <- data.frame(
  Symbol = keep_mm,
  Agg = mcg_mm[keep_mm, gene_mm],
  Topk = topk_l$Mouse[, gene_mm],
  N_na = mcg_na_mm[keep_mm, gene_mm])

topk_df_mm <- filter(topk_df_mm, Symbol != gene_mm)

plot(topk_df_mm$Agg, topk_df_mm$Topk)
cor(topk_df_mm$Agg, topk_df_mm$Topk, method = "spearman", use = "pairwise.complete.obs")
