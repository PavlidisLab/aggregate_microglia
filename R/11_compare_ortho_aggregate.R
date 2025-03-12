## Basic exploration of aggregate profiles between species
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(ggrepel)
source("R/00_config.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")

k <- 200
set.seed(5)

mcg_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)
ids_mm <- filter(mcg_meta, Species == "Mouse")$ID
ids_hg <- filter(mcg_meta, Species == "Human")$ID

# Gene tables
pc_hg <- read.delim(ref_hg_path)
pc_mm <- read.delim(ref_mm_path)
tfs_hg <- read.delim(tfs_hg_path)
tfs_mm <- read.delim(tfs_mm_path)
pc_ortho <- read.delim(pc_ortho_path)

# Using FZ aggregates
fz_hg <- readRDS(mcg_fz_hg_path)
fz_mm <- readRDS(mcg_fz_mm_path)


# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)


keep_ortho <- filter(pc_ortho, 
                     Symbol_hg %in% count_summ$Human$Filter_genes &
                     Symbol_mm %in% count_summ$Mouse$Filter_genes)


keep_tfs <- filter(keep_ortho, 
                   Symbol_hg %in% tfs_hg$Symbol & 
                   Symbol_mm %in% tfs_mm$Symbol)


# Subsetting agg matrices to ortho genes
mat_mm <- fz_mm$Agg_mat[keep_ortho$Symbol_mm, keep_ortho$Symbol_mm]
rownames(mat_mm) <- colnames(mat_mm) <- keep_ortho$Symbol_hg
mat_hg <- fz_hg$Agg_mat[keep_ortho$Symbol_hg, keep_ortho$Symbol_hg]

# Topk overlap between species
topk <- pair_colwise_topk(mat_mm, mat_hg, k = k, ncores = ncore)

# Null overlap all ortho genes
null_topk <- pair_shuffle_topk(mat_mm, mat_hg, k = k, ncores = ncore)

# Null overlap ortho TF genes
null_tfs_topk <- pair_shuffle_topk(mat_mm[, keep_tfs$Symbol_hg], 
                                   mat_hg[, keep_tfs$Symbol_hg], 
                                   k = k, 
                                   ncores = ncore)

# Spearman cor between ortho profiles  
scor <- pair_colwise_cor(mat_mm, mat_hg, ncores = ncore)

# Null spearman all ortho genes
null_scor <- pair_shuffle_cor(mat_mm, mat_hg, ncores = ncore)

# Null spearman ortho TF genes
null_tfs_scor <- pair_shuffle_cor(mat_mm[, keep_tfs$Symbol_hg], mat_hg[, keep_tfs$Symbol_hg], ncores = ncore)

# Organize similarities
sim_df <- data.frame(
  Symbol = keep_ortho$Symbol_hg,
  Is_TF = keep_ortho$Symbol_hg %in% keep_tfs$Symbol_hg,
  Topk = topk,
  Scor = scor
)

# Most consistent ortho TF profiles
top_tfs <- sim_df %>% filter(Is_TF) %>% slice_max(Topk, n = 20) %>% pull(Symbol)

sim_df$Is_top_TF <- sim_df$Symbol %in% top_tfs


# Hist of topk overlaid with null
p1 <- ggplot(sim_df, aes(x = Topk)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = median(null_topk)) +
  theme_classic() +
  theme(text = element_text(size = 25),
        plot.margin = margin(c(10, 20, 10, 10)))


# Scatter of topk vs scor, emphasizing TFs and labeling top TFs
p2 <- ggplot(sim_df, aes(x = Topk, y = Scor)) +
  geom_point(data = filter(sim_df, !Is_TF), 
             shape = 21, size = 2.4, alpha = 0.2) +
  geom_point(data = filter(sim_df, Is_TF), 
             shape = 21, size = 3.4, colour = "firebrick") +
  geom_text_repel(data = filter(sim_df, Is_top_TF),
                  aes(label = Symbol),
                  size = 6, max.overlaps = length(top_tfs)) +
  xlab("Top200") +
  ylab("Spearman's cor") +
  theme_classic() +
  theme(text = element_text(size = 25),
        plot.margin = margin(c(10, 20, 10, 10)))
