## This script performs a saturation analysis
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
library(ggrepel)
library(parallel)
library(egg)
library(pheatmap)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")

n_iter <- 100  # How many samples to run per step
k <- 200  # TopK overlap
set.seed(5)

# Dataset meta and data IDs
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)

# Relevant gene tables
pc_mm <- read.delim(ref_mm_path)
tfs_mm <- read.delim(tfs_mm_path)

# The final aggregate/global profile used for comparison
agg_mm <- readRDS("/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/aggregate_cormat_FZ_mm.RDS")

# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# Directory and file pattern for saved raw correlation matrices 
cmat_dir_mm <- "/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/"
pattern <- "_cormat.tsv"

# Output paths
ind_topk_mm_path <- file.path(data_out_dir, "microglia_individual_dataset_topk_mm.RDS")
saturation_mm_path <- file.path(data_out_dir, "microglia_coexpr_saturation_mm.RDS")


# Functions
# ------------------------------------------------------------------------------



# Loads aggregate correlation matrices or NA count matrices for the given dataset
# ids into a list. 
# Pattern: "_NA_mat.tsv" for NA counts

load_mat_to_list  <- function(ids,
                              dir,
                              pattern,  
                              genes,
                              sub_genes = NULL,
                              verbose = TRUE) {
  
  mat_l <- lapply(ids, function(id) {
    if (verbose) message(paste(id, Sys.time()))
    path <- file.path(dir, paste0(id, pattern))
    fread_to_mat(path, genes, sub_genes)
  })
  
  names(mat_l) <- ids
  gc(verbose = FALSE)
  
  return(mat_l)
}



# Call load_agg_mat_list to load/save a local copy

load_or_generate_agg <- function(path, 
                                 ids, 
                                 dir,
                                 pattern,
                                 genes, 
                                 sub_genes = NULL) {
  
  if (!file.exists(path)) {
    
    agg_l <- load_mat_to_list(ids = ids, 
                              dir = dir,
                              genes = genes, 
                              pattern = pattern,
                              sub_genes = sub_genes)
    
    saveRDS(agg_l, path)
  } else {
    agg_l <- readRDS(path)
  }
  
  return(agg_l)
}




# Subset cormat to min thresholded genes, NA->0 and FZ transform

ready_cmat <- function(cmat, keep_genes, keep_tfs) {
  
  cmat <- cmat[keep_genes, keep_tfs]
  cmat[is.na(cmat)] <- 0
  cmat[cmat > 1] <- 1
  cmat[cmat < -1] <- -1
  cmat <- fisherz(cmat)
  
}


# Summarize across columns of mat

colwise_summary <- function(mat) {
  
  summ_mat <- apply(mat, 2, function(x) {
    data.frame(
      Min = min(x, na.rm = TRUE),
      QR1 = quantile(x, 0.25, na.rm = TRUE),
      Median = median(x, na.rm = TRUE),
      Mean = mean(x, na.rm = TRUE),
      QR3 = quantile(x, 0.75, na.rm = TRUE),
      Max = max(x, na.rm = TRUE)
    )
  })
  
  summ_mat <- do.call(rbind, summ_mat)
  return(summ_mat)
}


# Get matrix for dividing each agg element by its count measured

count_msr <- function(na_mat_l) {
  
  na_mat <- Reduce("+", na_mat_l)
  msr_mat <- length(na_mat_l) - na_mat
  
  return(msr_mat)
}



# Running
# ------------------------------------------------------------------------------


# Ready aggregate used for comparison
keep_mm <- count_summ$Mouse$Filter_genes
keep_tfs_mm <- intersect(keep_mm, tfs_mm$Symbol)
agg_mat_mm <- agg_mm$Agg_mat[keep_mm, keep_tfs_mm]
msr_mm <- count_summ$Mouse$Summ_df


# Load and ready all correlation matrices
cmat_l <- load_mat_to_list(ids = ids_mm,
                           dir = cmat_dir_mm,
                           pattern = pattern,
                           genes = pc_mm$Symbol,
                           sub_genes = tfs_mm$Symbol)


cmat_l <- lapply(cmat_l, ready_cmat, keep_mm, keep_tfs_mm)


# Tracking NA/0 status for dividing aggregate by measured elements
na_mat_l <- lapply(cmat_l, function(x) {
  x <- x[keep_mm, keep_tfs_mm]
  x == 0
})



if (!file.exists(ind_topk_mm_path)) {
  
  # TopK overlap of each individual dataset with the global
  ind_topk_l <- lapply(cmat_l, function(mat) {
    pair_colwise_topk(mat1 = mat,
                      mat2 = agg_mat_mm,
                      k = k,
                      ncores = ncore)
  })
  
  ind_topk_mat <- do.call(cbind, ind_topk_l)
  
  saveRDS(ind_topk_mat, ind_topk_mm_path)
  
}



 

# TODO: function

steps <- 2:(length(ids_mm) - 1)



# By TF, keeping only datasets that measure the TF

step_l <- mclapply(keep_tfs_mm, function(tf) {
  
  message(paste(tf, Sys.time()))
  
  msr_ids <- names(which(!sapply(na_mat_l, function(x) x[tf, tf])))
  tf_mat <- as.matrix(bind_cols(lapply(cmat_l[msr_ids], function(x) x[, tf])))
  tf_mat[tf_mat == 0] <- NA  # TODO: rm and skip na->0 in ready mat if keeping
  
  rownames(tf_mat) <- keep_mm
  agg_vec <- topk_sort(agg_mat_mm[, tf], k = k)
  
  
  steps <- 2:(length(msr_ids) - 1)
  
  step_l <- lapply(steps, function(step) {
    
    iter_l <- lapply(1:n_iter, function(iter) {
      
      sample_ids <- sample(msr_ids, size = step, replace = FALSE)
      avg_sample <- rowMeans(tf_mat[, sample_ids], na.rm = TRUE)
      sample_vec <- topk_sort(avg_sample, k = k)
      topk_intersect(agg_vec, sample_vec)
      
    })
    
    iter_summ <- summary(unlist(iter_l))
    
  })
  
  step_summ <- cbind(do.call(rbind, step_l), N_step = steps)
  
}, mc.cores = ncore)
names(step_l) <- keep_tfs_mm


saveRDS(step_l, saturation_mm_path)





ind_topk_mat <- readRDS(ind_topk_mm_path)
step_l <- readRDS(saturation_mm_path)



# Summarize 
# ------------------------------------------------------------------------------


# Datasets that most resemble the global profiles


rank_topk <- t(colrank_mat(t(-ind_topk_mat)))
rank_topk <- rank_topk / max(rank_topk)

rank_topk <- rank_topk[order(rowMeans(rank_topk)), order(colMeans(rank_topk))]


avg_topk <- colMeans(ind_topk_mat)
avg_rank_topk <- colMeans(rank_topk)



is_10x <- c("10x 3' v1", "10x 3' v2", "10x 3' v2/v3", "10x 3' v3", "10x 5' v1", "10x 5' v2")


plot_df <- data.frame(
  ID = names(avg_rank_topk),
  Agree_global = avg_rank_topk
) %>% 
  left_join(mcg_meta, by = "ID") %>% 
  mutate(Is_10X = Platform %in% is_10x)



pheatmap(rank_topk,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE)


pheatmap(t(avg_rank_topk),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = "black",
         fontsize = 20,
         cellheight = 20)



ggplot(plot_df, aes(x = log10(N_cells), y = Agree_global)) +
  geom_point(shape = 21, size = 3) +
  geom_smooth(method = "lm") +
  xlab("Log10 count of cells") +
  ylab("Global agreement") +
  theme_classic() +
  theme(text = element_text(size = 30),
        plot.margin = margin(10, 20, 10, 10))


ggplot(plot_df, aes(x = Is_10X, y = Agree_global)) +
  geom_boxplot(width = 0.2) +
  geom_jitter(shape = 21, size = 2.4, fill = "slategrey", width = 0.05) +
  xlab("10X platform") +
  ylab("Global agreement") +
  theme_classic() +
  theme(text = element_text(size = 30),
        plot.margin = margin(10, 20, 10, 10))
  



# plot(log10(plot_df$N_cells), plot_df$Agree_global)
# plot(log10(plot_df$Median_UMI), plot_df$Agree_global)
# plot(log10(plot_df$N_msr_postfilt), plot_df$Agree_global)
# plot(log10(plot_df$N_msr_prefilt), plot_df$Agree_global)
# boxplot(plot_df$Agree_global ~ plot_df$Is_10X)
# 





# Proportion of experiments needed to achieve recovery


rec_df <- lapply(names(step_l), function(gene) {
  
  gene_df <- as.data.frame(step_l[[gene]])
  n_msr <- max(gene_df$N_step) + 1
  
  min_step <- gene_df %>% 
    filter(Mean >= (k * 0.8)) %>% 
    slice_min(N_step, n = 1) %>%
    pull(N_step)
  
  if (length(min_step) == 0) min_step <- 0
  
  data.frame(Symbol = gene,
             N_msr = n_msr,
             Min_n = min_step, 
             Min_prop = min_step / n_msr)
  
})


rec_df <- do.call(rbind, rec_df) %>% 
  as.data.frame() %>% 
  left_join(msr_mm, by = "Symbol")




hist(rec_df$Min_prop, breaks = 100)
hist(rec_df$Min_n, breaks = 100)
plot(rec_df$Min_prop, rec_df$QN_avg)
plot(rec_df$Min_n, rec_df$QN_avg)


ggplot(rec_df, aes(x = Min_prop)) +
  geom_histogram(bins = 100, fill = "slategrey") +
  ggtitle("Mouse") +
  xlab("Proportion recovery") +
  ylab("Count of TFs") +
  theme_classic() +
  theme(text = element_text(size = 30),
        plot.margin = margin(10, 20, 10, 10))



ggplot(rec_df, aes(x = QN_avg, y = Min_prop)) +
  geom_point(shape = 21, size = 3) +
  geom_smooth(method = "lm") +
  xlab("Mean log2 CPM") +
  ylab("Proportion recovery") +
  theme_classic() +
  theme(text = element_text(size = 20))



# Plotting
# ------------------------------------------------------------------------------


# Showing TopK overlap of individual experiments to global

plot_topk_ind <- function(gene, ind_topk_mat, k) {
  
  # Create a group label for the single experiment that had the best TopK
  plot_df <- data.frame(
    ID = colnames(ind_topk_mat),
    Topk = ind_topk_mat[gene, ]
  ) %>% 
    mutate(Label = ID %in% slice_max(., Topk, n = 1, with_ties = FALSE)$ID,
           Group = factor(1, levels = 1))
  
  
  ggplot(plot_df, aes(y = Topk, x = Group)) +
    geom_boxplot(width = 0.1) +
    geom_jitter(aes(x = Group, y = Topk), 
                shape = 21, colour = "slategrey", width = 0.02, height = 0) +
    geom_label_repel(data = filter(plot_df, Label),
                     aes(x = Group, y = Topk, label = ID),
                     nudge_x = -0.3, nudge_y = 5, size = 6) +
    ylim(c(0, k)) +
    ylab(expr("Top"[!!k])) +
    xlab("Individual experiments") +
    ggtitle(gene) +
    theme_classic() +
    theme(text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}


# Showing spread of TopK across iterations per step

plot_topk_steps <- function(gene, step_l, k) {
  
  # Isolating steps as well as minimum steps to achieve recovery from summary df
  plot_df <- as.data.frame(step_l[[gene]])
  
  min_step <- plot_df %>% 
    filter(Mean >= (k * 0.8)) %>% 
    slice_min(N_step, n = 1) %>%
    pull(N_step)
  
  ggplot(plot_df, aes(x = N_step, y = Mean)) +
    geom_crossbar(aes(x = N_step, ymin = `1st Qu.`, ymax = `3rd Qu.`)) +
    geom_point(shape = 19, colour = "firebrick") +
    geom_vline(xintercept = min_step, linetype = "dashed", colour = "black") +
    geom_hline(yintercept = (k * 0.8), linetype = "dashed", colour = "black") +
    ylim(c(0, k)) +
    ylab(expr("Top"[!!k])) +
    xlab("Count of sampled experiments") +
    ggtitle(gene) +
    theme_classic() +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
          plot.margin = margin(10, 20, 10, 10))
}




gene <- "Runx1"
p1a <- plot_topk_ind(gene, ind_topk_mat, k)
p1b <- plot_topk_steps(gene, step_l, k) + ylab(NULL)
p1 <- ggarrange(p1a, p1b, nrow = 1, widths = c(0.33, 1))

ggsave(p1, height = 6, width = 16, device = "png", dpi = 300,
    filename = file.path(plot_dir, paste0(gene, "_dataset_saturation.png")))
