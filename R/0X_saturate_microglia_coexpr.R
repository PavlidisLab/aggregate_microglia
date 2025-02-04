## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(aggtools)
library(ggrepel)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")

n_iter <- 100  # How many samples to run per step
k <- 200  # TopK overlap

# Dataset meta and data IDs
mcg_meta <- read.delim(mcg_meta_dedup_path)
ids_hg <- unique(filter(mcg_meta, Species == "Human")$ID)
ids_mm <- unique(filter(mcg_meta, Species == "Mouse")$ID)

# Relevant gene tables
pc_mm <- read.delim(ref_mm_path)
tfs_mm <- read.delim("/home/amorin/Data/Metadata/AnimalTFDB_mouse_V4.tsv")

# The final aggregate/global profile used for comparison
agg_mm <- readRDS("/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/aggregate_cormat_FZ_mm.RDS")

# List of measurement info to keep filtered genes
count_summ <- readRDS(mcg_count_summ_list_path)

# Directory and file pattern for saved raw correlation matrices 
cmat_dir_mm <- "/space/scratch/amorin/aggregate_microglia/Cormats/Mm_pcor/"
pattern <- "_cormat.tsv"

# Output paths
ind_topk_path <- file.path(data_out_dir, "microglia_individual_dataset_topk_mm.RDS")
saturation_path <- file.path(data_out_dir, "microglia_coexpr_saturation_mm.RDS")


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
  
  cmat <- cmat[keep_mm, keep_tfs_mm]
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



# Ready aggregate used for comparison
keep_mm <- count_summ$Mouse$Filter_genes
keep_tfs_mm <- intersect(keep_mm, tfs_mm$Symbol)
agg_mat <- agg_mm$Agg_mat[keep_mm, keep_tfs_mm]


# Load and ready all correlation matrices
cmat_l <- load_mat_to_list(ids = ids_mm,
                           dir = cmat_dir_mm,
                           pattern = pattern,
                           genes = pc_mm$Symbol,
                           sub_genes = tfs_mm$Symbol)


cmat_l <- lapply(cmat_l, ready_cmat, keep_mm, keep_tfs_mm)



# TopK overlap of each individual dataset with the global
ind_topk_l <- lapply(cmat_l, function(mat) {
  pair_colwise_topk(mat1 = mat,
                    mat2 = agg_mat,
                    k = k,
                    ncores = ncore)
})



ind_topk_mat <- do.call(cbind, ind_topk_l)
 


# TODO: need to account for NA when dividing if want to be consistent with FZ agg

steps <- 2:(length(ids_mm)-1)


step_l <- lapply(steps, function(step) {
  
  message(paste("Step ==", step, Sys.time()))
  
  iter_l <- lapply(1:n_iter, function(iter) {
    
    sample_ids <- sample(ids_mm, size = step, replace = FALSE)
    avg_sample <- Reduce("+", cmat_l[sample_ids]) / length(sample_ids)
    pair_colwise_topk(mat1 = avg_sample, mat2 = agg_mat, k = k, ncores = ncore)
    
  })
  
  iter_mat <- do.call(rbind, iter_l)
  iter_summ <- colwise_summary(iter_mat)

  return(iter_summ)
  
})



saveRDS(ind_topk_mat, ind_topk_path)
saveRDS(step_l, saturation_path)


ind_topk_mat <- readRDS(ind_topk_path)
step_l <- readRDS(saturation_path)





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
    theme(text = element_text(size = 30),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}


# Showing spread of TopK across iterations per step

plot_topk_steps <- function(gene, step_l, k) {
  
  # Isolating steps as well as minimum steps to achieve recovery from summary df
  plot_df <- bind_rows(lapply(step_l, function(x) x[gene, ])) %>% 
    mutate(N_step = steps)
  
  min_step <- plot_df %>% 
    filter(Mean >= (k * 0.8)) %>% 
    slice_min(N_step, n = 1) %>%
    pull(N_step)
  
  ggplot(plot_df, aes(x = N_step, y = Mean)) +
    geom_crossbar(aes(x = N_step, ymin = QR1, ymax = QR3)) +
    geom_point(shape = 19, colour = "firebrick") +
    geom_vline(xintercept = min_step, linetype = "dashed", colour = "black") +
    geom_hline(yintercept = (k * 0.8), linetype = "dashed", colour = "black") +
    ylim(c(0, k)) +
    ylab(expr("Top"[!!k])) +
    xlab("Count of sampled experiments") +
    ggtitle(gene) +
    theme_classic() +
    theme(text = element_text(size = 30),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
          plot.margin = margin(10, 20, 10, 10))
}


gene <- "Runx1"
plot_topk_ind(gene, ind_topk_mat, k)
plot_topk_steps(gene, step_l, k)









