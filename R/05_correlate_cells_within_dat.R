## Calculate the cell to cell correlation within each microglia dataset, to 
## compare how "consistent" the microglial cells are across data
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
library(egg)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")

set.seed(5)

# Dataset meta
mcg_meta <- read.delim(mcg_meta_dedup_path)


# Note: I inspected taking cor after filtering away 0 genes, but saw trivial 
# change to summary stats so I'm using the full set of genes

# Note: I was getting errors with some of the larger datasets, so I sample
# 1000 cells X10 and summarize this info for datasets with >= 10,000 cells

summarize_cell_cor <- function(dat_l, 
                               n_cell_threshold = 10e3, 
                               n_sample_cells = 1e3,
                               n_samples = 10,
                               verbose = TRUE,
                               ncores = 1) {
  
  ids <- names(dat_l)
  
  summ_cor_l <- mclapply(ids, function(x) {
    
    if (verbose) message(paste(x, Sys.time()))
    
    mat <- dat_l[[x]]$Mat
    ncells <- ncol(mat)
    
    if (ncells < n_cell_threshold) {
      
      cell_cor <- aggtools::sparse_pcor(dat_l[[x]]$Mat)
      cor_df <- mat_to_df(cell_cor, value_name = "Cor")
      
    } else {
      
      sample_l <- lapply(1:n_samples, function(i) {
        sample_cells <-  sample(colnames(mat), n_sample_cells, replace = FALSE)
        submat <- mat[, sample_cells]
        cell_cor <- aggtools::sparse_pcor(submat)
        mat_to_df(cell_cor, value_name = "Cor")
      })
      
      cor_df <- dplyr::bind_rows(sample_l)
      
    }
    
    summary(cor_df$Cor)
    
  }, mc.cores = ncores)
  names(summ_cor_l) <- ids
  
  return(summ_cor_l)
}



# Run/save/load
if (!file.exists(cell_cor_path)) {
  dat_l <- readRDS(mcg_dat_path) # List of count matrices and their metadata
  cor_l <- summarize_cell_cor(dat_l, ncores = ncore)
  saveRDS(cor_l, cell_cor_path)
} else {
  cor_l <- readRDS(cell_cor_path)
}



# Setting up summary df for plotting

plot_df <- do.call(rbind, cor_l) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(., mcg_meta, by = "ID") %>% 
  arrange(Median) %>%
  mutate(ID = factor(ID, levels = unique(ID)),
         Platform = str_replace_all(Platform, " ", ""),
         Platform = ifelse(is.na(Platform), "Mixed", Platform),
         Platform = ifelse(sapply(str_split(Platform, ","), length) != 1, "Mixed", Platform))


# Spread of intracellular correlations per dataset, colouring text by species

species_col <- ifelse(plot_df$Species == "Human", "royalblue", "goldenrod")

p1 <- 
  ggplot(plot_df) +
  geom_segment(aes(y = ID, x = `1st Qu.`, xend = `3rd Qu.`)) +
  geom_point(aes(y = ID, x = Median), shape = 21, colour = "black", fill = "firebrick", size = 3.4) +
  xlab("Cell to cell Pearson's correlation") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = species_col, size = 15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        plot.margin = margin(10, 30, 10, 10))


ggsave(p1, dpi = 300, width = 9, height = 12, device = "png",
       filename = file.path(plot_dir, "microglia_cell_to_cell_cor.png"))


# Scatter of median intracellular cor with UMI and count of cells

p2a <- 
  ggplot(plot_df, aes(x = log10(N_cells), y = Median)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = lm, colour = "darkblue") +
  xlab("Log10 count of cells") +
  ylab("Median intra-cellular correlation") +
  theme_classic() +
  theme(text = element_text(size = 20))



p2b <- 
  ggplot(plot_df, aes(x = log10(Median_UMI), y = Median)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = lm, colour = "darkblue") +
  xlab("Log10 median UMI/UMI-like") +
  ylab("Median intra-cellular correlation") +
  theme_classic() +
  theme(text = element_text(size = 20))
p2b_nolab <- p2b + theme(axis.title.y = element_blank())


p2 <- ggarrange(p2a, p2b_nolab, nrow = 1)

ggsave(p2, dpi = 300, width = 12, height = 5, device = "png",
       filename = file.path(plot_dir, "microglia_median_intracell_cor_scatter.png"))


# Boxplot of median intracellular cor with platform

p3 <- 
  ggplot(plot_df, aes(x = reorder(Platform, Median), y = Median)) +
  geom_boxplot() +
  geom_jitter(shape = 21, fill = "slategrey", width = 0.1) +
  xlab(NULL) +
  ylab("Median intra-cellular correlation") +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(p3, dpi = 300, width = 12, height = 5, device = "png",
       filename = file.path(plot_dir, "microglia_median_intracell_cor_platform_boxplot.png"))
