## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
library(egg)
source("R/00_config.R")
source("R/utils/functions.R")

mcg_meta <- read.delim(mcg_meta_dedup_path)
dat_l <- readRDS(mcg_dat_path)


# Note that I inspected taking cor after filtering away 0 genes, but saw
# trivial change to summary stats

# TODO: check cell count and subsample at threshold

calc_cell_cor <- function(dat_l) {
  
  cor_l <- lapply(names(dat_l), function(x) {
    
    message(paste(x, Sys.time()))
    
    result <- tryCatch({
      cell_cor <- aggtools::sparse_pcor(dat_l[[x]]$Mat)
      cell_cor <- mat_to_df(cell_cor, value_name = "Cor")
      summary(cell_cor$Cor)
    }, 
    error = function(e) NULL)
    
    return(result)
    
  })
  names(cor_l) <- names(dat_l)
  
  return(cor_l)
}



calc_sample_cell_cor <- function(dat_l, 
                                 # mcg_meta, 
                                 n_cell_threshold = 10e3, 
                                 n_sample_cells = 1e3,
                                 n_samples = 10) {
  
  cor_l <- lapply(names(dat_l), function(x) {
    
    mat <- dat_l[[x]]$Mat
    
    sample_l <- lapply(1:n_samples, function(i) {
      sample_cells <-  sample(colnames(mat), n_sample_cells, replace = FALSE)
      submat <- mat[, sample_cells]
      cell_cor <- aggtools::sparse_pcor(submat)
      cell_cor <- mat_to_df(cell_cor, value_name = "Cor")[["Cor"]]
    })
    
    summary(unlist(sample_l))

  })
  
  names(cor_l) <- names(dat_l)
  return(cor_l)
}





if (!file.exists(cell_cor_path)) {
  cor_l <- calc_cell_cor(dat_l)
  saveRDS(cor_l, cell_cor_path)
} else {
  cor_l <- readRDS(cell_cor_path)
}


stop()


# TODO: some elements failed in lapply call, but can be regen isolated...
# figure out how to do one pass

null_cor <- mcg_meta %>% 
  filter(ID %in% names(which(sapply(cor_l, function(x) is.null(x))))) %>% 
  arrange(N_cells)


rerun_cor_ids <- filter(null_cor, N_cells < 10e3) %>% pull(ID) %>% unique()
sample_cor_ids <- setdiff(null_cor$ID, rerun_cor_ids)


rerun_cor <- calc_cell_cor(dat_l[rerun_cor_ids])
sample_cor <- calc_sample_cell_cor(dat_l[sample_cor_ids])


input_l <- c(cor_l[!sapply(cor_l, is.null)], rerun_cor, sample_cor)


# Setting up summary df for plotting

plot_df <- do.call(rbind, input_l) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(., mcg_meta, by = "ID") %>% 
  arrange(Median) %>%
  mutate(ID = factor(ID, levels = unique(ID)),
         Platform = str_replace_all(Platform, " ", ""),
         Platform = ifelse(is.na(Platform), "Mixed", Platform),
         Platform = ifelse(sapply(str_split(Platform, ","), length) != 1, "Mixed", Platform))




# cor.test(plot_df$Median, plot_df$N_cells, method = "spearman")
# cor.test(plot_df$Median, plot_df$Median_UMI, method = "spearman")
# table(plot_df$Platform)



# Spread of intracellular correlations per dataset, colouring text by species

# TODO: This throws a warning, need to figure out native support
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
