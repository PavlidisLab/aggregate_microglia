## Load single cell data and filter count matrix and meta to just microglia
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

mcg_meta <- read.delim(mcg_meta_path)


load_and_filter_dat <- function(input_df) {
  
  ids <- unique(input_df$ID)
  
  dat_l <- lapply(ids, function(id) {
    
    message(paste(id, Sys.time()))
    df <- filter(input_df, ID == id)
    path <- unique(df$Path)
    
    dat <- aggtools::load_scdat(path)
    meta <- filter(dat$Meta, Cell_type %in% df$Cell_type)
    mat <- dat$Mat[, meta$ID]
    list(Mat = mat, Meta = meta)
    
  })
  names(dat_l) <- ids
  
  return(dat_l)
}



# TODO: run if not saved logic
dat_l <- load_and_filter_dat(mcg_meta)

saveRDS(dat_l, mcg_dat_path)



# Get counts of measured genes per dataset. Pre filtering means any gene with
# at least one count in any cell is considered. Post filtering uses the min
# measurement filtering for calculating coexpression: gene must have at least
# one count in at least 20 cells

n_msr <- lapply(dat_l, function(x) {
  
  pre <- t(x$Mat)
  post <- zero_sparse_cols(t(x$Mat))
  
  data.frame(
    N_gene_prefilt = sum(colSums(pre) != 0),
    N_gene_postfilt = sum(colSums(post) != 0)
  )
  
})



# Combine with meta

mcg_meta <- do.call(rbind, n_msr) %>% 
  as.data.frame() %>% 
  mutate(ID = names(dat_l)) %>% 
  left_join(., mcg_meta, by = "ID") %>% 
  relocate(c(N_gene_prefilt, N_gene_postfilt), .after = N_cells)


write.table(mcg_meta, sep = "\t", quote = FALSE, row.names = FALSE, 
            file = mcg_meta_path)





# Datasets that have multiple microglial annotations (each gets collapsed)
dupl_ids <- unique(mcg_meta$ID[which(duplicated(mcg_meta$ID))])

n_celltype <- count(mcg_meta, Cell_type) # this includes duplicates


n_species <- mcg_meta %>% distinct(ID, .keep_all = TRUE) %>% count(Species)


# Plotting


p1 <- ggplot(mcg_meta, aes(x = log10(N_cells))) +
    geom_histogram(bins = 20, colour = "slategrey", fill = "slategrey") +
    ggtitle("Count of microglial cells") +
    xlab("Log10 count") +
    ylab("Frequency") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))



p2 <- ggplot(mcg_meta, 
             aes(x = N_gene_prefilt, y = N_gene_postfilt)) +
  geom_point(shape = 21, size = 3.4, fill = "slategrey", colour = "white") +
  ggtitle("Count of microglial genes") +
  xlab("Pre-filter count") +
  ylab("Post-filter count") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))

