## Organizing metadata and the underlying count matrices for microglia
## -----------------------------------------------------------------------------

library(tidyverse)
library(aggtools)
source("R/00_config.R")
source("R/utils/functions.R")

# Table of all assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# List of cell types
ct_list <- readRDS(celltype_list_path)

# Based on manual inspection of the cell type list
ct_str <- ".*microglia.*|^microg.*|^mgs$|^mg$|^micro$"

# Create a metadata table of the relevant datasets
ct_df <- create_celltype_df(pattern = ct_str, ct_list, sc_meta)

# Ensure pathing to the corresponding data objects
ct_df$Path <- file.path(amat_dir, ct_df$ID, paste0(ct_df$ID, "_clean_mat_and_meta_CPM.RDS"))
stopifnot(all(map_lgl(ct_df$Path, file.exists)))


# Load data and filter to cell type, then save out

save_function_results(
  path = mcg_dat_path,
  fun = load_and_filter_celltype,
  args = list(input_df = ct_df)
)


dat_l <- readRDS(mcg_dat_path)


# Add gene count data to cell type metadata
ct_df <- add_counts_to_meta(dat_l, ct_df)


# Create meta that collapses when an ID has multiple cell types
ct_df_dedup <- collapse_dupl_ids(ct_df)


# Tally of annotated cell types (including from IDs with multiple cell types)
n_ct <- dplyr::count(ct_df, Cell_type) 


# Tally of species for unique/dedup IDs
n_species <- dplyr::count(ct_df_dedup, Species)



# Saving out celltype metadata

write.table(ct_df, sep = "\t", quote = FALSE, row.names = FALSE, 
            file = mcg_meta_path)

write.table(ct_df_dedup, sep = "\t", quote = FALSE, row.names = FALSE, 
            file = mcg_meta_dedup_path)



# Plots


p1 <- ggplot(ct_df_dedup, aes(x = log10(N_cells), fill = Species)) +
  geom_histogram(bins = 20) +
  ggtitle("Microglia") +
  xlab("Log10 count of microglial cells") +
  ylab("Frequency") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = c(0.9, 0.5),
        plot.margin = margin(c(10, 20, 10, 10)))



p2 <- ggplot(ct_df_dedup, 
             aes(x = N_msr_prefilt, y = N_msr_postfilt, fill = Species)) +
  geom_point(shape = 21, size = 3.4, colour = "black") +
  ggtitle("Count of measured microglial genes") +
  xlab("Pre-filter count") +
  ylab("Post-filter count") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = c(0.9, 0.5),
        plot.margin = margin(c(10, 20, 10, 10)))




p3 <- ggplot(ct_df_dedup, aes(x = log10(Median_UMI), fill = Species)) +
  geom_histogram(bins = 20) +
  ggtitle("Median UMI/UMI-like") +
  xlab("Log10 count") +
  ylab("Frequency") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = c(0.9, 0.5),
        plot.margin = margin(c(10, 20, 10, 10)))
