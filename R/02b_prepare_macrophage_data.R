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


# For inspecting cell types
all_ct_df <- do.call(rbind, lapply(names(ct_list), function(x) {
  data.frame(ct_list[[x]]$Ct_count, ID = x)
}))

tally_ct <- count(all_ct_df, Cell_type)


# These cell types needed inspection to confirm if macrophages
check <- c(
  "Macro",  # keep: GSE134003, GSE217368, and GSE183310
  "MAC",  # keep: GSE222510
  "Mac",  # keep: GSE198582
  "Mac-I",  # keep: GSE198582 (Mac-I to V and Mac-a to e)
  "Mac_",  # keep: GSE198582
  "Mon/Mac",  # rm: unclear if distinct from monocytes
  "Monocyte/Macrophage",  # rm: unclear if distinct from monocytes
  "endo-Mac",  # keep: GSE198582
  "epi-Mac",  # keep: GSE198582
  "macrophage dendritic cell progenitor",  # rm: unclear if distinct from monocytes
  "Macrophages_progenitors",  # rm: unclear if distinct from monocytes
  "MonocyteMacrophages",  # rm: unclear if distinct from monocytes
  "prol.Mac",  # keep: GSE198582
  "recMac",  # keep: GSE225664 (recruited macrophage)
  "IMs",  # keep: GSE225667, GSE225664, GSE225662, GSE225666, GSE225663, GSE225665 (interstitial macrophages)
  "Macrophage/cDC"  # rm: unclear if distinct from conventional dendritic cells
)


check_df <- filter(all_ct_df, Cell_type %in% check)


# This was based off of manual inspection of the cell type list
ct_str <- ".*macrophage.*|.*kupffer.*|^macro$|^mac$|mac-.*|mac_.*|.*-mac$|^recmac$|^ims$"


# Remove these cell types
rm_str <- c("macrophage dendritic cell progenitor",
            "MonocyteMacrophages",
            "Mon/Mac",
            "Monocyte/Macrophage",
            "Macrophages_progenitors")


ct_df <- create_celltype_df(pattern = ct_str, ct_list, sc_meta)
ct_df <- filter(ct_df, Cell_type %!in% rm_str)


# Ensure pathing to the corresponding data objects
ct_df$Path <- file.path(amat_dir, ct_df$ID, paste0(ct_df$ID, "_clean_mat_and_meta_CPM.RDS"))
stopifnot(all(map_lgl(ct_df$Path, file.exists)))


# Load data and filter to cell type, then save out

save_function_results(
  path = macro_dat_path,
  fun = load_and_filter_celltype,
  args = list(input_df = ct_df)
)


dat_l <- readRDS(mcg_dat_path)


# Add gene count data to cell type metadata
ct_df <- add_counts_to_meta(dat_l, ct_df)


# Create meta that collapses when an ID has multiple cell types
ct_df_dedup <- collapse_dupl_ids(ct_df)


# Tally of annotated cell types (including from IDs with multiple cell types)
n_ct <- count(ct_df, Cell_type) 


# Tally of species for unique/dedup IDs
n_species <- count(ct_df_dedup, Species)



# Saving out celltype metadata

write.table(ct_df, sep = "\t", quote = FALSE, row.names = FALSE, 
            file = macro_meta_path)

write.table(ct_df_dedup, sep = "\t", quote = FALSE, row.names = FALSE, 
            file = macro_meta_dedup_path)



# Plots


p1 <- ggplot(ct_df_dedup, aes(x = log10(N_cells), fill = Species)) +
  geom_histogram(bins = 20) +
  ggtitle("Macrophage") +
  xlab("Log10 count of macrophgage cells") +
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
  ggtitle("Count of measured macrophage genes") +
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
