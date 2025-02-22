## TODO
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

tally_ct <- dplyr::count(all_ct_df, Cell_type)


# These cell types needed inspection to confirm if neurons
check <- c()


check_df <- filter(all_ct_df, Cell_type %in% check)


# This was based off of manual inspection of the cell type list
ct_str <- ".*astro.*"


# Remove these cell types
rm_str <- c()


ct_df <- create_celltype_df(pattern = ct_str, ct_list, sc_meta)
ct_df <- filter(ct_df, Cell_type %!in% rm_str)


# Ensure pathing to the corresponding data objects
ct_df$Path <- file.path(amat_dir, ct_df$ID, paste0(ct_df$ID, "_clean_mat_and_meta_CPM.RDS"))
stopifnot(all(map_lgl(ct_df$Path, file.exists)))


# Load data and filter to cell type, then save out

save_function_results(
  path = astro_dat_path,
  fun = load_and_filter_celltype,
  args = list(input_df = ct_df)
)


dat_l <- readRDS(astro_dat_path)


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
            file = astro_meta_path)

write.table(ct_df_dedup, sep = "\t", quote = FALSE, row.names = FALSE, 
            file = astro_meta_dedup_path)
