library(tidyverse)

curation_path <- "/home/amorin/Data/Metadata/scRNA_metadata_PavLab_25-09-2024.tsv"
source("R/00_config.R")

sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)
curation <- read.delim(meta_path, stringsAsFactors = FALSE)

candidates <- curation %>% 
  filter(str_detect(`brain.`, ".*yes.*"),
         str_detect(`single.cell.`, ".*yes.*"),
         str_detect(`cell.type.annotations.`, ".*yes.*")) %>% 
  select(Acc, Taxa, PubMed, number_of_cells, platform., additional_supplementary_files.x)


new <- setdiff(candidates$Acc, sc_meta$ID)
common <- intersect(candidates$Acc, sc_meta$ID)


table(candidates$Taxa)
