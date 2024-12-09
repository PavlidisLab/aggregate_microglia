## Cleaning up/isolating the Pax6 curated targets
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)

pc_hg <- read.delim("/home/amorin/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
pc_mm <- read.delim("/home/amorin/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)
pc_ortho <- read.delim("/home/amorin/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", stringsAsFactors = FALSE)
curated <- read.delim("/home/amorin/Data/Metadata/Curated_targets_all_Sept2023.tsv", stringsAsFactors = FALSE)

# Normalize casing to get all Pax6 targets
pax6_curated_original <- filter(curated, str_to_lower(TF_Symbol) == "pax6")
pax6_curated <- pax6_curated_original

# Symbols with no match in either mouse or human protein coding
nomatch_species <- setdiff(pax6_curated$Target_Symbol, c(pc_hg$Symbol, pc_mm$Symbol))

# Symbols with no match in the 1:1 orthologous set
nomatch_ortho <- setdiff(pax6_curated$Target_Symbol, c(pc_ortho$Symbol_hg, pc_ortho$Symbol_mm))

# Isolate Pax6 target and match protein coding symbols to species 
pax6_targets <- unique(pax6_curated$Target_Symbol) %>%
  tibble(ID = .) %>%
  mutate(
    Symbol_hg = case_when(
      ID %in% pc_ortho$Symbol_hg ~ ID,
      ID %in% pc_ortho$Symbol_mm ~ pc_ortho$Symbol_hg[match(ID, pc_ortho$Symbol_mm)],
      ID %in% pc_hg$Symbol ~ ID,
      TRUE ~ NA_character_
    ),
    Symbol_mm = case_when(
      ID %in% pc_ortho$Symbol_mm ~ ID,
      ID %in% pc_ortho$Symbol_hg ~ pc_ortho$Symbol_mm[match(ID, pc_ortho$Symbol_hg)],
      ID %in% pc_mm$Symbol ~ ID,
      TRUE ~ NA_character_
    )
  )


# Remove unmatched IDs
pax6_targets <- filter(pax6_targets, !is.na(Symbol_hg) | !is.na(Symbol_mm))


# Some IDs will be "duplicated" as same symbol but different casing. Also
# instances where species-specific casing doesn't match symbol.
# Deduplicate ID on a case-insensitive basis, remove Mouse Fgf15 == Human FGF19,
# then coerce ID to human symbol if available
pax6_targets <- pax6_targets %>%
  mutate(ID_lower = tolower(ID)) %>%
  group_by(ID_lower) %>%
  slice(1) %>%
  ungroup() %>%
  select(-ID_lower) %>% 
  filter(ID != "Fgf15") %>% 
  mutate(ID = ifelse(!is.na(Symbol_hg), Symbol_hg, ID))


# Join IDs back with master table of curated targets
pax6_curated <- left_join(
  mutate(pax6_curated, Norm_symbol = str_to_lower(Target_Symbol)),
  mutate(pax6_targets, Norm_ID = str_to_lower(ID)),
  by = c("Norm_symbol" = "Norm_ID")
) %>% 
  dplyr::select(-Norm_symbol)


# Attempting to reduce the duplicated curations found between Chu2021 and Pavlab continued curation
pax6_curated <- pax6_curated %>% 
  mutate(
    Cell_Type = str_replace(Cell_Type, "UBERON_|UBERON:|CL_|CLO_|CL:|CLO:", ""),
    Exp_ID = paste(ID, Experiment_Type, Target_Species, Cell_Type, sep = "_")
    ) %>%
  distinct(Exp_ID, .keep_all = TRUE)




# Get counts/binary status of experiment type and species for each unique target
summarized_targets <- pax6_curated %>%
  group_by(ID) %>%
  summarize(
    Perturbation = any(Experiment_Type == "Perturbation", na.rm = TRUE),
    Reporter = any(Experiment_Type == "Reporter", na.rm = TRUE),
    Binding = any(Experiment_Type == "Binding", na.rm = TRUE),
    Mouse = any(Target_Species == "Mouse", na.rm = TRUE),
    Human = any(Target_Species == "Human", na.rm = TRUE),
    Total_Experiments = n()
  ) %>%
  mutate(across(Perturbation:Human, as.integer)) %>% 
  filter(!is.na(ID)) %>% 
  arrange(desc(Total_Experiments))



# Bar plot of total experiments per target
ggplot(summarized_targets, 
       aes(x = reorder(ID, -Total_Experiments), y = Total_Experiments)) +
  geom_bar(stat = "identity", colour = "black", fill = "slategrey") +
  xlab("Target genes") +
  ylab("Count of experiments") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 20))


# Binary heatmap of experiment type/species status for top targets
plot_cols <- c("Perturbation", "Reporter", "Binding", "Mouse", "Human")
plot_mat <- as.matrix(summarized_targets[1:25, plot_cols])
rownames(plot_mat) <- summarized_targets$ID[1:25]
plot_mat <- t(plot_mat)

pheatmap(plot_mat,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = c("white", "black"),
         border_color = NA,
         fontsize = 20)
