## Take cCRE-gene connections in microglia from catlas, and use the aggregate
## coexpression and binding data to link TF binding, ccRE status, and coexpr.
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(GenomicRanges)
library(ggrepel)
source("R/00_config.R")
source("R/utils/vector_comparison_functions.R")

# Gene tables
tfs_hg <- read.delim(tfs_hg_path)
tfs_mm <- read.delim(tfs_mm_path)
pc_ortho <- read.delim(pc_ortho_path)

# Aggregate microglia coexpression: focus on FZ
agg_hg <- readRDS(mcg_fz_hg_path)
agg_mm <- readRDS(mcg_fz_mm_path)

# Microglia cCREs
ccre_hg <- read.delim(mcg_ccre_path_hg)
ccre_mm <- read.delim(mcg_ccre_path_mm)

# Unibind chip-seq experiments as GR objects
bind_gr_hg <- readRDS(bind_gr_path_hg)
bind_gr_mm <- readRDS(bind_gr_path_mm)

# Bind metadata
bind_meta <- readRDS(bind_meta_path)

# List of gene count measurement summaries
count_summ <- readRDS(mcg_count_summ_list_path)



# Setup
# ------------------------------------------------------------------------------


# Keep only cCRE-associated genes also found in filtered aggregate
common_hg <- intersect(ccre_hg$Symbol, count_summ$Human$Filter_genes)
common_mm <- intersect(ccre_mm$Symbol, count_summ$Mouse$Filter_genes)

# For inspecting genes found in cCRE set only
ccre_only_hg <- setdiff(ccre_hg$Symbol, count_summ$Human$Filter_genes)
ccre_only_mm <- setdiff(ccre_mm$Symbol, count_summ$Mouse$Filter_genes)

# TFs that are measured in microglia
keep_tfs_hg <- intersect(tfs_hg$Symbol, count_summ$Human$Filter_genes)
keep_tfs_mm <- intersect(tfs_mm$Symbol, count_summ$Mouse$Filter_genes)

# The same cCRE ID can be associated to multiple genes AND the same gene more
# than once (differing start sites). Enforce that a cCRE can be associated to 
# gene only once
ccre_filt_hg <- ccre_hg %>% 
  filter(Symbol %in% common_hg) %>% 
  distinct(Symbol, cCRE_ID, .keep_all = TRUE)

ccre_filt_mm <- ccre_mm %>% 
  filter(Symbol %in% common_mm) %>% 
  distinct(Symbol, cCRE_ID, .keep_all = TRUE)

# Which cCRE genes exist in both species
ccre_ortho <- filter(pc_ortho, Symbol_hg %in% common_hg & Symbol_mm %in% common_mm)
tfs_ortho <- filter(pc_ortho, Symbol_hg %in% keep_tfs_hg & Symbol_mm %in% keep_tfs_mm)

# Ensure bind GR objects align with metadata (for extracting TF symbols)
stopifnot(identical(names(bind_gr_hg), bind_meta$Permissive_hg$File))
stopifnot(identical(names(bind_gr_mm), bind_meta$Permissive_mm$File))



# Associating each cCRE gene with its top coexpressed TFs.
# ------------------------------------------------------------------------------


# Binary matrix of whether or not a TF was in the top 10 coexpr. TFs per cCRE genes

build_top_tf_matrix <- function(genes, tfs, agg_mat, k = 10) {
  
  agg_mat <- agg_mat[tfs, ]
  agg_mat[is.na(agg_mat) | is.infinite(agg_mat)] <- 0
  
  top_tf_l <- lapply(genes, function(gene) {
    gene_vec <- agg_mat[, gene]
    top_tfs <- topk_sort(gene_vec, k)
    ifelse(tfs %in% top_tfs, 1, 0)
  })
  
  top_tf_mat <- do.call(cbind, top_tf_l)
  rownames(top_tf_mat) <- tfs
  colnames(top_tf_mat) <- genes
  
  return(top_tf_mat)
}


top_tfs_mat_hg <- build_top_tf_matrix(genes = common_hg, tfs = keep_tfs_hg, agg_mat = agg_hg$Agg_mat)
top_tfs_mat_mm <- build_top_tf_matrix(genes = common_mm, tfs = keep_tfs_mm, agg_mat = agg_mm$Agg_mat)


# Most common TFs that are top partners across cCRE genes
order_tfs_hg <- sort(rowSums(top_tfs_mat_hg), decreasing = TRUE)
order_tfs_mm <- sort(rowSums(top_tfs_mat_mm), decreasing = TRUE)


# Remove TFs that were never in the topK
never_tfs_hg <- names(which(order_tfs_hg == 0))
never_tfs_mm <- names(which(order_tfs_mm == 0))

top_tfs_mat_hg <- top_tfs_mat_hg[setdiff(names(order_tfs_hg), never_tfs_hg), ]
top_tfs_mat_mm <- top_tfs_mat_mm[setdiff(names(order_tfs_mm), never_tfs_mm), ]



# cCRE as Granges. Proximal == gene promoter, distal == assumed enhancer.
# ------------------------------------------------------------------------------


ccre_gr_hg <- ccre_filt_hg %>% 
  dplyr::rename(
    Chromosome = Chromosome_distal,
    Start = Start_distal,
    End = End_distal) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


ccre_gr_mm <- ccre_filt_mm %>% 
  dplyr::rename(
    Chromosome = Chromosome_distal,
    Start = Start_distal,
    End = End_distal) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)



# Summarize distance between midpoints of promoter/distal links
ccre_link_distance <- function(ccre) {
  
  # Distal
  ccre_distal <- ccre %>% 
    dplyr::rename(
      Chromosome = Chromosome_distal,
      Start = Start_distal,
      End = End_distal) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # Proximal
  ccre_proximal <- ccre %>% 
    dplyr::rename(
      Chromosome = Chromosome_proximal,
      Start = Start_proximal,
      End = End_proximal) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  dist <- distance(ccre_distal, ccre_proximal)
  names(dist) <- ccre$Symbol
  
  return(dist)
}


link_dist_hg <- ccre_link_distance(ccre_filt_hg)
link_dist_mm <- ccre_link_distance(ccre_filt_mm)



# Counting number of cCRE-gene associations
n_ccre_genes_hg <- n_distinct(ccre_filt_hg$Symbol)
n_ccre_genes_mm <- n_distinct(ccre_filt_mm$Symbol)
n_ccre_per_gene_hg <- dplyr::count(ccre_filt_hg, Symbol) %>% arrange(desc(n))
n_ccre_per_gene_mm <- dplyr::count(ccre_filt_mm, Symbol) %>% arrange(desc(n))
n_gene_per_ccre_hg <- dplyr::count(ccre_filt_hg, cCRE_ID) %>% arrange(desc(n))
n_gene_per_ccre_mm <- dplyr::count(ccre_filt_mm, cCRE_ID) %>% arrange(desc(n))



# Overlap ChIP-seq and cCREs
# ------------------------------------------------------------------------------


# Focus on SPI1: mouse has microglia data, human does not

tf_hg <- "SPI1"
tf_mm <- "Spi1"


all_tf_ids_hg <- bind_meta$Permissive_hg %>% 
  filter(Symbol == tf_hg) %>% 
  pull(File)


mcg_tf_ids_mm <- bind_meta$Permissive_mm %>% 
  filter(str_to_title(Symbol) == tf_mm & str_detect(ID, ".*microglia.*")) %>% 
  pull(File)


all_tf_ids_mm <- bind_meta$Permissive_mm %>% 
  filter(str_to_title(Symbol) == tf_mm) %>% 
  pull(File)


nonmcg_tf_ids_mm <- setdiff(all_tf_ids_mm, mcg_tf_ids_mm)



tf_ol_hg <- findOverlaps(ccre_gr_hg, GRangesList(bind_gr_hg[all_tf_ids_hg]))
tf_ol_mm <- findOverlaps(ccre_gr_mm, GRangesList(bind_gr_mm[all_tf_ids_mm]))
tf_mcg_ol_mm <- findOverlaps(ccre_gr_mm, GRangesList(bind_gr_mm[mcg_tf_ids_mm]))
tf_nonmcg_ol_mm <- findOverlaps(ccre_gr_mm, GRangesList(bind_gr_mm[nonmcg_tf_ids_mm]))


# Get the cCRE-genes links that overlap each bind experiment
extract_ol_ccre_genes <- function(ol, ccre, ids) {
  
  ix_ccre <- as.integer(ol@from)
  ix_bind <- as.integer(ol@to)
  
  data.frame(Symbol = ccre$Symbol[ix_ccre], 
             cCRE = ccre$cCRE_ID[ix_ccre],
             Bind = ids[ix_bind]) 
  
}


# Summarize how many cCREs were bound for each gene
summarize_ccre_binding <- function(ol, ccre, ids) {
  
  ol_table <- extract_ol_ccre_genes(ol, ccre, ids)
  n_ccre <- dplyr::count(ccre, Symbol)
  
  summ_df <- ol_table %>% 
    group_by(Symbol) %>%
    dplyr::summarise(
      N_cCRE_bound = n_distinct(cCRE),
      N_total_ol = length(Bind)) %>% 
    left_join(n_ccre, by = "Symbol") %>% 
    dplyr::rename("N_cCRE_total" = "n") %>% 
    mutate(Avg_ol_per_cCRE = N_total_ol / N_cCRE_total) %>% 
    relocate(N_cCRE_total, .after = Symbol) %>% 
    mutate(Prop_cCRE_bound = N_cCRE_bound / N_cCRE_total)
  
  return(list(Overlap_table = ol_table, Summary = summ_df))
  
}


ol_summ_hg <- summarize_ccre_binding(tf_ol_hg, ccre_filt_hg, all_tf_ids_hg)
ol_summ_mm <- summarize_ccre_binding(tf_ol_mm, ccre_filt_mm, all_tf_ids_mm)
ol_summ_mcg_mm <- summarize_ccre_binding(tf_mcg_ol_mm, ccre_filt_mm, mcg_tf_ids_mm)


# cCRE genes w/o any overlap
no_ol_hg <- setdiff(common_hg, ol_summ_hg$Summary$Symbol)
no_ol_mm <- setdiff(common_mm, ol_summ_mm$Summary$Symbol)



# Relate bound cCRE-genes to agg coexpr
# ------------------------------------------------------------------------------


# Joins tables for provided TF and creates a grouping for cCRE, ortho, and 
# binding status

ready_all_df <- function(agg_mat, 
                         keep_genes, 
                         ccre_genes, 
                         ortho_ccre_genes, 
                         tf, 
                         ol_summ) {
  
  data.frame(
    Symbol = keep_genes,
    FZ = agg_mat[keep_genes, tf]
  ) %>%
    filter(Symbol != tf) %>%
    left_join(ol_summ, by = "Symbol") %>%
    mutate(
      cCRE_gene = Symbol %in% ccre_genes,
      cCRE_ortho_gene = Symbol %in% c(ortho_ccre_genes$Symbol_hg, 
                                      ortho_ccre_genes$Symbol_mm),
      Has_bind = !is.na(N_cCRE_bound),
      Group = case_when(
        cCRE_ortho_gene & Has_bind ~ "cCRE(+) Ortho(+) Bind(+)",
        cCRE_ortho_gene & !Has_bind ~ "cCRE(+) Ortho(+) Bind(-)",
        !cCRE_ortho_gene & Has_bind ~ "cCRE(+) Ortho(-) Bind(+)",
        !cCRE_ortho_gene & cCRE_gene & !Has_bind ~ "cCRE(+) Ortho(-) Bind(-)",
        .default = "cCRE(-)")
      ) %>%
    replace(is.na(.), 0)
}



# Using info from all human spi1 chip datasets
all_df_hg <- ready_all_df(
  agg_mat = agg_hg$Agg_mat,
  keep_genes = count_summ$Human$Filter_genes,
  ccre_genes = common_hg,
  ortho_ccre_genes = ccre_ortho,
  tf = tf_hg,
  ol_summ = ol_summ_hg$Summary
)


# Using info from all mouse spi chip datasets
all_df_mm <- ready_all_df(
  agg_mat = agg_mm$Agg_mat,
  keep_genes = count_summ$Mouse$Filter_genes,
  ccre_genes = common_mm,
  ortho_ccre_genes = ccre_ortho,
  tf = tf_mm,
  ol_summ = ol_summ_mm$Summary
)


# Using info from microglia mouse spi datasets
mcg_df_mm <- ready_all_df(
  agg_mat = agg_mm$Agg_mat,
  keep_genes = count_summ$Mouse$Filter_genes,
  ccre_genes = common_mm,
  ortho_ccre_genes = ccre_ortho,
  tf = tf_mm,
  ol_summ = ol_summ_mcg_mm$Summary
)



# Isolating the genes that had the strongest evidence

top_targets_all_hg <- all_df_hg %>% 
  filter(Group == "cCRE(+) Ortho(+) Bind(+)") %>% 
  arrange(desc(FZ)) 
rownames(top_targets_all_hg) <- top_targets_all_hg$Symbol

top_targets_all_mm <- all_df_mm %>% 
  filter(Group == "cCRE(+) Ortho(+) Bind(+)") %>% 
  arrange(desc(FZ))
rownames(top_targets_all_mm) <- top_targets_all_mm$Symbol

top_targets_mcg_mm <- mcg_df_mm %>% 
  filter(Group == "cCRE(+) Ortho(+) Bind(+)") %>% 
  arrange(desc(FZ))
rownames(top_targets_mcg_mm) <- top_targets_mcg_mm$Symbol


# Organizing top ortho targets and their respective spi1 aggregate coexpr

top_targets_ortho <- ccre_ortho %>% 
  filter(Symbol_hg %in% top_targets_all_hg$Symbol & 
         Symbol_mm %in% top_targets_mcg_mm$Symbol)


top_targets_ortho <- data.frame(
  Symbol = top_targets_ortho$Symbol_hg,
  FZ_hg = top_targets_all_hg[top_targets_ortho$Symbol_hg, "FZ"],
  FZ_mm = top_targets_all_mm[top_targets_ortho$Symbol_mm, "FZ"]
) %>% 
  mutate(
    FZ_ortho = rowMeans(.[, c("FZ_hg", "FZ_mm")]),
    Group = Symbol %in% slice_max(., FZ_ortho, n = 10)$Symbol
  )




# Plots
# ------------------------------------------------------------------------------


# Showing density of how many times a TF was in the top 10 of ccre genes
plot(density(order_tfs_hg[rownames(top_tfs_mat_hg)]), main = "Human", xlab = "Count of genes")
plot(density(order_tfs_mm[rownames(top_tfs_mat_mm)]), main = "Mouse", xlab = "Count of genes")


# Binary heatmap (ccre gene x tf) showing if a TF was top 10 coexpr for a ccre gene
pheatmap(top_tfs_mat_hg,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         col = c("white", "black"),
         show_rownames = FALSE,
         show_colnames = FALSE)


pheatmap(top_tfs_mat_mm,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         col = c("white", "black"),
         show_rownames = FALSE,
         show_colnames = FALSE)


# Distribution of log10 bp distance between promoter-distal cCREs
hist(log10(link_dist_hg + 1), breaks = 100, xlim = c(0, 7))
hist(log10(link_dist_mm + 1), breaks = 100, xlim = c(0, 7))


# Histograms of the counts of ccre-gene associations
ggplot(n_ccre_per_gene_hg, aes(x = n)) +
  geom_histogram(bins = 30) +
  xlab("Count of cCREs per gene") +
  ylab("Count of genes") +
  ggtitle("Human") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(n_ccre_per_gene_mm, aes(x = n)) +
  geom_histogram(bins = 30) +
  xlab("Count of cCREs per gene") +
  ylab("Count of genes") +
  ggtitle("Mouse") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(n_gene_per_ccre_hg, aes(x = n)) +
  geom_histogram(bins = 30) +
  xlab("Count of genes per cCRE") +
  ylab("Count of cCREs") +
  ggtitle("Human") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(n_gene_per_ccre_mm, aes(x = n)) +
  geom_histogram(bins = 30) +
  xlab("Count of genes per cCRE") +
  ylab("Count of cCREs") +
  ggtitle("Mouse") +
  theme_classic() +
  theme(text = element_text(size = 25))


# Scatterplots of the count of cCRE associations per genes and their expression
left_join(n_ccre_per_gene_hg, count_summ$Human$Summ_df, by = "Symbol") %>% 
  filter(Symbol %in% count_summ$Human$Filter_genes) %>% 
  ggplot(., aes(x = n, y = QN_avg)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Count of linked cCREs") +
  ylab("Mean log2 CPM") +
  ggtitle("Human") +
  theme_classic() +
  theme(text = element_text(size = 25))

left_join(n_ccre_per_gene_mm, count_summ$Mouse$Summ_df, by = "Symbol") %>% 
  filter(Symbol %in% count_summ$Mouse$Filter_genes) %>% 
  ggplot(., aes(x = n, y = QN_avg)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Count of linked cCREs") +
  ylab("Mean log2 CPM") +
  ggtitle("Mouse") +
  theme_classic() +
  theme(text = element_text(size = 25))


# Scatter plots of agg coexpr ~ average bound ccres
ggplot(all_df_hg, aes(x = Avg_ol_per_cCRE, y = FZ)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = "lm") +
  xlab("Average overlaps per cCRE") +
  ylab("Aggregate SPI1 coexpression") +
  ggtitle("Human") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(all_df_mm, aes(x = Avg_ol_per_cCRE, y = FZ)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = "lm") +
  xlab("Average overlaps per cCRE") +
  ylab("Aggregate SPI1 coexpression") +
  ggtitle("Mouse (all SPI ChIP-seq)") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(mcg_df_mm, aes(x = Avg_ol_per_cCRE, y = FZ)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = "lm") +
  xlab("Average overlaps per cCRE") +
  ylab("Aggregate SPI1 coexpression") +
  ggtitle("Mouse (microglia ChIP-seq)") +
  theme_classic() +
  theme(text = element_text(size = 25))



# Boxplots binarizing binding-ccre status
ggplot(all_df_hg, aes(x = Has_bind, y = FZ)) +
  geom_violin(fill = "slategrey") +
  geom_boxplot(width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab(">=1 SPI1-bound cCRE") +
  ylab("Aggregate SPI1 coexpression") +
  ggtitle("Human") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(all_df_mm, aes(x = cCRE_ortho_gene, y = FZ)) +
  geom_violin(fill = "slategrey") +
  geom_boxplot(width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Aggregate SPI1 coexpression") +
  ggtitle("Mouse (all ChIP-seq)") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(mcg_df_mm, aes(x = Group, y = FZ)) +
  geom_violin(fill = "slategrey") +
  geom_boxplot(width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Aggregate SPI1 coexpression") +
  ggtitle("Mouse (microglia ChIP-seq)") +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 18))


# Scatter of the aggr coexpr of top targets, labeling top in both species
px <- 
  ggplot(top_targets_ortho, aes(x = FZ_hg, y = FZ_mm)) +
  geom_point(shape = 21, size = 2.4, colour = "black", fill = "slategrey", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "firebrick") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "firebrick") +
  geom_text_repel(data = filter(top_targets_ortho, Group),
                  aes(label = Symbol),
                  size = 9) +
  ylab("Aggregate mouse coexpression (n=44)") +
  xlab("Aggregate human coexpression (n=18)") +
  theme_classic() +
  theme(text = element_text(size = 25))




ggsave(px, dpi = 300, device = "png", height = 8, width = 8,
       filename = file.path(plot_dir, "microglia_ccre_spi_coexpr.png"))
