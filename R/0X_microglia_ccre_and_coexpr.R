## TODO 
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(GenomicRanges)
source("R/00_config.R")
source("R/utils/vector_comparison_functions.R")

tfs_hg <- read.delim(tfs_hg_path)
tfs_mm <- read.delim(tfs_mm_path)

agg_hg <- readRDS("~/mcgdat/Cormats/Hg_pcor/aggregate_cormat_FZ_hg.RDS")
agg_mm <- readRDS("~/mcgdat/Cormats/Mm_pcor/aggregate_cormat_FZ_mm.RDS")


ccre_hg <- read.delim(mcg_ccre_path_hg)
ccre_mm <- read.delim(mcg_ccre_path_mm)

# List of gene count measurement summaries
count_summ <- readRDS(mcg_count_summ_list_path)

# Keep only cCRE-associated genes also found in aggregate (from refseq protein coding)
ccre_only_hg <- setdiff(ccre_hg$Symbol, colnames(agg_hg$Agg_mat))
ccre_only_mm <- setdiff(ccre_mm$Symbol, colnames(agg_mm$Agg_mat))

common_hg <- intersect(ccre_hg$Symbol, colnames(agg_hg$Agg_mat))
common_mm <- intersect(ccre_mm$Symbol, colnames(agg_mm$Agg_mat))

ccre_hg <- filter(ccre_hg, Symbol %in% common_hg)
ccre_mm <- filter(ccre_mm, Symbol %in% common_mm)


keep_tfs_hg <- intersect(tfs_hg$Symbol, count_summ$Human$Filter_genes)
keep_tfs_mm <- intersect(tfs_mm$Symbol, count_summ$Mouse$Filter_genes)


# count(ccre_hg, Symbol) %>% arrange(desc(n))
# count(ccre_mm, Symbol) %>% arrange(desc(n))




# Binary matrix of whether or not a TF was in the top 10 coexpressed TFs per gene
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



# Most common TFs that are top partners
order_tfs_hg <- sort(rowSums(top_tfs_mat_hg), decreasing = TRUE)
order_tfs_mm <- sort(rowSums(top_tfs_mat_mm), decreasing = TRUE)


plot(density(order_tfs_hg))
plot(density(order_tfs_mm))


# Remove TFs that were never in the topK
never_tfs_hg <- names(which(order_tfs_hg == 0))
never_tfs_mm <- names(which(order_tfs_mm == 0))

top_tfs_mat_hg <- top_tfs_mat_hg[setdiff(names(order_tfs_hg), never_tfs_hg), ]
top_tfs_mat_mm <- top_tfs_mat_mm[setdiff(names(order_tfs_mm), never_tfs_mm), ]


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



# Or just use mat directly...

agg_mat_hg <- agg_hg$Agg_mat[keep_tfs_hg, common_hg]
agg_mat_mm <- agg_mm$Agg_mat[keep_tfs_mm, common_mm]

agg_mat_hg[is.na(agg_mat_hg) | is.infinite(agg_mat_hg)] <- 0
agg_mat_mm[is.na(agg_mat_mm) | is.infinite(agg_mat_mm)] <- 0


# Most commonly coexpressed TF
avg_tfs_hg <- sort(rowMeans(agg_mat_hg), decreasing = TRUE)
avg_tfs_mm <- sort(rowMeans(agg_mat_mm), decreasing = TRUE)


plot(order_tfs_hg[keep_tfs_hg], avg_tfs_hg[keep_tfs_hg])
plot(order_tfs_mm[keep_tfs_mm], avg_tfs_mm[keep_tfs_mm])


# cCRE as Granges. Proximal == gene promoter, distal == assumed enhancer.

ccre_gr_hg <- ccre_hg %>% 
  dplyr::rename(
    Chromosome = Chromosome_distal,
    Start = Start_distal,
    End = End_distal) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


ccre_gr_mm <- ccre_mm %>% 
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


link_dist_hg <- ccre_link_distance(ccre_hg)
link_dist_mm <- ccre_link_distance(ccre_mm)



hist(log10(link_dist_hg + 1), breaks = 100, xlim = c(0, 7))
hist(log10(link_dist_mm + 1), breaks = 100, xlim = c(0, 7))



# TODO: this workflow is too slow...
# Helper to extract the mid point of a GR, returned as a GR
# get_midpoint <- function(gr) {
#   chr <- seqnames(gr)
#   mid <- start(gr) + floor(width(gr) / 2)
#   GRanges(paste0(chr, ":", mid))
# }
# 
# 
# # Helper to isolate the proximal ranges, returned as a GR
# get_proximal <- function(gr) {
#   GRanges(paste0(gr$Chromosome_proximal, ":", 
#                  gr$Start_proximal, "-",
#                  gr$End_proximal))
# }
# 
# 
# # Summarize distance between midpoints of promoter/distal links
# ccre_link_distance <- function(ccre_gr) {
#   
#   dist_l <- lapply(1:length(ccre_gr), function(i) {
#     message(i)
#     distal <- get_midpoint(ccre_gr[i])
#     proximal <- get_midpoint(get_proximal(ccre_gr[i]))
#     distance(proximal, distal)
#   })
# }