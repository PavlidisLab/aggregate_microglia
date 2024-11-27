# https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/
# https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/
# https://www.bioconductor.org/packages/release/bioc/vignettes/RcisTarget/inst/doc/RcisTarget_MainTutorial.html

library(RcisTarget)
library(tidyverse)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")

motif_rank_hg <- importRankings("/space/scratch/amorin/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
motif_score_hg <- importRankings("/space/scratch/amorin/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather")
motif_rank_mm <- importRankings("/space/scratch/amorin/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")

data(motifAnnotations_hgnc)
motif_anno_hg <- motifAnnotations; rm(motifAnnotations)
data(motifAnnotations_mgi)
motif_anno_mm <- motifAnnotations; rm(motifAnnotations)

tfs_hg <- read.delim("/home/amorin/Data/Metadata/TFs_human.tsv", stringsAsFactors = FALSE)
tfs_mm <- read.delim("/home/amorin/Data/Metadata/TFs_mouse.tsv", stringsAsFactors = FALSE)


tf_rank_hg <- readRDS("/space/scratch/amorin/R_objects/ranking_agg_integrated_TF_hg.RDS")
tf_rank_mm <- readRDS("/space/scratch/amorin/R_objects/ranking_agg_integrated_TF_mm.RDS")


mcg_rank_hg <- readRDS(file.path(cmat_dir_hg, "aggregate_cormat_FZ_hg.RDS"))
mcg_rank_mm <- readRDS(file.path(cmat_dir_mm, "aggregate_cormat_FZ_mm.RDS"))




rcistarget_over_tf_list <- function(rank_list,
                                    motif_rank,
                                    motif_anno,
                                    n = 200) {
  
  tfs <- names(rank_list)
  
  res_l <- lapply(tfs, function(tf) {
    
    
    tryCatch({
      
      coexpr_genes <- slice_min(rank_list[[tf]], Rank_aggr_coexpr, n = n)$Symbol
      bind_genes <- slice_min(rank_list[[tf]], Rank_bind, n = n)$Symbol
    
      res_coexpr <- cisTarget(coexpr_genes, 
                              motifRankings = motif_rank, 
                              motifAnnot = motif_anno)
      
      res_bind <- cisTarget(bind_genes, 
                            motifRankings = motif_rank, 
                            motifAnnot = motif_anno)
      
      list(Coexpr = res_coexpr, Bind = res_bind)
      
    }, 
    error = function(e) NULL)
    
  })
  names(res_l) <- tfs
  
  return(res_l)
}





rcistarget_over_mat <- function(mat,
                                tfs,
                                motif_rank, 
                                motif_anno,
                                n = 200) {
  
  res_l <- mclapply(tfs, function(tf) {
    
    pos_coexpr <- names(head(sort(mat[, tf, drop = TRUE], decreasing = TRUE), n+1))
    pos_coexpr <- setdiff(pos_coexpr, tf)
    neg_coexpr <- names(head(sort(mat[, tf, drop = TRUE]), n))
    
    if (length(pos_coexpr) <= 1) {
      return(NA)
    }
    

    res_pos <- cisTarget(pos_coexpr, 
                         motifRankings = motif_rank, 
                         motifAnnot = motif_anno)
    
    res_neg <- cisTarget(neg_coexpr,
                         motifRankings = motif_rank, 
                         motifAnnot = motif_anno)
    
    list(Pos_coexpr = res_pos, Neg_coexpr = res_neg)
    
  })
  names(res_l) <- tfs
  res_l <- res_l[!is.na(res_l)]
  
  return(res_l)
}



rcistarget_tfrank_hg_path <- "/space/scratch/amorin/R_objects/Rcistarget_tfrank_hg.RDS" 
rcistarget_tfrank_mm_path <- "/space/scratch/amorin/R_objects/Rcistarget_tfrank_mm.RDS" 
rcistarget_mcgtf_hg_path <- "/space/scratch/amorin/R_objects/Rcistarget_mcgtf_hg.RDS" 
rcistarget_mcgtf_mm_path <- "/space/scratch/amorin/R_objects/Rcistarget_mcgtf_mm.RDS" 


# Human TF rank
save_function_results(
  path = rcistarget_tfrank_hg_path,
  fun = rcistarget_over_tf_list,
  args = list(
    rank_list = tf_rank_hg,
    motif_rank = motif_rank_hg, 
    motif_anno = motif_anno_hg,
    n = 200)
)


# Mouse TF rank
save_function_results(
  path = rcistarget_tfrank_mm_path,
  fun = rcistarget_over_tf_list,
  args = list(
    rank_list = tf_rank_mm,
    motif_rank = motif_rank_mm, 
    motif_anno = motif_anno_mm,
    n = 200)
)


# Human microglia TF
save_function_results(
  path = rcistarget_mcgtf_hg_path,
  fun = rcistarget_over_mat,
  args = list(
    mat = mcg_rank_hg$Agg_mat,
    tfs = tfs_hg$Symbol,
    motif_rank = motif_rank_hg, 
    motif_anno = motif_anno_hg,
    n = 200)
)



# Mouse microglia TF
save_function_results(
  path = rcistarget_mcgtf_mm_path,
  fun = rcistarget_over_mat,
  args = list(
    mat = mcg_rank_mm$Agg_mat,
    tfs = tfs_mm$Symbol,
    motif_rank = motif_rank_mm, 
    motif_anno = motif_anno_mm,
    n = 200)
)



rcistarget_tfrank_hg <- readRDS(rcistarget_tfrank_hg_path)
rcistarget_tfrank_mm <- readRDS(rcistarget_tfrank_mm_path)
rcistarget_mcgtf_hg <- readRDS(rcistarget_mcgtf_hg_path) 
rcistarget_mcgtf_mm <- readRDS(rcistarget_mcgtf_mm_path)


all_motifs_hg <- unique(motif_anno_hg$TF)
filter_motifs_hg <- unique(filter(motif_anno_hg, keptInRanking)$TF)


rcistarget_mcgtf_hg <- rcistarget_mcgtf_hg[!is.na(rcistarget_mcgtf_hg)]



tfrank_check_if_motif_present <- function(list) {
  
  tfs <- names(list)
  
  df_tfrank <- bind_rows(lapply(tfs, function(tf) {
    
    c(Pos_coexpr_high = any(str_detect(list[[tf]]$Coexpr$TF_highConf, tf)),
      Pos_coexpr_low = any(str_detect(list[[tf]]$Coexpr$TF_lowConf, tf)),
      Bind_high = any(str_detect(list[[tf]]$Bind$TF_highConf, tf)),
      Bind_low = any(str_detect(list[[tf]]$Bind$TF_lowConf, tf)))
  }))
  
  
  df_tfrank <- data.frame(Symbol = tfs, df_tfrank)
  
}


mcgrank_check_if_motif_present <- function(list) {
  
  tfs <- names(list)
  
  df_mcgtf <- bind_rows(lapply(tfs, function(tf) {
    
    c(Pos_high = any(str_detect(list[[tf]]$Pos_coexpr$TF_highConf, tf)),
      Pos_low = any(str_detect(list[[tf]]$Pos_coexpr$TF_lowConf, tf)),
      Neg_high = any(str_detect(list[[tf]]$Neg_coexpr$TF_highConf, tf)),
      Neg_low = any(str_detect(list[[tf]]$Neg_coexpr$TF_lowConf, tf)))
  }))
  
  
  df_mcgtf <- data.frame(Symbol = tfs, df_mcgtf)
  
}



df_tfrank_hg <- tfrank_check_if_motif_present(rcistarget_tfrank_hg)
df_tf_rank_hg <- filter(df_tf_rank_hg, Symbol %in% filter_motifs_hg)

df_tfrank_mm <- tfrank_check_if_motif_present(rcistarget_tfrank_mm)


sum(df_tfrank_hg$Pos_coexpr_high) / nrow(df_tfrank_hg)
sum(df_tfrank_hg$Bind_high) / nrow(df_tfrank_hg)

# sum(df_tfrank_hg$Pos_coexpr_high | df_tfrank_hg$Pos_coexpr_low) / nrow(df_tfrank_hg)
# sum(df_tfrank_hg$Bind_high | df_tfrank_hg$Bind_low) / nrow(df_tfrank_hg)


sum(df_tfrank_mm$Pos_coexpr_high) / nrow(df_tfrank_mm)
sum(df_tfrank_mm$Bind_high) / nrow(df_tfrank_mm)


df_mcgtf_hg <- mcgrank_check_if_motif_present(rcistarget_mcgtf_hg)
df_mcgtf_mm <- mcgrank_check_if_motif_present(rcistarget_mcgtf_hg)





sum(df_mcgtf_hg$Pos_high) / nrow(df_mcgtf_hg)
sum(df_mcgtf_hg$Neg_high) / nrow(df_mcgtf_hg)


sum(df_mcgtf_mm$Pos_high) / nrow(df_mcgtf_mm)
sum(df_mcgtf_mm$Neg_high) / nrow(df_mcgtf_mm)


# sum(df_mcgtf_hg$Pos_high | df_mcgtf_hg$Pos_low) / nrow(df_mcgtf_hg)
# sum(df_mcgtf_hg$Neg_high | df_mcgtf_hg$Neg_low) / nrow(df_mcgtf_hg)

view(df_mcgtf_hg[df_mcgtf_hg$Pos_high, ])







# Inspecting top TF/motifs for a given gene
check_gene <- "SPP1"

df <- left_join(
  motif_rank_hg@rankings[, c("motifs", check_gene)],
  motif_score_hg@rankings[, c("motifs", check_gene)],
  by = "motifs",
  suffix = c("_rank", "_score")
) %>%
  left_join(motif_anno_hg, by = c("motifs" = "motif")) %>%
  arrange(!!sym(paste0(check_gene, "_rank")))

