# https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/

library(RcisTarget)

motif_rank <- importRankings("/space/scratch/amorin/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
motif_score <- importRankings("/space/scratch/amorin/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather")
data(motifAnnotations_hgnc)

tf <- "MEF2C"
gene_list <- list(geneSetName = names(head(sort(cormat[, tf], decreasing = TRUE), 100)))



# view(data.frame(motifRankings@rankings[, c("motifs", "ASCL1")]))

motifEnrichmentTable_wGenes <- cisTarget(gene_list, 
                                         motif_rank,
                                         motifAnnot = motifAnnotations)


# Inspecting top TF/motifs for a given gene
check_gene <- "DLL1"

df <- left_join(
  motif_rank@rankings[, c("motifs", check_gene)],
  motif_score@rankings[, c("motifs", check_gene)],
  by = "motifs",
  suffix = c("_rank", "_score")
) %>%
  left_join(motifAnnotations, by = c("motifs" = "motif")) %>%
  arrange(!!sym(paste0(check_gene, "_rank")))
