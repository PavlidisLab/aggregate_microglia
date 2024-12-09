## Isolated download of Sepp2023 mouse cerebellum data
## https://cellxgene.cziscience.com/collections/72d37bc9-76cc-442d-9131-da0e273862db

library(tidyverse)
library(Seurat)
library(aggtools)
library(pheatmap)
source("R/00_config.R")
source("R/utils/vector_comparison_functions.R")
source("/home/amorin/Projects/TR_singlecell/R/utils/functions.R")
source("/home/amorin/Projects/TR_singlecell/R/utils/plot_functions.R")

# Ensembl protein coding table
pc <- read.delim("/home/amorin/Data/Metadata/ensembl_mouse_protein_coding_105.tsv")

# Existing evidence
trsc <- readRDS("/space/scratch/amorin/R_objects/ranking_agg_integrated_TF_mm.RDS")
gr2023_rank <- readRDS("/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS")
gr2023_dat <- readRDS("/space/scratch/amorin/R_objects/all_data_list_Apr2022.RDS")
curated_all_path <- "/home/amorin/Data/Metadata/Curated_targets_all_Sept2023.tsv"
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)


# ID and input/output
id <- "Sepp2023"
rds_link <- "https://datasets.cellxgene.cziscience.com/18fb432a-be41-4677-9396-e1680a0852bf.rds"
rds_path <- file.path("/home/amorin/sc_input", id, paste0(id, "_cellxgene_seurat.RDS"))
out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))


if (!file.exists(rds_path)) {
  dir.create(file.path("/home/amorin/sc_input", id))
  download.file(rds_link, rds_path)
}


if (!file.exists(processed_path)) {
  
  dir.create(out_dir)
  
  dat <- readRDS(rds_path)
  
  # Discrete counts -- counts slot is empty
  mat <- GetAssayData(dat, layer = "data")  
  
  # Isolate meta and add counts
  meta <- dat[[]] %>% 
    dplyr::rename(Cell_type = cell_type) %>% 
    rownames_to_column(var = "ID") %>% 
    add_count_info(mat = mat)
  
  # Remove cells failing QC, keep only protein coding genes, and normalize
  mat <- rm_low_qc_cells(mat, meta) %>%
    ensembl_to_symbol(ensembl_df = pc) %>% 
    get_pcoding_only(pcoding_df = pc) %>% 
    Seurat::NormalizeData(., normalization.method = "RC", scale.factor = 1e6, verbose = FALSE)
  
  saveRDS(list(Mat = mat, Meta = meta), file = processed_path)
  rm(dat)
  
}


dat <- readRDS(processed_path)
# mat <- dat$Mat
mat <- Matrix::Matrix(log(dat$Mat + 1), sparse = TRUE)
meta <- dat$Meta


change_from_gc <- c("Theiler stage 18",
                    "Theiler stage 20",
                    "Theiler stage 21",
                    "Theiler stage 22")

change_to_gc <- "Theiler stage 18-22"


order_gc <- c("Theiler_stage_18_22",
              "Theiler_stage_23",
              "Theiler_stage_25",
              "Theiler_stage_27",
              "4_day_old_stage",
              "7_day_old_stage",
              "2_week_old_stage",
              "2_month_old_stage")


change_from_gcp <- c("Theiler stage 21",
                     "Theiler stage 22")

change_to_gcp <- "Theiler stage 21-22"


rm_gcp <- c("Theiler stage 18",
            "Theiler stage 20",
            "2-month-old stage",
            "2-week-old stage")


order_gcp <- c("Theiler_stage_21_22",
               "Theiler_stage_23",
               "Theiler_stage_25",
               "Theiler_stage_27",
               "4_day_old_stage",
               "7_day_old_stage")



meta_gc <- meta %>% 
  filter(Cell_type == "cerebellar granule cell") %>% 
  dplyr::select(-Cell_type) %>% 
  dplyr::rename(Cell_type = development_stage) %>% 
  mutate(
    Cell_type = as.character(Cell_type),
    Cell_type = ifelse(Cell_type %in% change_from_gc, change_to_gc, Cell_type),
    Cell_type = str_replace_all(Cell_type, " |-", "_"))


meta_gcp <- meta %>% 
  filter(Cell_type == "cerebellar granule cell precursor") %>% 
  dplyr::select(-Cell_type) %>% 
  dplyr::rename(Cell_type = development_stage) %>% 
  filter(Cell_type %!in% rm_gcp) %>% 
  mutate(
    Cell_type = as.character(Cell_type),
    Cell_type = ifelse(Cell_type %in% change_from_gcp, change_to_gcp, Cell_type),
    Cell_type = str_replace_all(Cell_type, " |-", "_"))
  
  


generate_cmat_list <- function(mat, meta, pc_df, ncore) {
  
  cts <- unique(meta$Cell_type)
  
  cmat_l <- mclapply(cts, function(ct) {
    
    tryCatch({
      
      ct_mat <- prepare_celltype_mat(mat = mat,
                                     meta = meta,
                                     pc_df = pc_df,
                                     cell_type = ct)
      
      cmat <- calc_sparse_correlation(ct_mat, cor_method = "pearson")
      cmat[cmat > 1] <- 1
      cmat[is.na(cmat)] <- 0
      fzmat <- fisherz(cmat)
      diag(fzmat) <- NA
      
      fzmat
  
    })

    
  }, mc.cores = ncore)
  
  names(cmat_l) <- cts
  
  return(cmat_l)
}




cor_gc <- generate_cmat_list(mat = mat, meta = meta_gc, pc_df = pc, ncore = ncore)
cor_gcp <- generate_cmat_list(mat = mat, meta = meta_gcp, pc_df = pc, ncore = ncore)


pax6_gc <- as.matrix(bind_cols(lapply(cor_gc, function(x) x[, "Pax6"])))
rownames(pax6_gc) <- pc$Symbol
pax6_gc <- pax6_gc[, order_gc]


pax6_gcp <- as.matrix(bind_cols(lapply(cor_gcp, function(x) x[, "Pax6"])))
rownames(pax6_gcp) <- pc$Symbol
pax6_gcp <- pax6_gcp[, order_gcp]





pax6_curated <- curated %>% 
  filter(str_to_lower(TF_Symbol) == "pax6")

pax6_targets <- unique(pax6_curated$Target_Symbol)

dplyr::count(pax6_curated, Target_Symbol, sort = TRUE)
dplyr::count(pax6_curated, Experiment_Type, sort = TRUE)
table(pax6_curated$Target_Symbol, pax6_curated$Experiment_Type)



summ_df <- data.frame(
  Symbol = pc$Symbol,
  Cor_GC = rowMeans(pax6_gc, na.rm = TRUE),
  SD_GC = matrixStats::rowSds(pax6_gc),
  N_msr_GC = rowSums(pax6_gc != 0),
  Cor_GCP = rowMeans(pax6_gcp, na.rm = TRUE),
  SD_GCP = matrixStats::rowSds(pax6_gcp),
  N_msr_GCP = rowSums(pax6_gcp != 0)
)

summ_df$Cor_RP <- rank(rank(-summ_df$Cor_GC) * rank(-summ_df$Cor_GCP))
summ_df$Curated_target <- summ_df$Symbol %in% str_to_title(pax6_curated)

summ_df_filt <- filter(summ_df, (N_msr_GC != 0) | (N_msr_GCP != 0))
summ_df_filt$Cor_RP <- rank(rank(-summ_df_filt$Cor_GC) * rank(-summ_df_filt$Cor_GCP))




ggplot(summ_df, aes(x = Cor_GC, y = Cor_GCP)) +
  geom_point(shape = 21) +
  geom_hline(yintercept = 0, col = "firebrick", linetype = "dashed") +
  geom_vline(xintercept = 0, col = "firebrick", linetype = "dashed") +
  xlab("Average FZ(Pcor) Granule Cells") +
  ylab("Average FZ(Pcor) Granule Cells Progenitors") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))



plot_genes <- slice_min(summ_df, Cor_RP, n = 50)$Symbol
plot_mat <- cbind(pax6_gcp[plot_genes, ], pax6_gc[plot_genes, ])

cols <- gplots::bluered(100)
breaks <- seq(-max(plot_mat), max(plot_mat), length.out = 100)


pheatmap(plot_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = cols,
         breaks = breaks,
         border_color = NA,
         fontsize = 15,
         gaps_col = rep(length(order_gcp), 5))



# Cor ~ curation status
sum(summ_df$Curated_target, na.rm = TRUE)
boxplot(summ_df$Cor_GC ~ summ_df$Curated_target)
boxplot(summ_df$Cor_GCP ~ summ_df$Curated_target)

sum(summ_df_filt$Curated_target, na.rm = TRUE)
boxplot(summ_df_filt$Cor_GC ~ summ_df_filt$Curated_target)
boxplot(summ_df_filt$Cor_GCP ~ summ_df_filt$Curated_target)



# gr2023$Mouse$Pax6 %>% view


summ_df_gr2023 <- summ_df %>% 
  left_join(., gr2023$Mouse$Pax6, by = "Symbol") %>%
  filter(Symbol != "Pax6")

summ_df_gr2023_filt <- filter(summ_df_gr2023, Symbol %in% summ_df_filt$Symbol)

# Cor ~ binding for all
plot(summ_df_gr2023$Cor_GC, summ_df_gr2023$Mean_bind)
plot(summ_df_gr2023$Cor_GCP, summ_df_gr2023$Mean_bind)

# Cor ~ binding for filtered
plot(summ_df_gr2023_filt$Cor_GC, summ_df_gr2023_filt$Mean_bind)
plot(summ_df_gr2023_filt$Cor_GCP, summ_df_gr2023_filt$Mean_bind)

# Cor ~ perturb for all
plot(summ_df_gr2023$Count_DE, summ_df_gr2023$Cor_GC)
plot(summ_df_gr2023$Count_DE, summ_df_gr2023$Cor_GCP)

# Cor ~ perturb for filtered
plot(summ_df_gr2023_filt$Count_DE, summ_df_gr2023_filt$Cor_GC)
plot(summ_df_gr2023_filt$Count_DE, summ_df_gr2023_filt$Cor_GCP)





# trsc$Pax6 %>% view

summ_df_trsc <- summ_df %>% 
  left_join(., trsc$Pax6, by = "Symbol") %>% 
  filter(Symbol != "Pax6")

summ_df_trsc_filt <- filter(summ_df_trsc, Symbol %in% summ_df_filt$Symbol)

# Cor ~ agg coexpr for all
plot(summ_df_trsc$Cor_GC, summ_df_trsc$Avg_aggr_coexpr)
plot(summ_df_trsc$Cor_GCP, summ_df_trsc$Avg_aggr_coexpr)

# Cor ~ agg coexpr for filt
plot(summ_df_trsc_filt$Cor_GC, summ_df_trsc_filt$Avg_aggr_coexpr)
plot(summ_df_trsc_filt$Cor_GCP, summ_df_trsc_filt$Avg_aggr_coexpr)



# Inspecting Pax6 consistency

# Between time points within cell type
k <- 200
topk_mat <- pax6_gcp  # pax6_gc
topk_mat <- topk_mat[setdiff(rownames(topk_mat), "Pax6"), ]
topk_mat <- colwise_topk_intersect(topk_mat, k = k)
topk_df <- mat_to_df(topk_mat, value_name = "Topk")

diag(topk_mat) <- NA
pheatmap(topk_mat)


# Between time points across cell type
topk_mat <- cbind(pax6_gcp, pax6_gc)
topk_mat <- topk_mat[setdiff(rownames(topk_mat), "Pax6"), ]
colnames(topk_mat) <- c(paste0("GCP_", colnames(pax6_gcp)), paste0("GC_", colnames(pax6_gc)))
topk_mat <- colwise_topk_intersect(topk_mat, k = k)
topk_df <- mat_to_df(topk_mat, value_name = "Topk")

diag(topk_mat) <- NA
pheatmap(topk_mat)



# Between averaged cell types

topk_intersect(
  topk_sort(rowMeans(pax6_gc, na.rm = TRUE), k = k),
  topk_sort(rowMeans(pax6_gcp, na.rm = TRUE), k = k)
)



# For all measured genes

topk_l <- mclapply(pc$Symbol, function(gene) {
  
  gc <- as.matrix(bind_cols(lapply(cor_gc, function(mat) mat[, gene])))
  gcp <- as.matrix(bind_cols(lapply(cor_gcp, function(mat) mat[, gene])))
  rownames(gc) <- rownames(gcp) <- pc$Symbol
  
  if (sum(gc, na.rm = TRUE) == 0 || sum(gcp, na.rm = TRUE) == 0) {
    return(NA)
  }
  
  topk_intersect(
    topk_sort(rowMeans(gc, na.rm = TRUE), k = k),
    topk_sort(rowMeans(gcp, na.rm = TRUE), k = k)
  )
  
  
}, mc.cores = ncore)
names(topk_l) <- pc$Symbol
topk_l <- topk_l[!is.na(topk_l)]


view(data.frame(unlist(topk_l)))



view(data.frame(
  GC = rowMeans(gc, na.rm = TRUE),
  GCP = rowMeans(gcp, na.rm = TRUE))
)
