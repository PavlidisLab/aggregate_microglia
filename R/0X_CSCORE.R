library(tidyverse)
library(CSCORE)
library(Seurat)
library(parallel)
source("R/utils/functions.R")
source("/home/amorin/Projects/TR_singlecell/R/utils/vector_comparison_functions.R")
source("R/00_config.R")

# library(devtools)
# install_github("ChangSuBiostats/CS-CORE")


mcg_meta <- read.delim(mcg_meta_dedup_path)
ens_hg <- read.delim(ens_hg_path)
ens_mm <- read.delim(ens_mm_path)

cxg_hg <- filter(mcg_meta, Species == "Human" & str_detect(Data_link, "cellxgene"))
cxg_mm <- filter(mcg_meta, Species == "Mouse" & str_detect(Data_link, "cellxgene"))

multi <- sapply(c(cxg_hg$ID, cxg_mm$ID), function(x) {
  list.files(file.path(paste0("/home/amorin/sc_input/", x)), pattern = "seurat")
})

single <- multi[sapply(multi, length) == 1]
multi <- multi[setdiff(names(multi), names(single))]


outfile <- "/space/scratch/amorin/R_objects/jan2025_cxg_cscore.RDS"



ids <- names(single)




run_all <- function(id) {

  message(paste(id, Sys.time()))
  
  tryCatch({
    
    species <- filter(mcg_meta, ID == id)$Species
    
    if (species == "Mouse") {
      ens <- ens_mm
    } else {
      ens <- ens_hg
    }
    
    # Load cxg seurat object
    path <- file.path(sc_dir, id, paste0(id, "_cellxgene_seurat.RDS"))
    dat <- readRDS(path)
    
    # Ensembl -> refseq
    counts <- dat@assays$RNA@counts
    keep_ens <- intersect(rownames(counts), ens$Gene_ID)
    counts <- counts[keep_ens, ]
    rownames(counts) <- ens$Symbol[match(keep_ens, ens$Gene_ID)]
    
    dat <- CreateSeuratObject(counts = counts, meta.data = dat@meta.data)
    
    # Set raw counts
    DefaultAssay(dat) <- "RNA"
    
    
    # Microglia cells with raw counts, CPM, and log norm
    mcg_dat <- dat[, dat$cell_type == "microglial cell"]
    mcg_dat_cpm <- mcg_dat_ln <- mcg_dat
    mcg_dat_cpm <- NormalizeData(mcg_dat_cpm, normalization.method = "RC", scale.factor = 1e6)
    mcg_dat_ln <- NormalizeData(mcg_dat, normalization.method = "LogNormalize")
    
    
    # Ensure counts
    stopifnot(is.integer(mcg_dat[["RNA"]]$counts@i))
    
    
    # Top expressed genes for input
    avg_expr <- rowMeans(mcg_dat@assays$RNA$counts)
    keep_genes <- names(head(sort(avg_expr, decreasing = TRUE), 3000))
    
    # Run CSCORE
    res <- CSCORE(mcg_dat, genes = keep_genes)
    
    # Pearson's cor on raw counts, CPM, and log CPM
    counts <- t(mcg_dat[["RNA"]]$counts[keep_genes, ])
    cpm <- t(mcg_dat_cpm[["RNA"]]$data[keep_genes, ])
    ln <- t(mcg_dat_ln[["RNA"]]$data[keep_genes, ])
    
    cmat_count <- aggtools::sparse_pcor(counts)
    cmat_cpm <- aggtools::sparse_pcor(cpm)
    cmat_ln <- aggtools::sparse_pcor(ln)
    
    
    method_sim_list <- mclapply(keep_genes, function(x) {
      
      df <- data.frame(
        Symbol = keep_genes,
        CSCORE = res$est[keep_genes, x],
        Pcor_raw = cmat_count[keep_genes, x],
        Pcor_cpm = cmat_cpm[keep_genes, x],
        Pcor_ln = cmat_ln[keep_genes, x],
        Avg_expr = avg_expr[keep_genes],
        row.names = keep_genes
      )
      df <- filter(df, Symbol != x)
      
      scor <- cor(select_if(df, is.numeric), method = "spearman")
      colnames(scor) <- paste0("Scor_", colnames(scor))
      
      topk <- colwise_topk_intersect(as.matrix(select_if(df, is.numeric)), k = 200)
      colnames(topk) <- paste0("Topk_", colnames(topk))
      
      c(scor["CSCORE", 2:5], topk["CSCORE", 2:5])
      
    }, mc.cores = ncore)
    
    
    method_sim_df <- data.frame(
      Symbol = keep_genes, 
      do.call(rbind, method_sim_list),
      Avg_exp = avg_expr[keep_genes]
    )
    
    return(list(
      Summ_df = method_sim_df,
      CSCORE = res,
      Cor_cpm = cmat_cpm
    ))
    
  }, 
  error = function(e) NA)  
}



if (!file.exists(outfile)) {
  res_l <- mclapply(ids, run_all, mc.cores = ncore)
  names(res_l) <- ids
  saveRDS(res_l, outfile)
}



res_l <- readRDS(outfile)
names(res_l) <- ids
res_l <- res_l[!is.na(res_l)]



tt1 <- sapply(res_l, function(x) {
  median(x$Summ_df$Scor_Pcor_cpm, na.rm = TRUE)
}) %>% 
  data.frame(Median_scor = .) %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(mcg_meta, by = "ID")


tt2 <- plot(log10(tt1$N_cells), tt1$Median_scor)


scor_dens_l <- lapply(names(res_l), function(x) {
  
  ggplot(res_l[[x]]$Summ_df, aes(x = Scor_Pcor_cpm)) +
    geom_density() +
    xlab("Spearman's cor") +
    ggtitle(x) +
    theme_classic() +
    theme(text = element_text(size = 20),
          title = element_text(size = 15))
})




topk_dens_l <- lapply(names(res_l), function(x) {
  
  ggplot(res_l[[x]]$Summ_df, aes(x = Topk_Pcor_cpm)) +
    geom_density() +
    xlab("Top200 agreement") +
    ggtitle(x) +
    theme_classic() +
    theme(text = element_text(size = 20),
          title = element_text(size = 15))
})


egg::ggarrange(plots = scor_dens_l, nrow = 4)
egg::ggarrange(plots = topk_dens_l, nrow = 4)




# Load cxg seurat object
id <- "GSE207848"
ngenes <- 5000
path <- file.path(sc_dir, id, paste0(id, "_cellxgene_seurat.RDS"))
dat <- readRDS(path)

# Ensembl -> refseq
counts <- dat@assays$RNA@counts
keep_ens <- intersect(rownames(counts), ens_mm$Gene_ID)
counts <- counts[keep_ens, ]
rownames(counts) <- ens_mm$Symbol[match(keep_ens, ens_mm$Gene_ID)]

dat <- CreateSeuratObject(counts = counts, meta.data = dat@meta.data)

# Set raw counts
DefaultAssay(dat) <- "RNA"


# Microglia cells with raw counts, CPM, and log norm
mcg_dat <- dat[, dat$cell_type == "microglial cell"]
mcg_dat_cpm <- mcg_dat_ln <- mcg_dat
mcg_dat_cpm <- NormalizeData(mcg_dat_cpm, normalization.method = "RC", scale.factor = 1e6)
mcg_dat_ln <- NormalizeData(mcg_dat, normalization.method = "LogNormalize")


# Ensure counts
stopifnot(is.integer(mcg_dat[["RNA"]]$counts@i))


# Top expressed genes for input
avg_expr <- rowMeans(mcg_dat@assays$RNA$counts)
keep_genes <- names(head(sort(avg_expr, decreasing = TRUE), ngenes))

# Run CSCORE
res <- CSCORE(mcg_dat, genes = keep_genes)

# Pearson's cor on raw counts, CPM, and log CPM
counts <- t(mcg_dat[["RNA"]]$counts[keep_genes, ])
cpm <- t(mcg_dat_cpm[["RNA"]]$data[keep_genes, ])
ln <- t(mcg_dat_ln[["RNA"]]$data[keep_genes, ])

cmat_count <- aggtools::sparse_pcor(counts)
cmat_cpm <- aggtools::sparse_pcor(cpm)
cmat_ln <- aggtools::sparse_pcor(ln)



method_sim_list <- mclapply(keep_genes, function(x) {
  
  df <- data.frame(
    Symbol = keep_genes,
    CSCORE = res$est[keep_genes, x],
    Pcor_raw = cmat_count[keep_genes, x],
    Pcor_cpm = cmat_cpm[keep_genes, x],
    Pcor_ln = cmat_ln[keep_genes, x],
    Avg_expr = avg_expr[keep_genes],
    row.names = keep_genes
  )
  df <- filter(df, Symbol != x)
  
  scor <- cor(select_if(df, is.numeric), method = "spearman")
  colnames(scor) <- paste0("Scor_", colnames(scor))
  
  topk <- colwise_topk_intersect(as.matrix(select_if(df, is.numeric)), k = 200)
  colnames(topk) <- paste0("Topk_", colnames(topk))
  
  c(scor["CSCORE", 2:5], topk["CSCORE", 2:5])
  
}, mc.cores = ncore)


method_sim_df <- data.frame(
  Symbol = keep_genes, 
  do.call(rbind, method_sim_list),
  Avg_exp = avg_expr[keep_genes]
)


summary(method_sim_df)
cor(select_if(method_sim_df, is.numeric), method = "spearman", use = "pairwise.complete.obs")



ggplot(method_sim_df, aes(x = Scor_Pcor_cpm, y = Topk_Pcor_cpm)) +
  geom_point(shape = 21, alpha = 0.6, size = 3.2) +
  ylab("Top200 agreement") +
  xlab("Spearman's correlation") +
  ggtitle(id) +
  theme_classic() +
  theme(text = element_text(size = 25))


ggplot(method_sim_df, aes(x = Avg_exp, y = Scor_Pcor_cpm)) +
  geom_point(shape = 21, alpha = 0.6, size = 3.2) +
  ylab("Spearman's correlation") +
  xlab("Average counts") +
  ggtitle(id) +
  theme_classic() +
  theme(text = element_text(size = 25))




check_gene <- "Polr1d"


df <- data.frame(
  Symbol = keep_genes,
  CSCORE = res$est[keep_genes, check_gene],
  Pcor_raw = cmat_count[keep_genes, check_gene],
  Pcor_cpm = cmat_cpm[keep_genes, check_gene],
  Pcor_ln = cmat_ln[keep_genes, check_gene],
  Avg_expr = avg_expr[keep_genes],
  row.names = keep_genes
)
df <- filter(df, Symbol != check_gene)



plot(df$CSCORE, df$Pcor_raw, main = paste(id, check_gene))
plot(df$CSCORE, df$Pcor_cpm, main = paste(id, check_gene))
plot(df$CSCORE, df$Pcor_ln, main = paste(id, check_gene))

plot(df$Avg_expr, df$CSCORE, xlab = "Average raw counts", ylab = "CSCORE", main = paste(id, check_gene), cex.lab = 2)


gene1 <- check_gene
gene2 <- "Pabpn1"
plot(counts[, gene1], counts[, gene2], xlab = paste(gene1, "counts"), ylab = paste(gene2, "counts"), main = paste(id, gene1, gene2), cex.lab = 2)
plot(cpm[, gene1], cpm[, gene2],xlab = paste(gene1, "CPM"), ylab = paste(gene2, "CPM"), main = paste(id, gene1, gene2), cex.lab = 2)
plot(ln[, gene1], ln[, gene2], xlab = paste(gene1, "LN"), ylab = paste(gene2, "LN"), main = paste(id, gene1, gene2), cex.lab = 2)
