library(tidyverse)
library(CSCORE)
library(Seurat)
source("R/utils/functions.R")
source("/home/amorin/Projects/TR_singlecell/R/utils/vector_comparison_functions.R")
source("R/00_config.R")

# library(devtools)
# install_github("ChangSuBiostats/CS-CORE")


mcg_meta <- read.delim(mcg_meta_dedup_path)
ens_mm <- read.delim(ens_mm_path)

cxg_hg <- filter(mcg_meta, Species == "Human" & str_detect(Data_link, "cellxgene"))
cxg_mm <- filter(mcg_meta, Species == "Mouse" & str_detect(Data_link, "cellxgene"))

multi <- sapply(c(cxg_hg$ID, cxg_mm$ID), function(x) {
  list.files(file.path(paste0("/home/amorin/sc_input/", x)), pattern = "seurat")
})

single <- multi[sapply(multi, length) == 1]
multi <- multi[setdiff(names(multi), names(single))]


dat <- readRDS("~/sc_input/GSE207848/GSE207848_cellxgene_seurat.RDS")

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
keep_genes <- names(head(sort(avg_expr, decreasing = TRUE), 5000))

# Run CSCORE
res <- CSCORE(mcg_dat, genes = keep_genes)

# Pearson's cor on raw counts, CPM, and log CPM
counts <- t(mcg_dat[["RNA"]]$counts[keep_genes, ])
cpm <- t(mcg_dat_cpm[["RNA"]]$data[keep_genes, ])
ln <- t(mcg_dat_ln[["RNA"]]$data[keep_genes, ])

cmat_count <- aggtools::sparse_pcor(counts)
cmat_cpm <- aggtools::sparse_pcor(cpm)
cmat_ln <- aggtools::sparse_pcor(ln)


check_gene <- "Runx1"


df <- data.frame(
  Symbol = keep_genes,
  CSCORE = res$est[keep_genes, check_gene],
  Cor_raw = cmat_count[keep_genes, check_gene],
  Cor_cpm = cmat_cpm[keep_genes, check_gene],
  Cor_ln = cmat_ln[keep_genes, check_gene],
  Avg_expr = avg_expr[keep_genes],
  row.names = keep_genes
)
df <- filter(df, Symbol != check_gene)



cor(select_if(df, is.numeric), method = "spearman")
colwise_topk_intersect(as.matrix(select_if(df, is.numeric)), k = 200)


plot(df$CSCORE, df$Cor_raw)
plot(df$CSCORE, df$Cor_cpm)
plot(df$CSCORE, df$Cor_ln)
plot(df$CSCORE, df$Avg_expr)

plot(df$Cor_raw, df$Cor_cpm)
plot(df$Cor_raw, df$Cor_ln)
plot(df$Cor_raw, df$Avg_expr)

plot(df$Cor_cpm, df$Cor_ln)
plot(df$Cor_cpm, df$Avg_expr)

plot(df$Cor_ln, df$Avg_expr)

