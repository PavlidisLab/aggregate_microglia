# Checking "perfect" cors in clustered genes from 2 datasets.

library(aggtools)
source("R/00_config.R")

dat1 <- load_scdat("/space/scratch/amorin/TR_singlecell/GSE160512/GSE160512_clean_mat_and_meta_CPM.RDS")
dat2 <- load_scdat("/space/scratch/amorin/TR_singlecell/GSE212616/GSE212616_clean_mat_and_meta_CPM.RDS")
pc_df <- read.delim("/home/amorin/Data/Metadata/refseq_select_mm10.tsv")



keep_genes <- c("Rps24",   # Ribo as sanity check
                "Rps10",
                "Ugt1a5",
                "Ugt1a10",
                "Ugt1a9",
                "Ugt1a2",
                "Ugt1a1",
                "Ugt1a6b",
                "Ugt1a6a",
                "Ugt1a7c",
                "H3c13",
                "H3c14",
                "Pcdha1",
                "Pcdha2")



mat1 <- prepare_celltype_mat(mat = dat1$Mat,
                             meta = dat1$Meta,
                             pc_df = pc_df,
                             cell_type = "microglia")



mat2 <- prepare_celltype_mat(mat = dat2$Mat,
                             meta = dat2$Meta,
                             pc_df = pc_df,
                             cell_type = c("Microglia", "Microglia (Low-quality/activated)"))


mat1 <- mat1[, keep_genes]
mat2 <- mat2[, keep_genes]



# Create duplicated 
mat3 <- Matrix::rsparsematrix(2000, 4, density = 0.5)
colnames(mat3) <- paste0("Gene", 1:4)
mat3[, "Gene4"] <- mat3[, "Gene3"]



cmat1 <- sparse_pcor(mat1)
cmat2 <- sparse_pcor(mat2)
cmat3 <- sparse_pcor(mat3)


# So GSE160512 has multiple ostensible duplicates, some of which give perfect
# cor of 1, others that give "1" that has machine precision difference
cmat1["Ugt1a10", "Ugt1a9"]
cmat1["Pcdha1", "Pcdha2"]
identical(cmat1["Ugt1a10", "Ugt1a9"], 1)
identical(cmat1["Pcdha1", "Pcdha2"], 1)
diff1_ugt1a <- mat1[, "Ugt1a10"] - mat1[, "Ugt1a9"]
diff1_pcdha <- mat1[, "Pcdha1"] - mat1[, "Pcdha2"]
all(diff1_ugt1a == 0)
all(diff1_pcdha == 0)


# Duplicated gene counts in GSE212616 give exact cor of 1
cmat2["H3c13", "H3c14"]
identical(cmat2["H3c13", "H3c14"], 1)


# Sparse cor gives exact 1 for duplicate
cmat3["Gene3", "Gene4"]
identical(cmat3["Gene3", "Gene4"], 1)
