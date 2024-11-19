## Running GENIE3 and comparing to standard Pcor matrix
## https://bioconductor.org/packages/release/bioc/vignettes/GENIE3/inst/doc/GENIE3.html

library(tidyverse)
library(GENIE3)
library(doRNG)
source("R/00_config.R")
source("R/utils/functions.R")

set.seed(23)

tf_hg <- read.delim("~/Data/Metadata/TFs_human.tsv")
id <- "GSE180928"
dat_l <- readRDS(mcg_dat_path)
mat <- dat_l[[id]]$Mat

keep_genes <- which(rowSums(mat > 0) >= 20)
mat <- mat[keep_genes, ]
tf <- intersect(tf_hg$Symbol, rownames(mat))

mat <- as.matrix(mat)


start <- Sys.time()
print(start)

res <- GENIE3(mat, regulators = tf, nCores = ncore)

end <- Sys.time()
print(end - start)

saveRDS(res, "~/scratch/R_objects/GENIE3_microglia_GSE180928_TFonly.RDS")
