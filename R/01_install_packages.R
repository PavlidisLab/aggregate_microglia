## Collecting packages used. Note that most would have been installed individually.
## -----------------------------------------------------------------------------


options(repos = "http://cran.rstudio.com/")


packages <- c(
  "testthat",
  "parallel",
  "tidyverse",
  "BiocManager",
  "Seurat",
  "umap",
  "pheatmap",
  "Matrix",
  "matrixStats",
  "plyr",
  "ROCR",
  "data.table",
  "ggrepel",
  "qlcMatrix",
  "microbenchmark",
  "credentials",
  "rJava"
)


installed_packages <- packages %in% rownames(installed.packages())


if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}



BiocManager::install("preprocessCore")
# BiocManager::install("WGCNA")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("ComplexHeatmap")

devtools::install_github('PavlidisLab/ermineR')
devtools::install_github("PavlidisLab/aggtools")
