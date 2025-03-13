## Collecting packages used. Note that most would have been installed individually.
## -----------------------------------------------------------------------------


options(repos = "http://cran.rstudio.com/")


packages <- c(
  "devtools",
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
  "rJava",
  "viridis",
  "cluster"
)


installed_packages <- packages %in% rownames(installed.packages())


if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


BiocManager::install("preprocessCore")
BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")
BiocManager::install("GenomicRanges")


devtools::install_github("PavlidisLab/aggtools")
devtools::install_github("ChangSuBiostats/CS-CORE")
