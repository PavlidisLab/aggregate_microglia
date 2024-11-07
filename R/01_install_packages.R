## Collecting packages used. Note that most would have been installed individually.
## -----------------------------------------------------------------------------


options(repos = "http://cran.rstudio.com/")


packages <- c(
  # "testthat",
  "parallel",
  "tidyverse",
  "BiocManager",
  # "remotes",
  "Seurat",
  # "umap",
  "pheatmap",
  # "RcppTOML",
  # "igraph",
  # "Hmisc",
  "Matrix",
  # "patchwork",
  "plyr",
  "DescTools",
  # "googlesheets4",
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



# BiocManager::install("WGCNA")
BiocManager::install("preprocessCore")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("ComplexHeatmap")

devtools::install_github('PavlidisLab/ermineR')
devtools::install_github("PavlidisLab/aggtools")
