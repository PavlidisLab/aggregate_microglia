library(testthat)
source("R/utils/dev_functions.R")
source("R/00_config.R")


# Placeholder real data
id <- "GSE180928"
out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
pc <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
dat <- readRDS(processed_path)
meta <- dat$Meta
mat <- dat$Mat

cor_method = "sparse"
agg_method = "allrank"
min_cell = 20
ncores = ncore
sparse_arg <- (cor_method == "sparse")

amat <- init_agg_mat(mat)
na_mat <- amat  

ct <- "Endothelial"
ct_mat <- subset_celltype_and_filter(mat, meta, ct, min_cell, sparse = sparse_arg)
n_nonmsr <- max(sum(ct_mat == 0), sum(is.na(ct_mat)))
cmat1_raw <- calc_correlation(ct_mat, cor_method = cor_method, ncores = ncores)
na_mat <- increment_na_mat(cmat1_raw, na_mat)
cmat1_tr <- transform_correlation_mat(cmat1_raw, agg_method)
amat1 <- amat + cmat1_tr


ct <- "Microglia"
ct_mat <- subset_celltype_and_filter(mat, meta, ct, min_cell, sparse = sparse_arg)
n_nonmsr <- max(sum(ct_mat == 0), sum(is.na(ct_mat)))
cmat2_raw <- calc_correlation(ct_mat, cor_method = cor_method, ncores = ncores)
na_mat <- increment_na_mat(cmat2_raw, na_mat)
cmat2_tr <- transform_correlation_mat(cmat2_raw, agg_method)
amat2 <- amat1 + cmat2_tr


amat3 <- allrank_mat(amat2) / sum(!is.na(amat2))

amat3a <- allrank_mat(amat2, ties_arg = "max")
amat3b <- amat3a / sum(!is.na(amat3a))
amat3c <- diag_to_one(amat3b)
amat3d <- lowertri_to_symm(amat3c)

amat4a <- (cmat1_raw + cmat2_raw) / (2 - na_mat)
amat4b <- diag_to_one(amat4a)

plot(amat3d[, "RUNX1"], amat4b[, "RUNX1"])
view(data.frame(Allrank = amat3d[, "RUNX1"], Avg = amat4b[, "RUNX1"]))


# In raw, higher is more important
head(sort(cmat1_raw[, "RUNX1"], decreasing = TRUE), 10)

# In transformed, lower ranks are more important
head(sort(lowertri_to_symm(cmat1_tr)[, "RUNX1"]), 10)

# In the aggregate before finalization, lower sums are more important
head(sort(lowertri_to_symm(amat2)[, "RUNX1"]), 10)

# Finalization first re-ranks values, such that lower summed ranks (more important) now have a high rank
head(sort(lowertri_to_symm(amat3a)[, "RUNX1"], decreasing = TRUE), 10)

# Finally, it divides by the count of non-NA elements, which squeezes high rank (more important) to 1 and low rank to 0
head(sort(lowertri_to_symm(amat3b)[, "RUNX1"], decreasing = TRUE), 10)
