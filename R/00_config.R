## Specify pathing and global vars
## -----------------------------------------------------------------------------


# cores for parallel::
ncore <- 8

# A gene must be 'measured' in at least this fraction of datasets to keep
min_msr_frac <- 1/3


# Main data output
data_out_dir <- "/space/scratch/amorin/aggregate_microglia"


# Main plot output
plot_dir <- "/home/amorin/Plots/aggregate_microglia"


# Coexpression and aggregates matrices root directory
cmat_dir <- file.path(data_out_dir, "Cormats")
cmat_dir_mcg_hg <- file.path(cmat_dir, "Microglia_hg")
cmat_dir_mcg_mm <- file.path(cmat_dir, "Microglia_mm")
cmat_dir_macro_hg <- file.path(cmat_dir, "Macrophage_hg")
cmat_dir_macro_mm <- file.path(cmat_dir, "Macrophage_mm")


if (!dir.exists(data_out_dir)) dir.create(data_out_dir)
if (!dir.exists(plot_dir)) dir.create(plot_dir)
if (!dir.exists(cmat_dir)) dir.create(cmat_dir)
if (!dir.exists(cmat_dir_mcg_hg)) dir.create(cmat_dir_mcg_hg)
if (!dir.exists(cmat_dir_mcg_mm)) dir.create(cmat_dir_mcg_mm)
if (!dir.exists(cmat_dir_macro_hg)) dir.create(cmat_dir_macro_hg)
if (!dir.exists(cmat_dir_macro_mm)) dir.create(cmat_dir_macro_mm)



### !!!NOTE!!!: These paths are from external projects


## --- From the TRsc paper

# TRsc metadata
sc_meta_path <- "/space/grp/amorin/Metadata/single_cell_dataset_meta.tsv"

# Directory of initial objects to process scRNA-seq data
sc_dir <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets"


# Location of aggregate coexpression matrices
amat_dir <- "/space/scratch/amorin/TR_singlecell"

# List of cell types per dataset
celltype_list_path <- "/space/scratch/amorin/TRsc_output/celltype_list.RDS"

# ENSEMBL and Refseq Select protein coding paths
ref_hg_path <- "/space/grp/amorin/Metadata/refseq_select_hg38_jan2024.tsv"
ref_mm_path <- "/space/grp/amorin/Metadata/refseq_select_mm10_jan2024.tsv"
ens_hg_path <- "/space/grp/amorin/Metadata/ensembl_human_protein_coding_105.tsv"
ens_mm_path <- "/space/grp/amorin/Metadata/ensembl_mouse_protein_coding_105.tsv"

# AnimalTFDB paths for TF genes
tfs_hg_path <-  "/space/grp/amorin/Metadata/AnimalTFDB_human_V4.tsv"
tfs_mm_path <-  "/space/grp/amorin/Metadata/AnimalTFDB_mouse_V4.tsv"

# Orthology mapping
pc_ortho_path <- "/space/grp/amorin/Metadata/hg_mm_1to1_ortho_genes_DIOPT_V9.tsv"


## --- From Unibind analysis

# Unibind ChIP-seq experiments as GR objects
bind_gr_path_hg <- "/space/scratch/amorin/R_objects/unibind_grlist_perm_human.RDS"
bind_gr_path_mm <- "/space/scratch/amorin/R_objects/unibind_grlist_perm_mouse.RDS"

# Scored bind matrices and isolated metadata
bind_dat_path <- "/space/scratch/amorin/R_objects/processed_unibind_data.RDS"
bind_meta_path <- "/space/scratch/amorin/R_objects/Unibind_metadata.RDS"


## ---



# Cell type metadata paths
mcg_meta_path <- file.path(data_out_dir, "microglia_metadata.tsv")
mcg_meta_dedup_path <- file.path(data_out_dir, "microglia_metadata_dedup.tsv")
macro_meta_path <- file.path(data_out_dir, "macrophage_metadata.tsv")
macro_meta_dedup_path <- file.path(data_out_dir, "macrophage_metadata_dedup.tsv")
neuron_meta_path <- "/space/scratch/amorin/aggregate_microglia/neuron_metadata.tsv"
neuron_meta_dedup_path <- "/space/scratch/amorin/aggregate_microglia/neuron_metadata_dedup.tsv"
astro_meta_path <- "/space/scratch/amorin/aggregate_microglia/astrocyte_metadata.tsv"
astro_meta_dedup_path <- "/space/scratch/amorin/aggregate_microglia/astrocyte_metadata_dedup.tsv"


# Microglia list of count matrices and meta
mcg_dat_path <- file.path(data_out_dir, "microglia_dat_list.RDS")
macro_dat_path <- file.path(data_out_dir, "macrophage_dat_list.RDS")
neuron_dat_path <- file.path(data_out_dir, "neuron_dat_list.RDS")
astro_dat_path <- file.path(data_out_dir, "astrocyte_dat_list.RDS")


# Cell correlations per dataset
cell_cor_path <- file.path(data_out_dir, "microglia_cell_cors.RDS")

                            
# Summary of gene counts across experiments
mcg_count_summ_list_path <- file.path(data_out_dir, "microglia_gene_count_summary_list.RDS")
mcg_count_summ_table_hg <- file.path(data_out_dir, "microglia_gene_count_summary_table_hg.tsv")
mcg_count_summ_table_mm <- file.path(data_out_dir, "microglia_gene_count_summary_table_mm.tsv")
macro_count_summ_list_path <- file.path(data_out_dir, "macrophage_gene_count_summary_list.RDS")


# Lists of aggregate coexpression and NA tracking matrices
# Microglia
mcg_fz_hg_path <- file.path(cmat_dir_mcg_hg, "aggregate_cormat_FZ_microglia_hg.RDS")
mcg_allrank_hg_path <- file.path(cmat_dir_mcg_hg, "aggregate_cormat_allrank_microglia_hg.RDS")
mcg_colrank_hg_path <- file.path(cmat_dir_mcg_hg, "aggregate_cormat_colrank_microglia_hg.RDS")
mcg_fz_mm_path <- file.path(cmat_dir_mcg_mm, "aggregate_cormat_FZ_microglia_mm.RDS")
mcg_allrank_mm_path <- file.path(cmat_dir_mcg_mm, "aggregate_cormat_allrank_microglia_mm.RDS")
mcg_colrank_mm_path <- file.path(cmat_dir_mcg_mm, "aggregate_cormat_colrank_microglia_mm.RDS")
mcg_allrank_filt_hg_path <- file.path(cmat_dir_mcg_hg, "aggregate_cormat_allrank_filter_microglia_hg.RDS")
mcg_allrank_filt_mm_path <- file.path(cmat_dir_mcg_mm, "aggregate_cormat_allrank_filter_microglia_mm.RDS")
mcg_allrank_scor_mm_path <- file.path(cmat_dir_mcg_mm, "aggregate_cormat_allrank_scor_microglia_mm.RDS")
mcg_fz_scor_mm_path <- file.path(cmat_dir_mcg_mm, "aggregate_cormat_FZ_scor_microglia_mm.RDS")
# Macrophage
macro_fz_hg_path <- file.path(cmat_dir_macro_hg, "aggregate_cormat_FZ_macrophage_hg.RDS")
macro_allrank_hg_path <- file.path(cmat_dir_macro_hg, "aggregate_cormat_allrank_macrophage_hg.RDS")
macro_fz_mm_path <- file.path(cmat_dir_macro_mm, "aggregate_cormat_FZ_macrophage_mm.RDS")
macro_allrank_mm_path <- file.path(cmat_dir_macro_mm, "aggregate_cormat_allrank_macrophage_mm.RDS")


# CSCORE aggregate output
cscore_hg_path <- file.path(data_out_dir, "CSCORE_aggregate_microglia_hg.RDS")
cscore_mm_path <- file.path(data_out_dir, "CSCORE_aggregate_microglia_mm.RDS")


# GRNBoost2 aggregate output
grn_avg_hg_path <- file.path(data_out_dir, "GRNBoost2_average_mat_hg.RDS")
grn_avg_mm_path <- file.path(data_out_dir, "GRNBoost2_average_mat_mm.RDS")
grn_allrank_hg_path <- file.path(data_out_dir, "GRNBoost2_allrank_mat_hg.RDS")
grn_allrank_mm_path <- file.path(data_out_dir, "GRNBoost2_allrank_mat_mm.RDS")


# Microglia cCREs with region to gene predictions
mcg_ccre_path_hg <- file.path(data_out_dir, "human_microglia_ccre_gene_cors.txt")
mcg_ccre_path_mm <- file.path(data_out_dir, "mouse_microglia_ccre_gene_cors.txt")
