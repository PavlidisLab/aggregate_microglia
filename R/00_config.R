## Specify pathing and global vars
## -----------------------------------------------------------------------------


# cores for parallel::
ncore <- 8


# Main data output
data_out_dir <- "/space/scratch/amorin/aggregate_microglia"


# Main plot output
plot_dir <- "/home/amorin/Plots/aggregate_microglia"


# Coexpression and aggregates matrices to these directories
cmat_dir_hg <- file.path(data_out_dir, "Cormats", "Hg_pcor")
cmat_dir_mm <- file.path(data_out_dir, "Cormats", "Mm_pcor")


if (!dir.exists(data_out_dir)) dir.create(data_out_dir)
if (!dir.exists(plot_dir)) dir.create(plot_dir)
if (!dir.exists(cmat_dir_hg)) dir.create(cmat_dir_hg)
if (!dir.exists(cmat_dir_mm)) dir.create(cmat_dir_mm)



## TODO: These are paths from the TRsc paper, FIX dependency
## TODO: or, alternatively, Borealis data download link

# TRsc metadata
sc_meta_path <- "/space/grp/amorin/Metadata/single_cell_dataset_meta.tsv"

# Location of aggregate coexpression matrices
amat_dir <- "/space/scratch/amorin/TR_singlecell"

# List of cell types per dataset
celltype_list_path <- "/space/scratch/amorin/R_objects/celltype_list.RDS"

# ENSEMBL and Refseq Select protein coding paths
ref_hg_path <- "/space/grp/amorin/Metadata/refseq_select_hg38.tsv"
ref_mm_path <- "/space/grp/amorin/Metadata/refseq_select_mm10.tsv"
ens_hg_path <- "/space/grp/amorin/Metadata/ensembl_human_protein_coding_105.tsv"
ens_mm_path <- "/space/grp/amorin/Metadata/ensembl_mouse_protein_coding_105.tsv"


##



# Cell type metadata paths
mcg_meta_path <- file.path(data_out_dir, "microglia_metadata.tsv")
mcg_meta_dedup_path <- file.path(data_out_dir, "microglia_metadata_dedup.tsv")
macro_meta_path <-  file.path(data_out_dir, "macrophage_metadata.tsv")
macro_meta_dedup_path <-  file.path(data_out_dir, "macrophage_metadata_dedup.tsv")


# Microglia list of count matrices and meta
mcg_dat_path <-  file.path(data_out_dir, "microglia_dat_list.RDS")
macro_dat_path <-  file.path(data_out_dir, "macrophage_dat_list.RDS")


# Cell correlations per dataset
cell_cor_path <-  file.path(data_out_dir, "microglia_cell_cors.RDS")

                            
# Summary of gene counts across experiments
count_summ_path <-  file.path(data_out_dir, "microglia_gene_count_summaries.RDS")
count_summ_table_hg <-  file.path(data_out_dir, "microglia_gene_count_summary_hg.tsv")
count_summ_table_mm <-  file.path(data_out_dir, "microglia_gene_count_summary_mm.tsv")
