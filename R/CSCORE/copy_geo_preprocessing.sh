#!/bin/bash

input_file=$1
input_dir="/home/amorin/Projects/TR_singlecell/R/preprocessing_scripts/CPM"
output_dir="/home/amorin/Projects/aggregate_microglia/R/CSCORE/test"



mkdir -p "$output_dir"

while IFS=$'\t' read -r id; do

    # Input and output filenames
    input_filename="${input_dir}/${id}_cpm.R"
    output_filename="${output_dir}/${id}.R"

    # Read the input file
    content=$(<"$input_filename")

    # Apply the required modifications
    modified_content=$(echo "$content" | sed \
        -e '/if (!dir.exists(dat_dir)) dir.create(dat_dir)/a \
mcg_dat_meta <- read.delim(mcg_meta_dedup_path, stringsAsFactors = FALSE)\n\
ct <- mcg_dat_meta %>% filter(ID == id) %>% pull(Cell_type)\n\
outfile <- file.path(data_out_dir, "CSCORE", paste0(id, "_CSCORE.RDS"))' \
        -e 's/processed_path/outfile/' \
        -e 's/library(WGCNA)/library(CSCORE)/' \
        -e '/source("R\/utils\/plot_functions.R")/d' \
        -e '/out_dir <- file.path(amat_dir, id)/,/^$/d' \
        -e '/p1 <- all_hist(meta)/,/ggsave(p2.*QC_scatter.png"))/d' \
        -e '/get_pcoding_only(pcoding_df = pc)/,/Seurat::NormalizeData.*verbose = FALSE)/d' \
        -e '/message(paste("Count of cells:.*n_distinct(meta\$Cell_type)))/,/gc()/d' \
        -e '/else {/,/mat <- dat\$Mat/d' \
        -e '/fwrite(.*rsr_all\$Agg_mat.*)/,/file = namat_path/d')

    # Write the modified content to the output file
    echo "$modified_content" > "$output_filename"

done < "$input_file"




# while IFS=$'\t' read -r id; do
# 
#     in_file="${input_dir}/${id}_cpm.R"
#     out_file="${output_dir}/${id}.R"
# 
#     echo "Beginning $out_file"
# 
# done < "$input_file"