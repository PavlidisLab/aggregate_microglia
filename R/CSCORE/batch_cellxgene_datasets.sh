#!/bin/bash

# Running CSCORE over all cell x gene datasets

input_file=$1
ncore=4
preprocess="/home/amorin/Projects/aggregate_microglia/R/CSCORE/cellxgene_datasets.R"

cat "$input_file" | xargs -n 1 -P "$ncore" Rscript "$preprocess"

