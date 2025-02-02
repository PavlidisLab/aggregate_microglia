#!/bin/bash

# Running CSCORE over all GEO datasets.

input_file=$1
ncore=4
dir="/home/amorin/Projects/aggregate_microglia/R/CSCORE/"

cat "$input_file" | xargs -P "$ncore" -I {} Rscript "$dir"{}.R
# cat "$input_file" | xargs -P "$ncore" -I {} echo "$dir"{}.R
