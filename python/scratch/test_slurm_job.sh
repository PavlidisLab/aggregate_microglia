#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --output=test_slurm.out

source /etc/profile
source ~/.bashrc

conda activate agg
python3 /home/amorin/Projects/aggregate_microglia/R/misc/test_arboreto_grnboost2_TFonly.py