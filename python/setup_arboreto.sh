#!/bin/bash

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
source ~/miniconda3/bin/activate
conda init --all


conda create --name agg
conda activate agg
conda install -c bioconda arboreto

# https://github.com/aertslab/pySCENIC/issues/525
# https://github.com/aertslab/arboreto/issues/42
# Rollback installed dask (distributed?) to prevent localcluster issue
conda create -n pyscenic python=3.10.15
pip install dask-expr==0.5.3 distributed==2024.2.1

# https://github.com/aertslab/pySCENIC/issues/564
# Needed to rollback scipy to get back attribute from CSC objects
pip install scipy==1.13


# Needed to 're-install' (?) dask to get slurm + jobqueue
conda install dask  # actually not sure if this is needed
conda install dask-jobqueue -c conda-forge

# For loading .RDS objects
pip install rds2py