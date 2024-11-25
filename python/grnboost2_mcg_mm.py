import pandas as pd
import numpy as np
from pathlib import Path
from arboreto.algo import grnboost2
from dask.distributed import Client, LocalCluster
import utils
import rds2py
import datetime

n_workers=20  # How many workers for dask client
n_iter=100  # How many times to run GRNBoost2

mcg_dat_path = "/space/scratch/amorin/aggregate_microglia/microglia_dat_list.RDS"
mcg_meta_path = "/space/scratch/amorin/aggregate_microglia/microglia_metadata_dedup.tsv"
pc_mm_path = "/home/amorin/Data/Metadata/refseq_select_mm10.tsv"
tfs_mm_path = "/home/amorin/Data/Metadata/TFs_mouse.tsv"
out_dir = "/space/scratch/amorin/aggregate_microglia/GRNBoost2"

Path(out_dir).mkdir(exist_ok=True)

mcg_dat = rds2py.read_rds(mcg_dat_path)
mcg_meta = pd.read_table(mcg_meta_path)
pc_mm = pd.read_table(pc_mm_path)
tfs_mm = pd.read_table(tfs_mm_path)

# Isolate human IDs to iterate over dat_list
ids_mm = mcg_meta[mcg_meta.Species == "Mouse"]["ID"].tolist()

# Init dask local client
local_cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)
custom_client = Client(local_cluster)


if __name__ == '__main__':

    utils.iter_grnboost2_over_list(dat_list=mcg_dat, 
                                   ids = ids_mm, 
                                   gene_df=pc_mm, 
                                   tf_df=tfs_mm, 
                                   n_iter=n_iter,
                                   out_root=out_dir,
                                   client=custom_client)

    # Close connections
    custom_client.close()
    local_cluster.close()
