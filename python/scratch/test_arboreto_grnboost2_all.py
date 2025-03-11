import pandas as pd
import os.path
import time
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from arboreto.algo import genie3
from distributed import Client, LocalCluster


if __name__ == '__main__':

    tfs_path = "/home/amorin/Data/Metadata/TFs_human.tsv"
    mat_path = "/space/scratch/amorin/R_objects/GSE180928_mcg_filt.tsv"

    tfs = pd.read_table(tfs_path)["Symbol"].tolist()
    mat = pd.read_table(mat_path)
    
    tfs = set(tfs).intersection(mat.columns)

    local_cluster = LocalCluster(n_workers=8, threads_per_worker=1)
    custom_client = Client(local_cluster)

    start = time.time()

    network = grnboost2(expression_data=mat, 
                        tf_names=None,
                        seed=5,
                        client_or_address=custom_client)
    
    end = time.time()
    
    print(end - start)

    network.to_csv("/space/scratch/amorin/R_objects/GRNBoost2_all.tsv", 
                   sep="\t", header=False, index=False)
    
    custom_client.close()
    local_cluster.close()
