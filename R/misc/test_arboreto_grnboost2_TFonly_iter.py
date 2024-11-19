import pandas as pd
from pathlib import Path
import time
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from distributed import Client, LocalCluster


if __name__ == '__main__':

    tfs_path = "/home/amorin/Data/Metadata/TFs_human.tsv"
    mat_path = "/space/scratch/amorin/R_objects/GSE180928_mcg_filt.tsv"
    out_dir = "/space/scratch/amorin/R_objects/arboreto_test"
    
    tfs = pd.read_table(tfs_path)["Symbol"].tolist()
    mat = pd.read_table(mat_path)

    Path(out_dir).mkdir(exist_ok=True)

    tfs = set(tfs).intersection(mat.columns)

    local_cluster = LocalCluster(n_workers=8, threads_per_worker=1)
    custom_client = Client(local_cluster)

    n_iter = 100

    start = time.time()

    for i in range(n_iter):
        
        network = grnboost2(expression_data=mat, 
                            tf_names=tfs,
                            seed=i,
                            client_or_address=custom_client)

        file_name = f"GRNBoost2_TFonly_iter{i}.tsv"
        file_path = Path(out_dir, file_name)
        network.to_csv(file_path, sep="\t", header=False, index=False)


    end = time.time()
    print(end - start)

    custom_client.close()
    local_cluster.close()
