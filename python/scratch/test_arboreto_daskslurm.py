import pandas as pd
import time
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from dask.distributed import Client
from dask_jobqueue import SLURMCluster


if __name__ == '__main__':

    tfs_path = "/home/amorin/Data/Metadata/TFs_human.tsv"
    mat_path = "/space/scratch/amorin/R_objects/GSE180928_mcg_filt.tsv"

    tfs = pd.read_table(tfs_path)["Symbol"].tolist()
    mat = pd.read_table(mat_path)
    
    tfs = set(tfs).intersection(mat.columns)

    cluster = SLURMCluster(cores=1, 
                           processes=1,
                           memory="16GB",
                           account="amorin",
                           walltime="01:00:00")
    
    cluster.scale(jobs=4)

    custom_client = Client(cluster)

    print(custom_client)

    # print(f"Dashboard: {cluster.dashboard_link}")
    # print(f"Cluster jobs: {cluster.jobs}")

    # print('hello')
    time.sleep(30)
    # print(tfs)

    print(custom_client)
    print(cluster.scheduler_info())


    # start = time.time()

    # network = grnboost2(expression_data=mat, 
    #                     tf_names=tfs,
    #                     seed=5,
    #                     client_or_address=custom_client)
    
    # end = time.time()
    
    # print(end - start)

    # network.head()
    
    custom_client.close()
    cluster.close()
