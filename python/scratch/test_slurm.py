import pandas as pd
import time


if __name__ == '__main__':

    tfs_path = "/home/amorin/Data/Metadata/TFs_human.tsv"
    mat_path = "/space/scratch/amorin/R_objects/GSE180928_mcg_filt.tsv"

    tfs = pd.read_table(tfs_path)["Symbol"].tolist()
    mat = pd.read_table(mat_path)
    
    tfs = set(tfs).intersection(mat.columns)

    print(len(tfs))
