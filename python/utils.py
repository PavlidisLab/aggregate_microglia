def iter_grnboost2(mat, genes, tfs, n_iter, out_dir, client):
    
    for i in range(n_iter):
        
        network = grnboost2(expression_data=mat, 
                            tf_names=tfs,
                            gene_names=genes,
                            seed=i,
                            client_or_address=client)
        
        file_path = Path(out_dir, f"GRNBoost2_TFonly_iter{i}.tsv")
        network.to_csv(file_path, sep="\t", header=False, index=False)
        
    return None



def iter_grnboost2_over_list(dat_list, 
                             ids, 
                             gene_df, 
                             tf_df, 
                             n_iter, 
                             out_root, 
                             client):
    
    for id in ids:
        
        print(f"{id}", datetime.datetime.now())
        out_dir = Path(out_root, id)
        out_dir.mkdir(exist_ok=True)

        mat = dat_list[id]["Mat"]
        keep_ix = mat.getnnz(axis=1) >= 20  # keep genes measured in 20+ cells
        mat_sub = mat[keep_ix, ].T  # to cells by genes
        genes_sub = gene_df.iloc[keep_ix]["Symbol"].tolist()
        tfs_sub = list(set(genes_sub).intersection(tf_df["Symbol"]))

        iter_grnboost2(mat_sub, genes_sub, tfs_sub, n_iter, out_dir, client)

    return None