from utils.common import *


# using SpaGCN to identify high variable genes
def HDG_Analysis_by_SpaGCN(adata,
                           min_in_group_fraction=0.8,
                           min_in_out_group_ratio=0.5,
                           min_fold_change=1.5):
    import SpaGCN as spg
    array_row = adata.obs['array_row'].tolist()
    array_col = adata.obs['array_col'].tolist()
    array_pred = adata.obs['BayesSpace'].tolist()
    targets = list(set(array_pred))
    sc.pp.log1p(adata)
    adj_2d = spg.calculate_adj_matrix(x=array_row, y=array_col, histology=False)
    start, end = np.quantile(adj_2d[adj_2d != 0], q=0.001), np.quantile(adj_2d[adj_2d != 0], q=0.1)
    filtered_info = pd.DataFrame()
    for target in targets:
        r = spg.search_radius(target_cluster=target, cell_id=adata.obs.index.tolist(), x=array_row, y=array_col,
                              pred=array_pred, start=start, end=end, num_min=10, num_max=14, max_run=100)
        # Detect neighboring domains
        nbr_domains = spg.find_neighbor_clusters(target_cluster=target,
                                                 cell_id=adata.obs.index.tolist(),
                                                 x=array_row,
                                                 y=array_col,
                                                 pred=array_pred,
                                                 radius=r,
                                                 ratio=1 / 10)
        if nbr_domains is None:
            continue
        else:
            nbr_domains = nbr_domains[0:3]
        de_genes_info = spg.rank_genes_groups(input_adata=adata,
                                              target_cluster=target,
                                              nbr_list=nbr_domains,
                                              label_col="BayesSpace",
                                              adj_nbr=True,
                                              log=True)
        # Filter genes
        de_genes_info = de_genes_info[(de_genes_info["pvals_adj"] < 0.05) &
                                   (de_genes_info["in_out_group_ratio"] > min_in_out_group_ratio) &
                                   (de_genes_info["in_group_fraction"] > min_in_group_fraction) &
                                   (de_genes_info["fold_change"] > min_fold_change)]
        de_genes_info = de_genes_info.sort_values(by="in_group_fraction", ascending=False)
        de_genes_info["target_dmain"] = target
        de_genes_info["neighbors"] = str(nbr_domains)
        filtered_info = filtered_info.append(de_genes_info, ignore_index=True)
    adata.uns['filtered_info'] = filtered_info
    return adata