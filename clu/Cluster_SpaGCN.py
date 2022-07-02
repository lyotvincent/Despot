from utils.common import *


# Using SpaGCN to do dimension reduction and clustering
def Spatial_Cluster_Analysis_SpaGCN(adata):
    array_row = adata.obs['array_row'].tolist()
    array_col = adata.obs['array_col'].tolist()
    pixel_row = adata.obs['image_row'].tolist()
    pixel_col = adata.obs['image_col'].tolist()
    import SpaGCN as spg
    import torch
    adj = spg.calculate_adj_matrix(array_row, array_col, histology=False)
    adata.var_names_make_unique()
    # Spatial Domain Detection Using SpaGCN
    spg.prefilter_genes(adata)
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    p = 0.5
    l = spg.search_l(p, adj)
    clf = spg.SpaGCN()
    clf.set_l(l)
    random.seed(100)
    torch.manual_seed(100)
    np.random.seed(100)
    clf.train(adata, adj)
    pred, prob = clf.predict()
    adata.obs['clusters'] = pred
    adata.obs['clusters'] = adata.obs['clusters'].astype("category")
    adj_2d = spg.calculate_adj_matrix(array_row, array_col, histology=False)
    refine_pred = spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs['clusters'].tolist(), dis=adj_2d)
    adata.obs['clusters'] = refine_pred
    adata.obs['clusters'] = adata.obs['clusters'].astype("category")
    return adata
