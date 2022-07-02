from utils.common import *


# Using leiden to do dimension reduction and clustering
def Spatial_Cluster_Analysis(adata, n_top_genes=2000, n_comps=50):
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', inplace=True)
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata)
    # use umap and leiden for clustering
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added='clusters')
    return adata