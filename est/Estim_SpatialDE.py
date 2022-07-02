from utils.common import *

# Using SpatialDE to identify spatial variable genes
def HDG_Analysis_by_SpatialDE(adata,genes=200):
    import SpatialDE
    print(adata)
    counts = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs_names)
    coord = pd.DataFrame(adata.obsm["spatial"], columns=["x_coord", "y_coord"], index=adata.obs_names)
    results = SpatialDE.run(coord, counts)
    results.index = results["g"]
    adata.var = pd.concat([adata.var, results.loc[adata.var.index.values, :]], axis=1)
    return results.sort_values("qval").head(genes)