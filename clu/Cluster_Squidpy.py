from utils.common import *


# helper function returning a clustering
def squidpy_cluster_features(features: pd.DataFrame, like=None) -> pd.Series:
    """
    Calculate leiden clustering of features.

    Specify filter of features using `like`.
    """
    # filter features
    if like is not None:
        features = features.filter(like=like)
    # create temporary adata to calculate the clustering
    adata = ad.AnnData(features)
    # important - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata)
    # calculate leiden clustering
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)

    return adata.obs["leiden"]


# Using squidpy to do dimension reduction and clustering
def Spatial_Cluster_Analysis_squidpy(adata, tif):
    import squidpy as sq
    img = sq.im.ImageContainer(tif)
    # calculate features for different scales (higher value means more context)
    for scale in [1.0, 2.0]:
        feature_name = f"features_summary_scale{scale}"
        sq.im.calculate_image_features(
            adata,
            img.compute(),
            features="summary",
            key_added=feature_name,
            n_jobs=4,
            scale=scale,
        )

    # combine features in one dataframe
    adata.obsm["features"] = pd.concat(
        [adata.obsm[f] for f in adata.obsm.keys() if "features_summary" in f], axis="columns"
    )
    # make sure that we have no duplicated feature names in the combined table
    adata.obsm["features"].columns = ad.utils.make_index_unique(adata.obsm["features"].columns)

    # calculate feature clusters
    adata.obs["features_cluster"] = squidpy_cluster_features(adata.obsm["features"], like="summary")
    return adata
