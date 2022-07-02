from utils.common import *


# Using stSME to do dimension reduction and clustering
def Spatial_Cluster_Analysis_stSME(stdata):
    import stlearn as st
    # run PCA for gene expression data
    st.em.run_pca(stdata, n_comps=50)

    data_SME = stdata.copy()
    # apply stSME to normalise log transformed data
    st.spatial.SME.SME_normalize(data_SME, use_data="raw")
    data_SME.X = data_SME.obsm['raw_SME_normalized']
    st.pp.scale(data_SME)
    st.em.run_pca(data_SME, n_comps=50)

    # louvain clustering on stSME normalised data
    st.pp.neighbors(data_SME, n_neighbors=17, use_rep='X_pca')
    st.tl.clustering.louvain(data_SME, resolution=1.19)
    return data_SME
