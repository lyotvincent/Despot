from utils.common import *

# 数据质控
def Data_Qc(adata):
    sc.pl.highest_expr_genes(adata, n_top=20)


# 移除线粒体基因
def Remove_mito_genes(adata):
    mito_genes = adata.var_names.str.startswith('MT-')
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['mt_frac'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['total_counts'] = adata.X.sum(axis=1).A1
    return adata

# 各项数据的信息
def Preprocess_info(adata):
    print("Data Abstracts:")
    print(adata)
    # plot for Covariates for filtering
    fig, axs = plt.subplots(1, 2, figsize=(15, 4))
    fig.suptitle('Covariates for filtering')
    sb.distplot(adata.obs['total_counts'], kde=False, ax=axs[0])
    sb.distplot(adata.obs['total_counts'][adata.obs['total_counts'] < 10000], kde=False, bins=40, ax=axs[1])

# 基因过滤
def Filter_genes(adata, min_max_counts=(10, 350000), mt_frac=0.1, min_genes_cells=(3000, 10)):
    mincounts, maxcounts = min_max_counts[0], min_max_counts[1]
    mingenes, mincells = min_genes_cells[0], min_genes_cells[1]
    sc.pp.filter_cells(adata, min_counts=mincounts)
    print(f'Number of cells after min count filter: {adata.n_obs}')
    # sc.pp.filter_cells(adata, max_counts=maxcounts)
    # print(f'Number of cells after max count filter: {adata.n_obs}')
    # adata = adata[adata.obs['mt_frac'] < mt_frac]
    # print(f'Number of cells after MT filter: {adata.n_obs}')
    # sc.pp.filter_cells(adata, min_genes=mincells)
    # print(f'Number of cells after gene filter: {adata.n_obs}')
    sc.pp.filter_genes(adata, min_cells=mincells)
    print(f'Number of genes after cell filter: {adata.n_vars}')
    return adata


# 用stLearn进行预处理
def Preprocess_stLearn(stdata):
    import stlearn as st
    # spot tile is the intermediate result of image pre-processing
    TILE_PATH = Path("/tmp/tiles")
    TILE_PATH.mkdir(parents=True, exist_ok=True)

    # pre-processing for gene count table
    st.pp.filter_genes(stdata, min_cells=1)
    st.pp.normalize_total(stdata)
    st.pp.log1p(stdata)

    # pre-processing for spot image
    st.pp.tiling(stdata, TILE_PATH)

    # this step uses deep learning model to extract high-level features from tile images
    # may need few minutes to be completed
    st.pp.extract_feature(stdata)

    return stdata