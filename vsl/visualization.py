from utils.io import *
from utils.preprocess import *

# report the sptFile, including brief information, umap graphs, deconvolution composition
def spTRS_report(sptFile, count="matrix"):
    # report brief information
    adata = Load_spt_to_AnnData(sptFile, count)
    figd, axsd = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
    sc.pl.spatial(adata, color="ground_truth", show=False, ax = axsd)
    figd.savefig("spatial1.pdf")
    # report deconvolution heatmap
    W = adata.obsm['spacexr']
    W0 = W.copy()
    W0['annotation'] = adata.obs['ground_truth'].copy()
    df = W0.groupby(['annotation']).sum()
    from sklearn.preprocessing import MinMaxScaler

    scaler = MinMaxScaler()
    tran = scaler.fit_transform(df)
    df0 = pd.DataFrame(tran, columns=df.columns, index=df.index)
    figd, axsd = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
    import seaborn as sns
    sns.heatmap(df0.T, ax=axsd, cmap='Purples',linewidths=0.2)
    figd.savefig("spacexr1.pdf")

    # report deconvolution composition
    W0 = W.copy()
    W0['annotation'] = adata.obs['ground_truth'].copy()
    df1 = pd.DataFrame(W0.groupby(['annotation']).mean())
    fige, axse = plt.subplots(1, 1, figsize=(6, 3.5), constrained_layout=True)
    width = 0.5
    x = df1.index
    bottom_y = np.array([0] * len(df1))
    for cell,color in zip(df1.columns, mcolors.XKCD_COLORS):
        y = df1.loc[:, cell]
        axse.bar(x, y, width, bottom=bottom_y, label=W.columns, color=color)
        bottom_y = y + bottom_y
    axse.set_title('The cellular composition of each region')
    axse.set_xticks(x, x)
    axse.legend(W.columns, loc=2, bbox_to_anchor=(1.0, 1.0),borderaxespad=0.)


def Show_Origin_histology(sptFile, figname, hires=True, key=None):
    adata = Load_spt_to_AnnData(sptFile, hires=hires)
    fig, axs = plt.subplots(1, 1, figsize=(6, 3), constrained_layout=True)
    sc.pl.spatial(adata, img_key='hires', ax=axs, color=key)
    fig.savefig(figname, dpi=400)


def Show_Umap(sptFile, figname, is_spatial=True, key='ground_truth', title=None):
    if is_spatial:
        adata = Load_spt_to_AnnData(sptFile)
    else:
        adata = Load_spt_sc_to_AnnData(sptFile)
    adata = Remove_mito_genes(adata)
    adata = Filter_genes(adata)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata = Spatial_Cluster_Analysis(adata)
    fig, axs = plt.subplots(1, 1, figsize=(6, 3), constrained_layout=True)
    sc.pl.umap(adata, ax=axs, color=key, show=False)
    axs.set_title(title)
    fig.savefig(figname, dpi=400)
