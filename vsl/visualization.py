from utils.io import *
from utils.preprocess import *
import seaborn as sns

# report the sptFile, including brief information, umap graphs, deconvolution composition
# PVL mapping right
# B-cell mapping wrong
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


def Show_Spatial(sptFile, figname, key, h5data='matrix', do_cluster=True):
    adata = Load_spt_to_AnnData(sptFile, h5data)
    fig, axs = plt.subplots(1, 1, figsize=(4, 3), constrained_layout=True)
    if do_cluster:
        adata = Remove_mito_genes(adata)
        adata = Filter_genes(adata)
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata = Spatial_Cluster_Analysis(adata)
    sc.pl.spatial(adata, ax=axs, color=key)
    fig.savefig(figname, dpi=400)


def Show_Violin(sptFile, figname, key, h5data='matrix'):
    adata = Load_spt_to_AnnData(sptFile, h5data)
    adata.var_names_make_unique()
    fig, axs = plt.subplots(1, 1, figsize=(4, 3), constrained_layout=True)
    adata = Remove_mito_genes(adata)
    adata = Filter_genes(adata)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata = Spatial_Cluster_Analysis(adata)
    sc.pl.violin(adata, ax=axs, groupby="BayesSpace", keys=key)
    fig.savefig(figname, dpi=400)


def Show_Umap(sptFile, figname, h5data='matrix', is_spatial=True, key='ground_truth', title=None, do_cluster=True):
    if is_spatial:
        adata = Load_spt_to_AnnData(sptFile, h5data)
    else:
        adata = Load_spt_sc_to_AnnData(sptFile, h5data)
    if do_cluster:
        adata = Remove_mito_genes(adata)
        adata = Filter_genes(adata)
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata = Spatial_Cluster_Analysis(adata)
    fig, axs = plt.subplots(1, 1, figsize=(9, 4), constrained_layout=True)
    sc.pl.umap(adata, ax=axs, color=key, show=False, palette=sc.pl.palettes.default_20, size=12)
    axs.set_title(title)
    fig.savefig(figname, dpi=400)


def Show_Spatial_and_Umap(sptFile, h5data='matrix',key='BayesSpace'):
    adata = Load_spt_to_AnnData(sptFile, h5data)
    adata = Remove_mito_genes(adata)
    adata = Filter_genes(adata)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata = Spatial_Cluster_Analysis(adata)
    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    sc.pl.spatial(adata, ax=axs[0], color=key, legend_loc=None)
    axs[0].set_title("Tissue Plot")
    sc.pl.umap(adata, ax=axs[1], color=key, show=False)
    axs[1].set_title("Umap Plot")
    fig.show()


def Show_Domain_and_Composition(sptFile, cell_type):
    info = sptInfo(sptFile)
    # arrangement of map chain: ["cell-type", "F1-score", "domain", "dct", "dcv", "clu"]
    map_chain = info.get_map_chains(cell_type)
    adata = Load_spt_to_AnnData(sptFile, map_chain['dct'])
    adata = Remove_mito_genes(adata)
    adata = Filter_genes(adata)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata = Spatial_Cluster_Analysis(adata)
    adata.obs['domain'] = list(adata.obs[map_chain['clu']])
    adata.obs['domain'][adata.obs['domain'] != map_chain['domain']] = -1
    adata.obs['domain'] = pd.Series(adata.obs['domain'], dtype='category')
    adata.obs[map_chain['dcv']] = adata.obsm[map_chain['dcv']][cell_type]
    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    sc.pl.spatial(adata, ax=axs[0], color='domain', legend_loc=None, cmap='coolwarm', show=False)
    axs[0].set_title("Domain Plot")
    sc.pl.spatial(adata, color=map_chain['dcv'], ax=axs[1], show=False)
    axs[1].set_title("Composition Plot")
    fig.show()


def Show_intersection_genes(sptFile, gene, folder, sptFile1=None, h5data='matrix', cell_type=""):
    adata_sp = Load_spt_to_AnnData(sptFile, h5data)
    fig, axs = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
    adata_sp = Remove_mito_genes(adata_sp)
    adata_sp = Filter_genes(adata_sp)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata_sp = Spatial_Cluster_Analysis(adata_sp)
    adata_sp.var_names = pd.Series(adata_sp.var_names).apply(str.upper)
    sc.pl.spatial(adata_sp, ax=axs, color=gene, show=False)
    fig.savefig(folder+'/'+gene+"_"+cell_type+".eps", dpi=400)

    fig, axs = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
    adata_sc1 = Load_spt_sc_to_AnnData(sptFile, "scRNA_seq")
    adata_sc1 = Remove_mito_genes(adata_sc1)
    adata_sc1 = Filter_genes(adata_sc1)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata_sc1 = Spatial_Cluster_Analysis(adata_sc1)
    adata_sc1.var_names = pd.Series(adata_sc1.var_names).apply(str.upper)
    sc.pl.umap(adata_sc1, ax=axs, color=gene, show=False)
    fig.savefig(folder+'/'+gene+"_"+cell_type+"_sc1.eps", dpi=400)

    if sptFile1 is not None:
        fig, axs = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
        adata_sc2 = Load_spt_sc_to_AnnData(sptFile1, "scRNA_seq")
        adata_sc2 = Remove_mito_genes(adata_sc2)
        adata_sc2 = Filter_genes(adata_sc2)
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata_sc2 = Spatial_Cluster_Analysis(adata_sc2)
        adata_sc2.var_names = pd.Series(adata_sc2.var_names).apply(str.upper)
        sc.pl.umap(adata_sc2, ax=axs, color=gene, size=2, show=False)
        fig.savefig(folder + '/' + gene +"_" +cell_type + "_sc2.eps", dpi=400)


def Show_Heatmap(sptFile, figname, h5data='matrix', title=None):
    adata = Load_spt_to_AnnData(sptFile, h5data)

    import seaborn as sns
    import scanpy as sc
    sc.pp.log1p(adata)
    X = pd.DataFrame(adata.X.todense(), index=adata.obs_names, columns=adata.var_names)
    X = X.iloc[0:200, 0:200]
    fig, axs = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
    sns.heatmap(X, cmap="viridis", ax=axs)
    axs.set_title(title)
    fig.savefig(figname, dpi=400)

# downsampling comparison. this function draws a figure to compare cell expression patterns
# between downsampling and VAE.
def Show_DS_comparison(sptFile, figname, h5data='scRNA_seq'):
    from utils.stereoScope import Easy_Sample
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    from sklearn import metrics
    scdata = Load_spt_sc_to_AnnData(sptFile, 'scRNA_seq')
    scdata0 = Easy_Sample(scdata, 25, None)
    scdata1 = Load_spt_sc_to_AnnData(sptFile, h5data)
    scdata0 = Spatial_Cluster_Analysis(scdata0)
    scdata1 = Spatial_Cluster_Analysis(scdata1)
    scdata = Spatial_Cluster_Analysis(scdata)
    ARI = metrics.adjusted_rand_score(scdata.obs['clusters'], scdata.obs['annotation'])
    ARI0 = metrics.adjusted_rand_score(scdata0.obs['clusters'], scdata0.obs['annotation'])
    ARI1 = metrics.adjusted_rand_score(scdata1.obs['clusters'], scdata1.obs['annotation'])
    fig, axs = plt.subplots(2, 1, figsize=(3, 4), constrained_layout=True)
    sc.pl.umap(scdata0, color='clusters', palette=sc.pl.palettes.default_20, size=14, ax=axs[0], show=False)
    axs[0].set_title("Easy DownSample ARI={:.2f}".format(ARI0))
    sc.pl.umap(scdata1, color='clusters', palette=sc.pl.palettes.default_20, size=14, ax=axs[1], show=False)
    axs[1].set_title("VAE DownSample ARI={:.2f}".format(ARI1))
    fig.savefig(figname, dpi=400)


