from utils.io import *
import seaborn as sns

# generate figure 2-a
def Gen_figure_2a(sptFile, hires=False, dataPath = None):
    adata = Load_spt_to_AnnData(sptFile, "SPCS_mat", hires=hires, dataPath=dataPath)
    adata.obs["spTRS"] = adata.obs["BayesSpace"]
    adata0 = Load_spt_to_AnnData(sptFile, "matrix", hires=hires, dataPath=dataPath)
    figa, axsa = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    sc.pl.spatial(adata0, color="Seurat", size=1, ax=axsa[0, 0], show=False)
    sc.pl.spatial(adata0, color="Giotto", size=1, ax=axsa[0, 1], show=False)
    sc.pl.spatial(adata0, color="Squidpy", size=1, ax=axsa[1, 0], show=False)
    sc.pl.spatial(adata, color="spTRS", size=1, ax=axsa[1, 1], show=False)
    if not os.path.exists("figures"):
        os.mkdir("figures")
    figa.savefig("figures/fig2-a.svg")


def Gen_figure_2b(sptFile):
    adata = Load_spt_to_AnnData(sptFile)
    figa, axsa = plt.subplots(1, 2, figsize=(6.5, 3), constrained_layout=True)
    sc.pl.spatial(adata, img_key="lowres", ax=axsa[0], show=False)
    axsa[0].set_title("H&E image for 151673")
    sc.pl.spatial(adata, img_key="lowres", color="ground_truth", size=1, ax=axsa[1], show=False)
    if not os.path.exists("figures"):
        os.mkdir("figures")
    figa.savefig("figures/fig2-b.svg")

def Gen_figure_2c(tempdir):
    slices = ['151507', '151508', '151509', '151510',
              '151669', '151670', '151671', '151672',
              '151673', '151674', '151675', '151676']
    h5datas = ['matrix', 'SpotClean_mat', 'SPCS_mat']
    benchmarks = {}
    for slice in slices:
        sptFile = tempdir + "/All_" +slice +".h5spt"
        if os.path.isfile(sptFile):
            bmc = {}
            for h5data in h5datas:
                bm = Load_spt_to_Benchmark(sptFile, h5data)
                bmc[h5data] = bm
            benchmarks[slice] = bmc

    SPCS_Bayes = []
    for slice in slices:
        mk = benchmarks[slice]['SPCS_mat'].loc["BayesSpace", "ARI"]
        SPCS_Bayes.append(mk)

    bmlist = ['ARI', 'NMI', 'F.score']
    BM = {}
    for bm in bmlist:
        none_Bayes = []
        for slice in slices:
            mk = benchmarks[slice]['matrix'].loc["BayesSpace", bm]
            none_Bayes.append(mk)

        none_SpaGCN = []
        for slice in slices:
            mk = benchmarks[slice]['matrix'].loc["SpaGCN", bm]
            none_SpaGCN.append(mk)

        SPCS_SpaGCN = []
        for slice in slices:
            mk = benchmarks[slice]['SPCS_mat'].loc["SpaGCN", bm]
            SPCS_SpaGCN.append(mk)

        none_Giotto = []
        for slice in slices:
            mk = benchmarks[slice]['matrix'].loc["Giotto", bm]
            none_Giotto.append(mk)

        none_Seurat = []
        for slice in slices:
            mk = benchmarks[slice]['matrix'].loc["Seurat", bm]
            none_Seurat.append(mk)

        none_Squidpy = []
        for slice in slices:
            mk = benchmarks[slice]['matrix'].loc["Squidpy", bm]
            none_Squidpy.append(mk)

        bms = pd.DataFrame({ "Squidpy":none_Squidpy, "Giotto":none_Giotto, "Seurat":none_Seurat,"spTRS":SPCS_SpaGCN})
        BM[bm] = bms
    figa, axsa = plt.subplots(1, 3, figsize=(9, 3), constrained_layout=True)
    sns.boxplot(data=BM['ARI'], width=0.4, ax=axsa[0])
    axsa[0].set_title("ARI for 12 slices")
    sns.boxplot(data=BM['NMI'], width=0.4, ax=axsa[1])
    axsa[1].set_title("NMI for 12 slices")
    sns.boxplot(data=BM['F.score'], width=0.4, ax=axsa[2])
    axsa[2].set_title("F-score for 12 slices")
    figa.savefig("figures/fig2-c.svg")

def Gen_figure_2d(sptFile):
    bm = Load_spt_to_Benchmark(sptFile, h5data = "matrix", mode='cluster')
    bm0 = Load_spt_to_Benchmark(sptFile, h5data = "SPCS_mat", mode='cluster')
    bmari = pd.DataFrame({"ARI": bm['ARI'], 'attr':['raw'] * 8, 'method':bm.index})
    bmari = bmari.loc[["BayesSpace", "SpaGCN", "Giotto", "Seurat", "leiden"], :]
    bmari0 = pd.DataFrame({"ARI": bm0['ARI'], 'attr': ['imputed'] * 7, 'method': bm0.index})
    bmari0 = bmari0.loc[["BayesSpace", "SpaGCN", "Giotto", "Seurat", "leiden"], :]
    bmaris = pd.concat([bmari, bmari0])

    bmnmi = pd.DataFrame({"NMI": bm['NMI'], 'attr': ['raw'] * 8, 'method': bm.index})
    bmnmi = bmnmi.loc[["BayesSpace", "SpaGCN", "Giotto", "Seurat", "leiden"], :]
    bmnmi0 = pd.DataFrame({"NMI": bm0['NMI'], 'attr': ['imputed'] * 7, 'method': bm0.index})
    bmnmi0 = bmnmi0.loc[["BayesSpace", "SpaGCN", "Giotto", "Seurat", "leiden"], :]
    bmnmis = pd.concat([bmnmi, bmnmi0])

    # bmf = pd.DataFrame({"F.score": bm['F.score'], 'attr': ['raw'] * 8, 'method': bm.index})
    # bmf = bmf.loc[["BayesSpace", "SpaGCN", "Giotto", "Seurat", "leiden"], :]
    # bmf0 = pd.DataFrame({"F.score": bm0['F.score'], 'attr': ['imputed'] * 7, 'method': bm0.index})
    # bmf0 = bmf0.loc[["BayesSpace", "SpaGCN", "Giotto", "Seurat", "leiden"], :]
    # bmfs = pd.concat([bmf, bmf0])
    figd, axsd = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    sns.barplot(data=bmaris, x='method', y='ARI', hue='attr', ax=axsd[0])
    axsd[0].set_title("ARI")
    sns.barplot(data=bmnmis, x='method', y='NMI', hue='attr', ax=axsd[1])
    axsd[1].set_title("NMI")
    if not os.path.exists("figures"):
        os.mkdir("figures")
    figd.savefig("figures/fig2-d.svg")
    figd.savefig("figures/fig2-d.pdf")

def Gen_figure_2e(sptFile):
    h5data = "SPCS_mat"
    bme = Load_spt_to_Benchmark(sptFile, h5data, mode="estimate")
    BM = pd.DataFrame()
    for met in bme.keys():
        df = pd.DataFrame(bme[met]['moransI'])
        df.columns = [met]
        BM = pd.concat([BM, df], axis=1)
    showbm = BM[['SpaGCN',"SpatialDE", "Giotto"]]
    showbm.columns = ['spTRS', "SpatialDE", "Giotto"]
    fige, axse = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
    sns.boxplot(data=showbm, width=0.4, ax=axse)
    axse.set_title("Morans' I for 151673")
    fige.savefig("figures/fig2-e.svg")
    fige.savefig("figures/fig2-e.pdf")

def Gen_figure_2f(sptFile):
    h5data = "SPCS_mat"
    adata = Load_spt_to_AnnData(sptFile, h5data)
    HVG = adata.uns["HVGs"]['SpaGCN']
    HVG.sort_values('fold_change', ascending=False)
    group = HVG.groupby('cluster').head(2)
    hvgs = np.array(list(set(group['name'])), dtype='str')
    figa, axsa = plt.subplots(2, 3, figsize=(9, 6), constrained_layout=True)
    sc.pl.spatial(adata, img_key="lowres", color="COX6C", size=1, ax=axsa[0, 0], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="CAMK2N1", size=1, ax=axsa[0, 1], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="NEFL", size=1, ax=axsa[0, 2], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="DIRAS2", size=1, ax=axsa[1, 0], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="HPCAL1", size=1, ax=axsa[1, 1], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="MBP", size=1, ax=axsa[1, 2], show=False)
    if not os.path.exists("figures"):
        os.mkdir("figures")
    figa.savefig("figures/fig2-f.svg")
    figa.savefig("figures/fig2-f.pdf")


def Gen_figure_3a(sptFile):
    adata = Load_spt_to_AnnData(sptFile)
    figa, axsa = plt.subplots(1, 2, figsize=(6.5, 3), constrained_layout=True)
    sc.pl.spatial(adata, img_key="lowres", ax=axsa[0], show=False)
    axsa[0].set_title("H&E image for FFPE Kidney")
    adata.obs["spTRS"] = adata.obs["BayesSpace"]
    sc.pl.spatial(adata, img_key="lowres", color="spTRS", size=1, ax=axsa[1], show=False)
    if not os.path.exists("figures"):
        os.mkdir("figures")
    figa.savefig("figures/fig3-a.svg")
    figa.savefig("figures/fig3-a.pdf")

def Gen_figure_3bc(sptFile, Wfile):
    adata = Load_spt_to_AnnData(sptFile)
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    W.columns = ['CD45 B cell', 'CD45 T cell', 'CD45 macrophage',
                 'kidney distal convoluted tubule epithelial cell',
                 'brush cell', 'kidney collecting duct principal cell',
                 'kidney proximal convoluted tubule epithelial cell',
                 'podocyte', 'proximal tube epithelial cell',
                 'thick ascending tube S epithelial cell',
                 'Kidney cortex artery cell',
                 'fenestrated capillary endothelial cell',
                 'kidney capillary endothelial cell',
                 'Stroma fibroblast']
    obs = adata.obs
    adata.obs = pd.concat([obs, W], axis=1)
    W1 = W.loc[:,['kidney distal convoluted tubule epithelial cell',
                  'brush cell',
                  'fenestrated capillary endothelial cell']]
    figb1, axsb1 = plt.subplots(1, 3, figsize=(12, 3), constrained_layout=True)
    for i in range(3):
        sc.pl.spatial(adata, img_key="lowres", color=W1.columns[i], size=1, ax=axsb1[i], show=False)
    figb1.savefig("figures/fig3-b1.svg")
    figb1.savefig("figures/fig3-b1.pdf")

    W2 = W.loc[:, ['kidney collecting duct principal cell',
                   'kidney capillary endothelial cell',
                   'podocyte']]
    figb2, axsb2 = plt.subplots(1, 3, figsize=(12, 3), constrained_layout=True)
    for i in range(3):
        sc.pl.spatial(adata, img_key="lowres", color=W2.columns[i], size=1, ax=axsb2[i], show=False)
    figb2.savefig("figures/fig3-b2.svg")
    figb2.savefig("figures/fig3-b2.pdf")

    W3 = W.loc[:, ['kidney proximal convoluted tubule epithelial cell',
                   'proximal tube epithelial cell',
                   'thick ascending tube S epithelial cell']]
    figc, axsc = plt.subplots(1, 3, figsize=(12, 3), constrained_layout=True)
    for i in range(3):
        sc.pl.spatial(adata, img_key="lowres", color=W3.columns[i], size=1, ax=axsc[i], show=False)
    figc.savefig("figures/fig3-c.svg")
    figc.savefig("figures/fig3-c.pdf")

def Gen_figure_3d(sptFile, Wfile):
    adata = Load_spt_to_AnnData(sptFile)
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    W.columns = ['CD45 B cell', 'CD45 T cell', 'CD45 macrophage',
                 'kidney distal convoluted tubule epithelial cell',
                 'brush cell', 'kidney collecting duct principal cell',
                 'kidney proximal convoluted tubule epithelial cell',
                 'podocyte', 'proximal tube epithelial cell',
                 'thick ascending tube S epithelial cell',
                 'Kidney cortex artery cell',
                 'fenestrated capillary endothelial cell',
                 'kidney capillary endothelial cell',
                 'Stroma fibroblast']
    obs = adata.obs
    adata.obs = pd.concat([obs, W], axis=1)
    W0 = W.copy()
    W0['annotation'] = adata.obs['BayesSpace'].copy()
    df = W0.groupby(['annotation']).sum()
    from sklearn.preprocessing import MinMaxScaler

    scaler = MinMaxScaler()
    tran = scaler.fit_transform(df)
    df0 = pd.DataFrame(tran, columns=df.columns, index=df.index)
    figd, axsd = plt.subplots(1, 1, figsize=(7, 3), constrained_layout=True)
    sns.heatmap(df0.T, ax=axsd, cmap='Purples',linewidths=0.2)
    figd.savefig("figures/fig3-d.svg")
    figd.savefig("figures/fig3-d.pdf")

def Gen_figure_3e(sptFile, Wfile):
    adata = Load_spt_to_AnnData(sptFile)
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    W.columns = ['CD45 B cell', 'CD45 T cell', 'CD45 macrophage',
                 'kidney distal convoluted tubule epithelial cell',
                 'brush cell', 'kidney collecting duct principal cell',
                 'kidney proximal convoluted tubule epithelial cell',
                 'podocyte', 'proximal tube epithelial cell',
                 'thick ascending tube S epithelial cell',
                 'Kidney cortex artery cell',
                 'fenestrated capillary endothelial cell',
                 'kidney capillary endothelial cell',
                 'Stroma fibroblast']
    obs = adata.obs
    W0 = W.copy()
    W0['annotation'] = adata.obs['BayesSpace'].copy()
    df1 = pd.DataFrame(W0.groupby(['annotation']).mean())
    fige, axse = plt.subplots(1, 1, figsize=(7, 3.5), constrained_layout=True)
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
    fige.savefig("figures/fig3-e.svg")
    fige.savefig("figures/fig3-e.pdf")

def Gen_figure_3fg(sptFile):
    cntFile = "h5ads/references/Kidney/cnt_data.tsv"
    cnt = pd.read_csv(cntFile, sep='\t', index_col=0)
    mtaFile = "h5ads/references/Kidney/mta_data.tsv"
    mta = pd.read_csv(mtaFile, sep='\t', index_col=0)
    scdata = Load_spt_sc_to_AnnData(sptFile)
    from utils.stereoScope import Easy_Sample
    from utils.analysis import Spatial_Cluster_Analysis
    from sklearn import metrics
    scdata0 = Easy_Sample(scdata, 25, None)
    scdata1 = ad.AnnData(X=cnt)
    scdata0 = Spatial_Cluster_Analysis(scdata0)
    scdata1 = Spatial_Cluster_Analysis(scdata1)
    scdata = Spatial_Cluster_Analysis(scdata)
    ARI =  metrics.adjusted_rand_score(scdata.obs['clusters'], scdata.obs['annotation'])
    ARI0 = metrics.adjusted_rand_score(scdata0.obs['clusters'], scdata0.obs['annotation'])
    ARI1 = metrics.adjusted_rand_score(scdata1.obs['clusters'], mta['bio_celltype'])
    figf, axsf = plt.subplots(1, 1, figsize=(8.7, 4), constrained_layout=True)
    sc.pl.umap(scdata, color='annotation', palette=sc.pl.palettes.default_20, size=12, ax=axsf, show=False)
    axsf.set_title("Kidney scRNA-seq Data ARI={:.2f}".format(ARI))
    figf.savefig("figures/fig3-f.svg")
    figf.savefig("figures/fig3-f.pdf")

    obs = scdata1.obs
    mta.index = obs.index
    scdata1.obs = pd.concat([obs, mta], axis=1)
    figg, axsg = plt.subplots(1, 2, figsize=(7.5, 3.5), constrained_layout=True)
    sc.pl.umap(scdata0, color='clusters', palette=sc.pl.palettes.default_20, size=14, ax=axsg[0], show=False)
    axsg[0].set_title("Easy DownSample ARI={:.2f}".format(ARI0))
    sc.pl.umap(scdata1, color='clusters', palette=sc.pl.palettes.default_20, size=14, ax=axsg[1], show=False)
    axsg[1].set_title("VAE DownSample ARI={:.2f}".format(ARI1))
    figg.savefig("figures/fig3-g.svg")
    figg.savefig("figures/fig3-g.pdf")

def Generate_figure_4a(sptFile):
    adata = Load_spt_to_AnnData(sptFile, hires=True, dataPath=dataPath)

def Generate_figure_4d(sptFile):
    adata = Load_spt_to_AnnData(sptFile, hires=True, dataPath=dataPath)
    Wfile = "h5ads/res/spt_data/W.2022-04-14172940.888071.tsv"    # d14
    # Wfile = "h5ads/res/spt_data/W.2022-04-14171350.064126.tsv"    # d0
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    obs = adata.obs
    adata.obs = pd.concat([obs, W], axis=1)
    figb, axsb = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    for i in range(2):
        for j in range(3):
            if 3 * i + j >= 5:
                break
            sc.pl.spatial(adata, img_key="hires", color=W.columns[3 * i + j], size=1, ax=axsb[i, j], show=False)

def Generate_figure_5abcd():
    sptFiles = ["h5ads/FFPE_Human_Normal_Prostate.h5spt",
                "h5ads/FFPE_Human_Prostate_Acinar_Cell_Carcinoma.h5spt",
                "h5ads/FFPE_Human_Prostate_Cancer.h5spt",
                "h5ads/FFPE_Human_Prostate_IF.h5spt"]
    Adatas = []
    for sptFile in sptFiles:
        adata = Load_spt_to_AnnData(sptFile, count="SPCS_mat")
        Adatas.append(adata)
    figa, axsa = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    sc.pl.spatial(Adatas[0], color="BayesSpace", size=1, ax=axsa[0, 0], show=False)
    axsa[0, 0].set_title("Normal Prostate")
    sc.pl.spatial(Adatas[1], color="BayesSpace", size=1, ax=axsa[0, 1], show=False)
    axsa[0, 1].set_title("Acinar Cell Carcinoma")
    sc.pl.spatial(Adatas[2], color="BayesSpace", size=1, ax=axsa[1, 0], show=False)
    axsa[1, 0].set_title("Adeno & Invasive Carcinoma")
    sc.pl.spatial(Adatas[3], color="BayesSpace", size=1, ax=axsa[1, 1], show=False)
    axsa[1, 1].set_title("Adjacent Normal Section")
    figa.savefig("figures/fig5-a.svg")
    figa.savefig("figures/fig5-a.pdf")

    figb, axsb = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    sc.pl.spatial(Adatas[0], color="CDKN1A", size=1, ax=axsb[0, 0], show=False)
    axsb[0, 0].set_title("Normal Prostate")
    sc.pl.spatial(Adatas[1], color="CDKN1A", size=1, ax=axsb[0, 1], show=False)
    axsb[0, 1].set_title("Acinar Cell Carcinoma")
    sc.pl.spatial(Adatas[2], color="CDKN1A", size=1, ax=axsb[1, 0], show=False)
    axsb[1, 0].set_title("Adeno & Invasive Carcinoma")
    sc.pl.spatial(Adatas[3], color="CDKN1A", size=1, ax=axsb[1, 1], show=False)
    axsb[1, 1].set_title("Adjacent Normal Section")
    figb.suptitle("CDKN1A")
    figb.savefig("figures/fig5-b.svg")
    figb.savefig("figures/fig5-b.pdf")

    figb, axsb = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    sc.pl.spatial(Adatas[0], color="GSTP1", size=1, ax=axsb[0, 0], show=False)
    axsb[0, 0].set_title("Normal Prostate")
    sc.pl.spatial(Adatas[1], color="GSTP1", size=1, ax=axsb[0, 1], show=False)
    axsb[0, 1].set_title("Acinar Cell Carcinoma")
    sc.pl.spatial(Adatas[2], color="GSTP1", size=1, ax=axsb[1, 0], show=False)
    axsb[1, 0].set_title("Adeno & Invasive Carcinoma")
    sc.pl.spatial(Adatas[3], color="GSTP1", size=1, ax=axsb[1, 1], show=False)
    axsb[1, 1].set_title("Adjacent Normal Section")
    figb.suptitle("GSTP1")
    figb.savefig("figures/fig5-c.svg")
    figb.savefig("figures/fig5-c.pdf")

    figb, axsb = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    sc.pl.spatial(Adatas[0], color="APOD", size=1, ax=axsb[0, 0], show=False)
    axsb[0, 0].set_title("Normal Prostate")
    sc.pl.spatial(Adatas[1], color="APOD", size=1, ax=axsb[0, 1], show=False)
    axsb[0, 1].set_title("Acinar Cell Carcinoma")
    sc.pl.spatial(Adatas[2], color="APOD", size=1, ax=axsb[1, 0], show=False)
    axsb[1, 0].set_title("Adeno & Invasive Carcinoma")
    sc.pl.spatial(Adatas[3], color="APOD", size=1, ax=axsb[1, 1], show=False)
    axsb[1, 1].set_title("Adjacent Normal Section")
    figb.suptitle("APOD")
    figb.savefig("figures/fig5-d.svg")
    figb.savefig("figures/fig5-d.pdf")
    HVG0 = Adatas[0].uns["HVGs"]['SpaGCN']
    HVG0.sort_values('fold_change', ascending=False)
    group0 = HVG0.groupby('cluster')
    HVG1 = Adatas[1].uns["HVGs"]['SpaGCN']
    HVG1.sort_values('fold_change', ascending=False)
    group1 = HVG1.groupby('cluster')
    HVG2 = Adatas[2].uns["HVGs"]['SpaGCN']
    HVG2.sort_values('fold_change', ascending=False)
    group2 = HVG2.groupby('cluster')

def Generate_figure_5e(sptFile, Wfile):
    adata = Load_spt_to_AnnData(sptFile)
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    W = W[["B.cells", "NK", "Macrophage"]]
    obs = adata.obs
    adata.obs = pd.concat([obs, W], axis=1)
    figb1, axsb1 = plt.subplots(1, 3, figsize=(12, 3), constrained_layout=True)
    for i in range(3):
        sc.pl.spatial(adata, img_key="lowres", color=W.columns[i], size=1, ax=axsb1[i], show=False)
    figb1.savefig("figures/fig5-e1.svg")
    figb1.savefig("figures/fig5-e1.pdf")


def Generate_figure_5fghi(sptFile, Wfile):
    adata = Load_spt_to_AnnData(sptFile)
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    obs = adata.obs
    adata.obs = pd.concat([obs, W], axis=1)
    figb1, axsb1 = plt.subplots(3, 3, figsize=(12, 9), constrained_layout=True)
    for i in range(3):
        for j in range(3):
            sc.pl.spatial(adata, img_key="lowres", color=W.columns[i * 3 + j], size=1, ax=axsb1[i, j], show=False)
    figb1.savefig("figures/fig5-h.svg")
    figb1.savefig("figures/fig5-h.pdf")
    figb2, axsb2 = plt.subplots(3, 3, figsize=(12, 9), constrained_layout=True)
    for i in range(3):
        for j in range(3):
            sc.pl.spatial(adata, img_key="lowres", color=W.columns[i * 3 + j + 9], size=1, ax=axsb2[i, j], show=False)
    figb2.savefig("figures/fig5-i.svg")
    figb2.savefig("figures/fig5-i.pdf")

def Gen_figure_5j(sptFile, Wfile):
    adata = Load_spt_to_AnnData(sptFile, count="SPCS_mat")
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    obs = adata.obs
    W0 = W.copy()
    W0['annotation'] = adata.obs['SpaGCN'].copy()
    df1 = pd.DataFrame(W0.groupby(['annotation']).mean())
    fige, axse = plt.subplots(1, 1, figsize=(7, 4), constrained_layout=True)
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
    fige.savefig("figures/fig6-c.svg")
    fige.savefig("figures/fig6-c.pdf")

def Gen_figure_5k(sptFile, Wfile):
    adata = Load_spt_to_AnnData(sptFile)
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    obs = adata.obs
    adata.obs = pd.concat([obs, W], axis=1)
    W0 = W.copy()
    W0['annotation'] = adata.obs['BayesSpace'].copy()
    df = W0.groupby(['annotation']).sum()
    from sklearn.preprocessing import MinMaxScaler

    scaler = MinMaxScaler()
    tran = scaler.fit_transform(df)
    df0 = pd.DataFrame(tran, columns=df.columns, index=df.index)
    figd, axsd = plt.subplots(1, 1, figsize=(7, 3), constrained_layout=True)
    sns.heatmap(df0.T, ax=axsd, cmap='Purples',linewidths=0.2)
    figd.savefig("figures/fig6-d.svg")
    figd.savefig("figures/fig6-d.pdf")


def Gen_figure_5lm(sptFile):
    cntFile = "h5ads/references/FFPE_Human_Normal_Prostate/cnt_data.tsv"
    cnt = pd.read_csv(cntFile, sep='\t', index_col=0)
    mtaFile = "h5ads/references/FFPE_Human_Normal_Prostate/mta_data.tsv"
    mta = pd.read_csv(mtaFile, sep='\t', index_col=0)
    scdata = Load_spt_sc_to_AnnData(sptFile)
    from utils.stereoScope import Easy_Sample
    from utils.analysis import Spatial_Cluster_Analysis
    from sklearn import metrics
    scdata0 = Easy_Sample(scdata, 25, None)
    scdata1 = ad.AnnData(X=cnt)
    scdata0 = Spatial_Cluster_Analysis(scdata0)
    scdata1 = Spatial_Cluster_Analysis(scdata1)
    scdata = Spatial_Cluster_Analysis(scdata)
    ARI =  metrics.adjusted_rand_score(scdata.obs['clusters'], scdata.obs['annotation'])
    ARI0 = metrics.adjusted_rand_score(scdata0.obs['clusters'], scdata0.obs['annotation'])
    ARI1 = metrics.adjusted_rand_score(scdata1.obs['clusters'], mta['bio_celltype'])
    figf, axsf = plt.subplots(1, 1, figsize=(8.7, 4), constrained_layout=True)
    sc.pl.umap(scdata, color='annotation', palette=sc.pl.palettes.default_20, size=12, ax=axsf, show=False)
    axsf.set_title("Breast Cancer scRNA-seq Data ARI={:.2f}".format(ARI))
    figf.savefig("figures/fig6-e.svg")
    figf.savefig("figures/fig6-e.pdf")

    obs = scdata1.obs
    mta.index = obs.index
    scdata1.obs = pd.concat([obs, mta], axis=1)
    figg, axsg = plt.subplots(1, 2, figsize=(7.5, 2.8), constrained_layout=True)
    sc.pl.umap(scdata0, color='clusters', palette=sc.pl.palettes.default_20, size=14, ax=axsg[0], show=False)
    axsg[0].set_title("Easy DownSample ARI={:.2f}".format(ARI0))
    sc.pl.umap(scdata1, color='clusters', palette=sc.pl.palettes.default_20, size=14, ax=axsg[1], show=False)
    axsg[1].set_title("VAE DownSample ARI={:.2f}".format(ARI1))
    figg.savefig("figures/fig5-m.svg")
    figg.savefig("figures/fig5-m.pdf")

def Gen_figure_6a(sptFile, hires = False):
    adata = Load_spt_to_AnnData(sptFile, "SPCS_mat", hires=hires, dataPath=dataPath)
    adata.obs["spTRS"] = adata.obs["SpaGCN"]
    adata0 = Load_spt_to_AnnData(sptFile, "matrix", hires=hires, dataPath=dataPath)
    figa, axsa = plt.subplots(2, 2, figsize=(6, 6), constrained_layout=True)
    sc.pl.spatial(adata0, color="Seurat", size=1, ax=axsa[0, 0], show=False)
    sc.pl.spatial(adata0, color="Giotto", size=1, ax=axsa[0, 1], show=False)
    sc.pl.spatial(adata, color="spTRS", size=1, ax=axsa[1, 0], show=False)
    anno = mpimg.imread('data/Visium_FFPE_Human_Breast_Cancer/spatial/Visium_FFPE_Human_Breast_Cancer_Pathologist_Annotations.png')
    axsa[1, 1].imshow(anno)
    if not os.path.exists("figures"):
        os.mkdir("figures")
    figa.savefig("figures/fig6-a.pdf")
    figa.savefig("figures/fig6-a.svg")

def Generate_figure_6b(sptFile, Wfile):
    adata = Load_spt_to_AnnData(sptFile)
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    obs = adata.obs
    adata.obs = pd.concat([obs, W], axis=1)
    figb1, axsb1 = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    for i in range(2):
        for j in range(3):
            sc.pl.spatial(adata, img_key="lowres", color=W.columns[i * 3 + j], size=1, ax=axsb1[i, j], show=False)
    figb1.savefig("figures/fig6-b.pdf")

cfg = Load_json_configs("params.json")
tempdir = cfg['tempdir']
sptFile = cfg['sptFile']
dataPath = cfg['dataPath']
#Wfile = 'h5ads/res/spt_data/W.2022-04-27222940.309296.tsv' # prostate cancer
Wfile = 'h5ads/res/spt_data/W.2022-05-02112136.541108.tsv'  # breast cancer
Gen_figure_5j(sptFile, Wfile)
# # Gen_figure_2a(sptFile)
# # Gen_figure_2b(sptFile)
# # Gen_figure_2c(tempdir)
# # Gen_figure_2f(sptFile)
# Gen_figure_3a(sptFile)
# Gen_figure_3bc(sptFile, Wfile)
# Gen_figure_3d(sptFile, Wfile)
# Gen_figure_3e(sptFile, Wfile)
# Gen_figure_3fg(sptFile)
# Generate_figure_5abcd()
