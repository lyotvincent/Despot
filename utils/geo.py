import esda
import pandas as pd
import geopandas as gpd
import scvi
from geopandas import GeoDataFrame
import libpysal as lps
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from utils.io import *
from utils.purity import purity_score
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score
import alphashape
import networkx as nx


# adata = Load_spt_to_AnnData("h5ads/FFPE_Mouse_Kidney.h5spt")
# geometry = [Point(xy) for xy in zip(adata.obs['array_row'], adata.obs['array_col'])]
# ad_gdf = GeoDataFrame(adata.obsm['spacexr'], geometry=geometry)
# wq = lps.weights.Queen.from_dataframe(ad_gdf)
# wq.transform = 'r'
# MI = list(map(lambda x: esda.moran.Moran(ad_gdf.iloc[:, x], wq), list(range(14))))
# I = list(map(lambda x: MI[x].I, list(range(14))))
# print(I)
#
# bcell = ad_gdf.iloc[:, 0]
# bcell_lag = lps.weights.lag_spatial(wq, bcell)
# b, a = np.polyfit(x=np.array(bcell), y=np.array(bcell_lag), deg=1)
# plt.plot(bcell, bcell_lag, '.', color='firebrick')
#
#  # dashed vert at mean of the price
# plt.vlines(bcell.mean(), bcell_lag.min(), bcell_lag.max(), linestyle='--')
#  # dashed horizontal at mean of lagged price
# plt.hlines(bcell_lag.mean(), bcell.min(), bcell.max(), linestyle='--')
#
# # red line of best fit using global I as slope
# plt.plot(bcell, a + b*bcell, 'r')
# plt.title('Moran Scatterplot')
# plt.ylabel('Spatial Lag of Price')
# plt.xlabel('Price')
# plt.show()
# # y = ad_gdf.iloc[:, 0]
# # ylag = lps.weights.lag_spatial(wq, y)
# # mi = esda.moran.Moran(y, wq)
# li = esda.moran.Moran_Local(bcell, wq)
# sig = 1 * (li.p_sim < 0.05)
# hotspot = 1 * (sig * li.q==1)
# coldspot = 3 * (sig * li.q==3)
# doughnut = 2 * (sig * li.q==2)
# diamond = 4 * (sig * li.q==4)
# spots = hotspot + coldspot + doughnut + diamond
# spot_labels = [ '0 ns', '1 hot spot', '2 doughnut', '3 cold spot', '4 diamond']
# p_moransI = [spot_labels[i] for i in spots]
# from matplotlib import colors
# hmap = colors.ListedColormap([ 'lightgrey', 'red', 'lightblue', 'blue', 'pink'])
# f, ax = plt.subplots(1, figsize=(9, 9))
# ad_gdf.assign(cl=p_moransI).plot(column='cl', categorical=True,
#         k=2, cmap=hmap, linewidth=0.1, ax=ax,
#         edgecolor='white', legend=True)
# ax.set_axis_off()
# plt.show()


def p_moransI_purity(adata, dcv_mtd: str, clu: str):
    geometry = [Point(xy) for xy in zip(adata.obs['array_row'], adata.obs['array_col'])]
    # first calculate global moransI
    ad_gdf = GeoDataFrame(adata.obsm[dcv_mtd], geometry=geometry, copy=True)
    wq = lps.weights.Queen.from_dataframe(ad_gdf)
    wq.transform = 'r'
    perf = pd.DataFrame(columns=['moransI', 'ARI', 'NMI', 'purity'])
    for i in range(len(adata.obsm[dcv_mtd].columns) - 1):
        comp = adata.obsm[dcv_mtd].columns[i]
        if comp == 'geometry':
            break
        y = ad_gdf.loc[:, comp]
        MI = esda.moran.Moran(y, wq)
        y_lag = lps.weights.lag_spatial(wq, y)
        # then calculate local moransI
        li = esda.moran.Moran_Local(y, wq)
        # using lisa cluster to find hotpoint, coldpoint, doughnut and diamond
        sig = 1 * (li.p_sim < 0.05)
        hotspot = 1 * (sig * li.q == 1)
        coldspot = 3 * (sig * li.q == 3)
        doughnut = 2 * (sig * li.q == 2)
        diamond = 4 * (sig * li.q == 4)
        spots = hotspot + coldspot + doughnut + diamond
        # spot_labels = ['0 ns', '1 hot spot', '2 doughnut', '3 cold spot', '4 diamond']
        # p_moransI = [spot_labels[i] for i in spots]
        # select the point in lisa cluster
        adata.obs['p_moransI_' + comp] = spots
        adata0 = adata[adata.obs['p_moransI_' + comp] > 0]
        # calculate the purity with the segmentation
        pred = adata0.obs[clu]
        true = adata0.obs['p_moransI_' + comp]
        ps = purity_score(true, pred)
        ari = adjusted_rand_score(true, pred)
        nmi = normalized_mutual_info_score(true, pred)
        pf = pd.DataFrame(data=[MI.I, ari, nmi, ps], index=['moransI', 'ARI', 'NMI', 'purity'], columns=[comp])
        pf = pf.T
        perf = pd.concat([perf, pf], ignore_index=False)
    return perf


def suitable_location(adata, perf, clu, beta=1):
    cell_type = perf.index
    SL = pd.DataFrame()
    for ct in cell_type:
        spots = adata.obs['p_moransI_' + ct]
        clusters = adata.obs[clu]
        clu_num = np.unique(adata.obs[clu])
        temp = pd.concat([spots, clusters], axis=1)
        HH = []
        for cl in clu_num:
            temp0 = temp[temp[clu] == cl]
            precision = len(temp0[temp0['p_moransI_' + ct] == 1]) / len(temp0)
            temp1 = temp[temp['p_moransI_' + ct] == 1]
            recall = len(temp1[temp1[clu] == cl]) / len(temp1)
            if precision == 0 and recall == 0:
                HH.append(0)
            else:
                Fscore = (1 + beta ** 2) * (precision * recall) / (beta ** 2 * (precision + recall))
                HH.append(Fscore)
        HH = pd.DataFrame(HH, index=clu_num)
        SL = pd.concat([SL, HH], axis=1)
    SL.columns = cell_type
    return SL


def spTRS_findgroups(sptFile, beta=1):
    sptinfo = sptInfo(sptFile)
    spmats = sptinfo.get_spmatrix()
    Groups = {}
    Best_Fscore = {}
    for spmat in spmats:
        adata = Load_spt_to_AnnData(sptFile, spmat)
        clu_mtds = sptinfo.get_clu_methods(spmat)
        dcv_mtds = sptinfo.get_dcv_methods(spmat)
        print(clu_mtds)
        print(dcv_mtds)
        for dcv in dcv_mtds:
            for clu in clu_mtds:
                mtd_name = "{0}+{1}+{2}".format(spmat, dcv, clu)
                print("Finding Groups in the combination of {0}".format(mtd_name))
                perf = p_moransI_purity(adata, dcv_mtd=dcv, clu=clu)
                perf = perf[perf['moransI'] > 0.5]
                SL = suitable_location(adata, perf, clu, beta)
                BST = pd.concat([SL.max(axis=0), SL.idxmax(axis=0)], axis=1)
                BST.columns = ['F-score', 'domain']
                Groups[mtd_name] = SL
                Best_Fscore[mtd_name] = BST
    return Groups, Best_Fscore


def spTRS_Find_bestGroup(sptFile, beta=1):
    Groups, Best_Fscore = spTRS_findgroups(sptFile, beta)
    Best_dict = {}
    for mtd in Best_Fscore.keys():
        mtd_Fscore = Best_Fscore[mtd]
        # for each celltype, find the cluster which has the highest Fscore
        for ct in mtd_Fscore.index:
            if ct not in Best_dict.keys():
                Best_dict[ct] = (mtd_Fscore.loc[ct, 'F-score'], mtd_Fscore.loc[ct, 'domain'], mtd)
            else:
                if mtd_Fscore.loc[ct, 'F-score'] > Best_dict[ct][0]:
                    Best_dict[ct] = (mtd_Fscore.loc[ct, 'F-score'], mtd_Fscore.loc[ct, 'domain'], mtd)
    # find markers for each region
    # for cell_type in Best_dict.keys():
    #     region = Best_dict[cell_type][1]
    #     mtd = Best_dict[cell_type][2]
    #     clu_mtd = mtd.split('+')[2]
    #     dcv_mtd = mtd.split('+')[1]
    #     dct_mtd = mtd.split('+')[0]
    #     adata = Load_spt_to_AnnData(sptFile, count=dct_mtd)
    #     sc.tl.filter_rank_genes_groups(adata, key=region, groupby=clu_mtd)
    return Best_dict, Groups


def Fscore_Comparison(Best_dict: dict, Groups):
    # find f-score about ground_truth for each pipeline
    ground_truth_max = pd.DataFrame()
    Seurat_max = pd.DataFrame()
    Giotto_max = pd.DataFrame()
    # spTRS compare with Seurat and Giotto
    for mtd in Groups.keys():
        clu_mtd = mtd.split('+')[2]
        if clu_mtd == 'ground_truth':
            group = Groups[mtd]
            ground_truth_max = ground_truth_max.append(pd.DataFrame(group.max()).T)
        dcv_mtd = mtd.split('+')[1]
        dct_mtd = mtd.split('+')[0]
        if dct_mtd == 'matrix' and clu_mtd == 'Seurat' and dcv_mtd == 'Seurat':
            group = Groups[mtd]
            Seurat_max = Seurat_max.append(pd.DataFrame(group.max()).T)
        if dct_mtd == 'matrix' and clu_mtd == 'Giotto' and dcv_mtd == 'Giotto':
            group = Groups[mtd]
            Giotto_max = Giotto_max.append(pd.DataFrame(group.max()).T)
    ground_truth_max = pd.DataFrame(ground_truth_max.max()).T
    ground_truth_max.index = ['ground_truth']
    Seurat_max = pd.DataFrame(Seurat_max.max()).T
    Seurat_max.index = ['Seurat']
    ground_truth_max = ground_truth_max.append(Seurat_max)
    Giotto_max = pd.DataFrame(Giotto_max.max()).T
    Giotto_max.index = ['Giotto']
    ground_truth_max = ground_truth_max.append(Giotto_max)
    cell_type = []
    fscore = []
    for key in Best_dict.keys():
        cell_type.append(key)
        fscore.append(Best_dict[key][0])
    spTRS_max = pd.DataFrame(data=fscore, index=cell_type, columns=['spTRS']).T
    ground_truth_max = ground_truth_max.append(spTRS_max)

    # find f-score about ground_truth for each pipeline
    undct_max = pd.DataFrame()
    dct_max = pd.DataFrame()
    for mtd in Groups.keys():
        dct_mtd = mtd.split('+')[0]
        if dct_mtd == 'matrix':
            group = Groups[mtd]
            undct_max = undct_max.append(pd.DataFrame(group.max()).T)
        else:
            group = Groups[mtd]
            dct_max = dct_max.append(pd.DataFrame(group.max()).T)
    undct_max = pd.DataFrame(undct_max.max()).T
    undct_max.index = ['matrix']
    ground_truth_max = ground_truth_max.append(undct_max)
    dct_max = pd.DataFrame(dct_max.max()).T
    dct_max.index = ['decontaminated']
    ground_truth_max = ground_truth_max.append(dct_max)
    # comparison among vae, es and na
    vae_max = pd.DataFrame()
    es_max = pd.DataFrame()
    na_max = pd.DataFrame()
    for mtd in Groups.keys():
        dcv_mtd = mtd.split('+')[1]
        if dcv_mtd in ["StereoScope", "spacexr", "SPOTlight_vae"]:
            group = Groups[mtd]
            vae_max = vae_max.append(pd.DataFrame(group.max()).T)
        elif dcv_mtd in ["StereoScope_es", "spacexr_es", "SPOTlight_es"]:
            group = Groups[mtd]
            es_max = es_max.append(pd.DataFrame(group.max()).T)
        elif dcv_mtd in ['StereoScope_na', 'SPOTlight']:
            group = Groups[mtd]
            na_max = na_max.append(pd.DataFrame(group.max()).T)
    vae_max = pd.DataFrame(vae_max.max()).T
    vae_max.index = ['vae']
    ground_truth_max = ground_truth_max.append(vae_max)

    es_max = pd.DataFrame(es_max.max()).T
    es_max.index = ['es']
    ground_truth_max = ground_truth_max.append(es_max)

    na_max = pd.DataFrame(na_max.max()).T
    na_max.index = ['na']
    pip_res = ground_truth_max.append(na_max)
    return pip_res


def Show_best_group(sptFile, cell_type, fscore, domain, mtd_chain: str, figname):
    mtds = mtd_chain.split('+')
    dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
    abd_name = dct_mtd + '_' + dcv_mtd + '_' + clu_mtd + '_' + cell_type
    adata = Load_spt_to_AnnData(sptFile, count=dct_mtd)
    adata.obs[abd_name] = adata.obsm[dcv_mtd][cell_type]
    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    sc.pl.spatial(adata, color=clu_mtd, ax=axs[0], show=False)
    sc.pl.spatial(adata, color=abd_name, ax=axs[1], show=False)
    axs[1].set_title(dcv_mtd)
    fig.suptitle("find {0} locate in {1} using methods chain:"
                 "\n{2}, F-score={3}".format(cell_type, str(domain), abd_name, fscore))
    fig.savefig(figname, dpi=400)


def Show_bash_best_group(sptFile, Best_dict, folder):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    cell_types = list(Best_dict.keys())
    for cell_type in cell_types:
        Best_item = Best_dict[cell_type]
        Show_best_group(sptFile, cell_type, Best_item[0], Best_item[1], Best_item[2], folder + "/" + cell_type + ".eps")


def spTRS_self_correlation(sptFile, alpha=1):
    # using scipy to calculate pearson's r
    from scipy.stats import pearsonr, spearmanr
    sptinfo = sptInfo(sptFile)
    h5datas = sptinfo.get_spmatrix()
    for h5data in h5datas:
        adata = Load_spt_to_AnnData(sptFile, h5data)
        dcv_mtds = sptinfo.get_dcv_methods(h5data)
        for dcv in dcv_mtds:
            geometry = [Point(xy) for xy in zip(adata.obs['array_row'], adata.obs['array_col'])]
            # first calculate global moransI
            ad_gdf = GeoDataFrame(adata.obsm[dcv], geometry=geometry, copy=True)
            wq = lps.weights.Queen.from_dataframe(ad_gdf)
            wq.transform = 'r'
            # we combine the cell enrichment with the neighbor plot
            comp = pd.DataFrame(adata.obsm[dcv])
            y_lag = pd.DataFrame(list(map(lambda i: lps.weights.lag_spatial(wq, comp.loc[:, i]), comp.columns))).T
            y_lag.columns, y_lag.index = comp.columns, comp.index
            # set alpha=0 to cancel the efficiency of neighbor expressions
            comp = comp + alpha * y_lag
            # normalize the comp
            ct_mean = np.mean(comp, axis=0)
            ct_std = np.std(comp, axis=0)
            comp = (comp - ct_mean) / ct_std
            # comp = comp + lags
            coef = np.corrcoef(comp.T)
            corr = np.zeros((comp.shape[1], comp.shape[1]))
            p_val = np.zeros((comp.shape[1], comp.shape[1]))
            for x in range(p_val.shape[0]):
                for y in range(p_val.shape[1]):
                    corr[x, y], p_val[x, y] = spearmanr(comp.iloc[:, x], comp.iloc[:, y])
            coef = pd.DataFrame(coef, index=comp.columns, columns=comp.columns)
            Save_spt_from_corrcoef(sptFile, coef, dcv, h5data = h5data)


def spTRS_combine_best_group(sptFile, Best_dict):
    cell_types = list(Best_dict.keys())
    abundance = pd.DataFrame()
    sptinfo = sptInfo(sptFile)
    dct_mtds = sptinfo.get_spmatrix()
    adatas = list(map(lambda dct: Load_spt_to_AnnData(sptFile, count=dct), dct_mtds))
    adatas = dict(zip(dct_mtds, adatas))
    for cell_type in cell_types:
        Best_item = Best_dict[cell_type]
        mtds = Best_item[2].split('+')
        dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
        abundance = pd.concat([abundance, adatas[dct_mtd].obsm[dcv_mtd][cell_type]], axis=1)
    geometry = [Point(xy) for xy in zip(adatas['matrix'].obs['array_row'], adatas['matrix'].obs['array_col'])]
    ad_gdf = GeoDataFrame(abundance, geometry=geometry, copy=True)
    wq = lps.weights.Queen.from_dataframe(ad_gdf)
    wq.transform = 'r'
    comp = abundance
    y_lag = pd.DataFrame(list(map(lambda i: lps.weights.lag_spatial(wq, comp.loc[:, i]), comp.columns))).T
    y_lag.columns, y_lag.index = comp.columns, comp.index
    return comp, y_lag


def spTRS_group_correlation(sptFile, Best_dict, alpha=1):
    comp, y_lag = spTRS_combine_best_group(sptFile, Best_dict)
    # set alpha=0 to cancel the efficiency of neighbor expressions
    comp = comp + alpha * y_lag
    corr, p_val = do_corr(comp, "pearsonr")
    return comp, corr, p_val


def do_corr(comp, method="pearsonr"):
    # normalize the comp
    ct_mean = np.mean(comp, axis=0)
    ct_std = np.std(comp, axis=0)
    comp = (comp - ct_mean) / ct_std
    from scipy.stats import pearsonr, spearmanr
    if method == "pearsonr":
        f = pearsonr
    else:
        f = spearmanr
    corr = np.zeros((comp.shape[1], comp.shape[1]))
    p_val = np.zeros((comp.shape[1], comp.shape[1]))
    for x in range(p_val.shape[0]):
        for y in range(p_val.shape[1]):
            corr[x, y], p_val[x, y] = f(comp.iloc[:, x], comp.iloc[:, y])
    corr = pd.DataFrame(corr, index=comp.columns, columns=comp.columns)
    p_val = pd.DataFrame(p_val, index=comp.columns, columns=comp.columns)
    return corr, p_val


def Show_self_correlation(corr, folder):
    import seaborn as sns
    fig, axs = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    sns.heatmap(corr, cmap="viridis", ax=axs)
    fig.savefig(folder+'/corr1.eps', dpi=400)


def Generate_New_idents(sptFile, Best_dict):
    cell_types = list(Best_dict.keys())
    sptinfo = sptInfo(sptFile)
    dct_mtds = sptinfo.get_spmatrix()
    adatas = list(map(lambda dct: Load_spt_to_AnnData(sptFile, count=dct), dct_mtds))
    adatas = dict(zip(dct_mtds, adatas))
    new_idents = pd.DataFrame(index=adatas['matrix'].obs_names)
    # set all domain to NA
    new_idents['domain'] = ""
    new_idents["fscore"] = 0
    for cell_type in cell_types:
        Best_item = Best_dict[cell_type]
        fscore = Best_item[0]
        domain = Best_item[1]
        mtds = Best_item[2].split('+')
        dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
        adata = adatas[dct_mtd]
        for spot in new_idents[adata.obs[clu_mtd] == domain].index:
            if fscore > new_idents.loc[spot, 'fscore']:
                new_idents.loc[spot, 'domain'] = cell_type
                new_idents.loc[spot, 'fscore'] = fscore

    new_idents[new_idents['domain'] == ""] = "NA"


def Gen_venn(sptFile, Best_dict, top_genes=200):
    from matplotlib_venn import venn2, venn2_circles
    cell_types = list(Best_dict.keys())
    sptinfo = sptInfo(sptFile)
    dct_mtds = sptinfo.get_spmatrix()
    adatas = list(map(lambda dct: Load_spt_to_AnnData(sptFile, count=dct), dct_mtds))
    adatas = dict(zip(dct_mtds, adatas))

    # load the single cell data
    adata_sc = Load_spt_sc_to_AnnData(sptFile)
    sc.pp.log1p(adata_sc)
    for cell_type in cell_types:
        Best_item = Best_dict[cell_type]
        domain = Best_item[1]
        mtds = Best_item[2].split('+')
        dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
        adata = adatas[dct_mtd]
        sc.pp.log1p(adata)
        sc.tl.rank_genes_groups(adata, groupby=clu_mtd, groups=[domain])
        sp_gene = pd.DataFrame(adata.uns['rank_genes_groups']['names'][:top_genes])
        sp_gene = pd.DataFrame(sp_gene.iloc[:, 0].apply(str.upper))
        sc_gene = pd.DataFrame(adata_sc.uns['rank_genes_groups']['names'][:top_genes])
        sc_gene = pd.DataFrame(sc_gene.iloc[:, 0].apply(str.upper))
        f = venn2(
            subsets=[set(sp_gene.iloc[:,0]), set(sc_gene.iloc[:,0])],
            set_labels=('domain:' + str(domain), cell_type),
            set_colors=("#098154", "#c72e29"),
            alpha=0.6,  # 透明度
            normalize_to=1.0)

        f.get_patch_by_id('10').set_edgecolor('red')  # 左圈外框颜色
        f.get_patch_by_id('10').set_linestyle('--')  # 左圈外框线型
        f.get_patch_by_id('10').set_linewidth(2)  # 左圈外框线宽
        f.get_patch_by_id('01').set_edgecolor('green')  # 右圈外框颜色
        f.get_patch_by_id('11').set_edgecolor('blue')  # 中间圈外框颜色
        # plt.show()

        plt.annotate('intersection genes',
                     color='black',
                     xy=f.get_label_by_id('11').get_position() + np.array([0, 0.05]),
                     xytext=(20, 80),
                     ha='center', textcoords='offset points',
                     bbox=dict(boxstyle='round,pad=0.5', fc='grey', alpha=0.6),
                     arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3,rad=-0.5', color='black')
                     )
        plt.show()


def Gen_venn3(adata_sp, adata_sc1, adata_sc2, domain, clu_mtd, cell_type1, cell_type2, folder, top_genes=200):
    from matplotlib_venn import venn3, venn3_circles
    sc.pp.log1p(adata_sp)
    sc.pp.log1p(adata_sc1)
    sc.pp.log1p(adata_sc2)
    sc.tl.rank_genes_groups(adata_sp, groupby=clu_mtd, groups=[domain])
    sc.tl.rank_genes_groups(adata_sc1, groupby="annotation", groups=[cell_type1])
    sc.tl.rank_genes_groups(adata_sc2, groupby="annotation", groups=[cell_type2])
    sp_gene = pd.DataFrame(adata_sp.uns['rank_genes_groups']['names'][:top_genes])
    sp_gene = pd.DataFrame({"names": sp_gene.iloc[:, 0].apply(str.upper),
                            "scores": pd.DataFrame(adata_sp.uns['rank_genes_groups']['scores'][:top_genes]).iloc[:, 0]})
    sp_gene.index = sp_gene['names']
    sc_gene1 = pd.DataFrame(adata_sc1.uns['rank_genes_groups']['names'][:top_genes])
    sc_gene1 = pd.DataFrame({"names": sc_gene1.iloc[:, 0].apply(str.upper),
                             "scores": pd.DataFrame(adata_sc1.uns['rank_genes_groups']['scores'][:top_genes]).iloc[:, 0]})
    sc_gene1.index = sc_gene1['names']
    sc_gene2 = pd.DataFrame(adata_sc2.uns['rank_genes_groups']['names'][:top_genes])
    sc_gene2 = pd.DataFrame({"names": sc_gene2.iloc[:, 0].apply(str.upper),
                             "scores": pd.DataFrame(adata_sc2.uns['rank_genes_groups']['scores'][:top_genes]).iloc[:, 0]})
    sc_gene2.index = sc_gene2['names']

    fig, axs = plt.subplots(1, 1, figsize=(3.5, 3.5), constrained_layout=True)
    f = venn3(
        subsets=[set(sp_gene.iloc[:, 0]), set(sc_gene1.iloc[:, 0]), set(sc_gene2.iloc[:, 0])],
        set_labels=('domain:' + str(domain), cell_type1, cell_type2),
        alpha=0.8,  # 透明度
        normalize_to=1.0, ax=axs)

    f.get_patch_by_id('111').set_linestyle('--')  # 左圈外框线型
    f.get_patch_by_id('111').set_edgecolor('white')  # 中间圈外框颜色
    inter_genes = (set(sp_gene.iloc[:, 0]).intersection(set(sc_gene1.iloc[:, 0])).intersection(set(sc_gene2.iloc[:, 0])))
    inter_genes = pd.DataFrame(list(inter_genes), columns=['names'])
    inter_genes['scores'] = np.mean(np.array([list(sp_gene.loc[list(inter_genes['names']), 'scores']),
                                    list(sc_gene1.loc[inter_genes['names'], 'scores']),
                                    list(sc_gene2.loc[inter_genes['names'], 'scores'])]), axis=0)
    inter_genes = inter_genes.sort_values(by=['scores'], ascending=False)
    axs.annotate('intersection genes',
                 color='black',
                 xy=f.get_label_by_id('111').get_position() + np.array([0, 0.05]),
                 xytext=(20, 100),
                 ha='center', textcoords='offset points',
                 bbox=dict(boxstyle='round,pad=0.5', fc='grey', alpha=0.7),
                 arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3,rad=-0.5', color='black')
                 )
    fig.savefig(folder+"/"+cell_type1+"_"+cell_type2+".eps", dpi=400)
    return inter_genes


def distance(a, b):
    """
    计算a,b两点的直线距离
    :param a: a点坐标（大地坐标 m）
    :param b: b点坐标（大地坐标 m）
    :return: dis 直线距离
    """
    p1 = np.array(a)
    p2 = np.array(b)

    dis = sum(abs(p2 - p1))
    return dis


def Alpha_Shape6(points: pd.DataFrame):
    # 针对6邻居spot阵列的轮廓提取
    x = points['array_row']
    y = points['array_col']
    count = len(x)
    edge_x = []
    edge_y = []
    single_x = []
    single_y = []
    # 首先筛选边缘点
    for i in range(count):
        dis = abs(x - x[i]) + abs(y - y[i])
        neighbors = sum(dis == 2)
        if neighbors >= 6:
            continue
        elif neighbors <= 0:
            single_x.append(x[i])
            single_y.append(y[i])
        else:
            edge_x.append(x[i])
            edge_y.append(y[i])
    # 然后判断多边形数量，以及边界顺序
    # 首先建立点和边
    G = nx.Graph()
    bound_len = len(edge_x)
    nodes = range(bound_len)
    G.add_nodes_from(nodes)

    x = np.array(edge_x.copy())
    y = np.array(edge_y.copy())
    dis = np.array(list(map(lambda i: abs(x - x[i]) + abs(y - y[i]), range(bound_len))))
    dis0 = np.where((dis <= 6) * (dis > 0), dis, 0)
    G = nx.from_numpy_array(dis0)

    nx.is_eulerian(G)
    for i in range(bound_len):
        dis.append(list(abs(x - x[i]) + abs(y - y[i])))
        neighbors = np.where(dis <= 4 and dis > 0)
        for j in range(len(neighbors)):
            G.add_edge(i, neighbors[j])
    cc = nx.connected_components(G)

    return edge_x, edge_y


def Alpha_Shape(points: pd.DataFrame, alpha):
    x = points['array_row']
    y = points['array_col']
    count = len(x)
    i = 0
    temp_i = i
    edge_x = []
    edge_y = []
    edge = []
    edge_len = 0
    while i < count:
        i_range = np.array([t for t, v in enumerate(x < x[i] + 2 * alpha) if v == True])
        i_range = i_range[[t for t, v in enumerate(points.loc[i_range, 'array_row'] > x[i] - 2 * alpha) if v == True]]
        i_range = i_range[[t for t, v in enumerate(points.loc[i_range, 'array_col'] < y[i] + 2 * alpha) if v == True]]
        i_range = i_range[[t for t, v in enumerate(points.loc[i_range, 'array_col'] > y[i] - 2 * alpha) if v == True]]

        """ i_range 是可能与当前点的连线组成边界线的备选点集"""
        for k in i_range:
            # 避免重复选点（一个点不能组成三个边界线，最多只能组成两条边界）
            if edge_x.count(x[k]) > 0 and edge_y.count(y[k]) > 0:
                continue

            # 计算 当前点 与 备选点 的距离
            dis = distance((x[i], y[i]), (x[k], y[k]))

            # 因为当前点 和 备选点 要在一个圆上，所有距离不能大于圆直径, 或者 当前点 与 备选点 重合
            if dis > 2 * alpha or dis == 0:
                continue

            # 当前点与备选点的连线 L 线段中点
            center_x = (x[k] + x[i]) * 0.5
            center_y = (y[k] + y[i]) * 0.5

            # L 的 方向向量
            direct_x = x[k] - x[i]
            direct_y = y[k] - y[i]

            #  L 的 法向量K 及其单位化
            nomal_x = -direct_y / math.sqrt(pow(direct_x, 2) + pow(direct_y, 2))
            nomal_y = direct_x / math.sqrt(pow(direct_x, 2) + pow(direct_y, 2))

            # 圆心到直线L 距离
            disOfcenter = math.sqrt(pow(alpha, 2) - 0.25 * pow(dis, 2))

            if nomal_x != 0:
                a = math.sqrt(pow(disOfcenter, 2) / (1 + pow(nomal_y / nomal_x, 2)))
                b = a * (nomal_y / nomal_x)
            else:
                a = 0
                b = alpha

            # 圆心 坐标
            cicle1_x = center_x + a
            cicle1_y = center_y + b

            cicle2_x = center_x - a
            cicle2_y = center_y - b

            # 检测圆内是否有其他点
            b1 = False
            b2 = False

            # 筛选以圆心R1为质点，边长为 2*alpha 的方形区域内的点集
            inSquare = np.array([t for t, v in enumerate(x < cicle1_x + alpha) if v == True])
            inSquare = inSquare[
                [t for t, v in enumerate(points.loc[inSquare, 'array_row'] > cicle1_x - alpha) if v == True]]
            inSquare = inSquare[
                [t for t, v in enumerate(points.loc[inSquare, 'array_col'] < cicle1_y + alpha) if v == True]]
            inSquare = inSquare[
                [t for t, v in enumerate(points.loc[inSquare, 'array_col'] > cicle1_y - alpha) if v == True]]
            if len(inSquare) != 0:
                for j in inSquare:
                    if j == i or j == k:  # 点集内的点 除去当前点i 和 备选点k
                        continue
                    else:
                        d = distance((x[j], y[j]), (cicle1_x, cicle1_y))  # 计算备选点k与点集内的点的距离
                        if d <= alpha:  # 圆内有点 ，跳过
                            b1 = True
                            break
            # 筛选 圆R2
            inSquare = np.array([t for t, v in enumerate(x < cicle2_x + alpha) if v == True])
            inSquare = inSquare[
                [t for t, v in enumerate(points.loc[inSquare, 'array_row'] > cicle2_x - alpha) if v == True]]
            inSquare = inSquare[
                [t for t, v in enumerate(points.loc[inSquare, 'array_col'] < cicle2_y + alpha) if v == True]]
            inSquare = inSquare[
                [t for t, v in enumerate(points.loc[inSquare, 'array_col'] > cicle2_y - alpha) if v == True]]
            if len(inSquare) != 0:
                for j in inSquare:  # 与原两个点的坐标点一样
                    if j == i or j == k or distance((x[j], y[j]), (x[i], y[i])) == 0:
                        continue
                    else:
                        d = distance((x[j], y[j]), (cicle2_x, cicle2_y))
                        if d <= alpha:  # 圆内有点 ，跳过
                            b2 = True
                            break

            # 圆1 或 圆2 内没有其他的点，则备选点k是边界点
            if b1 == False or b2 == False:
                if edge_x.count(x[i]) < 1 or edge_y.count(y[i]) < 1:  # 当前点未加入边界点集
                    edge_x.append(x[i])
                    edge_y.append(y[i])
                    edge.append(i)
                if edge_x.count(x[k]) < 1 or edge_y.count(y[k]) < 1:  # 备选点k已未加入边界点集
                    edge_x.append(x[k])
                    edge_y.append(y[k])
                    edge.append(k)

                temp_k = k
                break
        # if edge_len < len(edge_x):  # 跳转到新的边界点
        #     temp_i = i
        #     i = temp_k
        #     edge_len = len(edge_x)
        # else:
        temp_i = i
        i = i + 1

    return edge_x, edge_y


def Create_Boundary(adata: ad.AnnData, clu_mtd: str):
    idents = adata.obs[clu_mtd]
    geometry = [xy for xy in zip(adata.obs['array_row'], adata.obs['array_col'])]
    # first calculate global moransI
    ad_gdf = GeoDataFrame(adata.obs[clu_mtd], geometry=geometry, copy=True)

    wq = lps.weights.Queen.from_dataframe(ad_gdf)
    wq.transform = 'r'
