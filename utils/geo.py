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
        adata.obs['p_moransI_'+comp] = spots
        adata0 = adata[adata.obs['p_moransI_'+comp] > 0]
        # calculate the purity with the segmentation
        pred = adata0.obs[clu]
        true = adata0.obs['p_moransI_'+comp]
        ps = purity_score(true, pred)
        ari = adjusted_rand_score(true, pred)
        nmi = normalized_mutual_info_score(true, pred)
        pf =pd.DataFrame(data=[MI.I, ari, nmi, ps], index=['moransI', 'ARI', 'NMI', 'purity'], columns=[comp])
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
                Fscore = (1 + beta**2) * (precision * recall) / (beta ** 2 * (precision + recall))
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
                mtd_name = "{0}_{1}_{2}".format(spmat, dcv, clu)
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
    return Best_dict


def Show_best_group(adata:ad.AnnData, dcv_mtd, clu_mtd, cell_type, domain, figname):
    abd_name = dcv_mtd + '_' + clu_mtd + '_' + cell_type
    adata.obs[abd_name] = adata.obsm[dcv_mtd][cell_type]
    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    sc.pl.spatial(adata, color=clu_mtd, ax=axs[0])
    sc.pl.spatial(adata, color=abd_name, ax=axs[1])
    axs[1].set_title(dcv_mtd)
    fig.suptitle("find {0} locate in {1} using methods chain:{2}".format(cell_type, str(domain), abd_name))
    fig.savefig(figname)
