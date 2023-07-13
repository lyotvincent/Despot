import esda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from geopandas import GeoDataFrame
import libpysal as lps
from shapely.geometry import Point
from utils.io import *
from utils.purity import purity_score
from utils.preprocess import *
from vsl.visualization import Show_intersection_genes, Show_matched_domains, Show_matched_genes, Show_prop_domains, \
    Ax_prop_domains, Ax_expr_spgenes
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score
from scipy.stats import hypergeom
from statistics import mean, stdev
from math import sqrt
import networkx as nx
import seaborn as sns


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
        # y_lag = lps.weights.lag_spatial(wq, y)
        # then calculate local moransI
        li = esda.moran.Moran_Local(y, wq)
        # using lisa cluster to find hot point, cold point, doughnut and diamond
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


def suitable_location(adata, perf, clu, beta=1, greedy=2):
    # find suitable location of cell-type in tissue
    cell_type = perf.index
    best_match = pd.DataFrame(columns=['F-score', 'domain'])

    for ct in cell_type:
        # get sig. spot
        spots = adata.obs['p_moransI_' + ct]
        clusters = adata.obs[clu]
        clu_num = np.unique(adata.obs[clu])
        temp = pd.concat([spots, clusters], axis=1)
        domain_list = []  # suitable domains
        # Use greedy strategy to integrate other domains
        max_Fscore = 0
        bst_cl = None
        for g in range(greedy):
            if g > len(clu_num):  # stop if greedy num larger than cluster num
                break
            for cl in clu_num:
                clu_spot = temp[clu] == cl
                for i in domain_list:
                    clu_spot += temp[clu] == i
                temp0 = temp[clu_spot]

                # calculate the F1-score of proper domains
                precision = len(temp0[temp0['p_moransI_' + ct] == 1]) / len(temp0)
                temp1 = temp[temp['p_moransI_' + ct] == 1]
                clu_spot1 = temp1[clu] == cl
                for i in domain_list:
                    clu_spot1 += temp1[clu] == i
                recall = len(temp1[clu_spot1]) / len(temp1)
                if precision == 0 and recall == 0:
                    Fscore = 0
                else:
                    Fscore = (2 * beta) * (precision * recall) / (precision + beta * recall)
                if Fscore > max_Fscore:
                    max_Fscore = Fscore
                    bst_cl = cl
            if bst_cl is not None:
                domain_list.append(bst_cl)
            bst_cl = None
        best_match = pd.concat([best_match, pd.DataFrame([max_Fscore, tuple(domain_list)],
                                                         index=['F-score', 'domain']).T])
    best_match.index = cell_type
    return best_match

# # set the single-cell marker background
# adata_sc = Load_spt_sc_to_AnnData(sptFile)
# sc.pp.normalize_total(adata_sc)
# sc.pp.log1p(adata_sc)
# sc.tl.rank_genes_groups(adata_sc, groupby='annotation', method='wilcoxon')
#
# cell_types = np.unique(adata_sc.obs['annotation'])
# sc_DEgenes = {}
# for ct in cell_types:
#     print("calculating {0}".format(ct))
#     DE_genes = pd.DataFrame(adata_sc.uns['rank_genes_groups']['names'])
#     DE_genes = pd.DataFrame({"names": DE_genes.loc[:, ct].apply(str.upper),
#                              "scores": pd.DataFrame(adata_sc.uns['rank_genes_groups']['scores']).loc[:, ct],
#                              "pvals": pd.DataFrame(adata_sc.uns['rank_genes_groups']['pvals']).loc[:, ct],
#                              "pvals_adj": pd.DataFrame(adata_sc.uns['rank_genes_groups']['pvals_adj']).loc[:, ct],
#                              "logfoldchanges": pd.DataFrame(adata_sc.uns['rank_genes_groups']['logfoldchanges']).loc[:, ct]
#                              })
#     DE_genes = DE_genes[DE_genes['pvals_adj'] > 0.05]
#     DE_genes = DE_genes[DE_genes['scores'] > 0]
#     sc_DEgenes[ct] = list(DE_genes['names'])


# calculate the sp_DEgenes
# sc.pp.normalize_total(adata)
# sc.pp.log1p(adata)
# for clu in clu_mtds:
#     sc.tl.rank_genes_groups(adata, groupby=clu, method='wilcoxon')
#     sp_DEgenes = {}
#     clu_name = np.unique(adata.obs[clu])
#     for cl in clu_name:
#         print("calculating {0}".format(cl))
#         DE_genes = pd.DataFrame(adata_sc.uns['rank_genes_groups']['names'])
#         DE_genes = pd.DataFrame({"names": DE_genes.loc[:, cl].apply(str.upper),
#                                  "scores": pd.DataFrame(adata_sc.uns['rank_genes_groups']['scores']).loc[:, cl],
#                                  "pvals": pd.DataFrame(adata_sc.uns['rank_genes_groups']['pvals']).loc[:, cl],
#                                  "pvals_adj": pd.DataFrame(adata_sc.uns['rank_genes_groups']['pvals_adj']).loc[
#                                               :, cl],
#                                  "logfoldchanges": pd.DataFrame(
#                                      adata_sc.uns['rank_genes_groups']['logfoldchanges']).loc[:, cl]
#                                  })
#         DE_genes = DE_genes[DE_genes['pvals_adj'] > 0.05]
#         DE_genes = DE_genes[DE_genes['scores'] > 0]
#         sp_DEgenes[cl] = list(DE_genes['names'])
#
# for ct in cell_types:
#     for cl in clu_mtds:
#         sc_genes = sc_DEgenes[ct]
#         sp_genes = sp_DEgenes[clu]
#         overlaps = list(set(sc_genes).intersection(set(sp_genes)))
#         precision = len(overlaps) / len(sc_genes)
#         recall = len(overlaps) / len(sp_genes)
#         marker_score = (precision * recall) / (precision + recall)


def spTRS_findgroups(sptFile, beta=1, greedy=1, spmat_subset=None, clu_subset=None, dcv_subset=None):
    sptinfo = sptInfo(sptFile)
    platform = sptinfo.configs['platform'][0]
    spmats = sptinfo.get_spmatrix()
    if spmat_subset:
        spmats = list(set(spmats).intersection(spmat_subset))
    Best_Fscore = {}

    for spmat in spmats:
        adata = Load_spt_to_AnnData(sptFile, spmat, platform=platform, loadDeconv=True)
        clu_mtds = sptinfo.get_clu_methods(spmat)
        dcv_mtds = sptinfo.get_dcv_methods(spmat)
        if clu_subset:
            clu_mtds = list(set(clu_mtds).intersection(clu_subset))
        if dcv_subset:
            dcv_mtds = list(set(dcv_mtds).intersection(dcv_subset))
        print(clu_mtds)
        print(dcv_mtds)
        for dcv in dcv_mtds:
            for clu in clu_mtds:
                if clu == 'ground_truth':
                    continue
                mtd_name = "{0}+{1}+{2}".format(spmat, dcv, clu)
                print("Finding Groups in the combination of {0}".format(mtd_name))
                perf = p_moransI_purity(adata, dcv_mtd=dcv, clu=clu)
                perf = perf[perf['moransI'] > 0.5]
                BST = suitable_location(adata, perf, clu, beta, greedy)
                Best_Fscore[mtd_name] = BST
    return Best_Fscore


def spTRS_Find_bestGroup(sptFile, beta=1, greedy=1, spmat_subset=None, clu_subset=None, dcv_subset=None):
    Best_Fscore = spTRS_findgroups(sptFile, beta, greedy, spmat_subset, clu_subset, dcv_subset)
    Best_dict = {}

    # find the best method chains
    for mtd in Best_Fscore.keys():
        mtd_Fscore = Best_Fscore[mtd]
        # for each cell-type, find the cluster which has the highest Fscore
        for ct in mtd_Fscore.index:
            if ct not in Best_dict.keys():
                Best_dict[ct] = (mtd_Fscore.loc[ct, 'F-score'], mtd_Fscore.loc[ct, 'domain'], mtd)
            else:
                if mtd_Fscore.loc[ct, 'F-score'] > Best_dict[ct][0]:
                    Best_dict[ct] = (mtd_Fscore.loc[ct, 'F-score'], mtd_Fscore.loc[ct, 'domain'], mtd)

    # check the domains if we need to extend
    sptinfo = sptInfo(sptFile)
    platform = sptinfo.configs['platform'][0]
    for ct in Best_dict.keys():
        f1score, domains, mtd = Best_dict[ct]
        spmat, dcv, clu = mtd.split('+')
        adata = Load_spt_to_AnnData(sptFile, spmat, platform=platform, loadDeconv=True)
        if type(domains) == int:
            domain_list = [domains]
        else:
            domain_list = [d for d in domains]
        perf = p_moransI_purity(adata, dcv_mtd=dcv, clu=clu)
        perf = perf[perf['moransI'] > 0.5]
        spots = adata.obs['p_moransI_' + ct]
        clusters = adata.obs[clu]
        clu_num = np.unique(adata.obs[clu])
        temp = pd.concat([spots, clusters], axis=1)
        max_Fscore = f1score
        bst_cl = None
        for i in clu_num:
            for cl in clu_num:
                clu_spot = temp[clu] == cl
                for i in domain_list:
                    clu_spot += temp[clu] == i
                temp0 = temp[clu_spot]
                precision = len(temp0[temp0['p_moransI_' + ct] == 1]) / len(temp0)

                temp1 = temp[temp['p_moransI_' + ct] == 1]
                clu_spot1 = temp1[clu] == cl
                for i in domain_list:
                    clu_spot1 += temp1[clu] == i
                recall = len(temp1[clu_spot1]) / len(temp1)
                if precision == 0 and recall == 0:
                    Fscore = 0
                else:
                    Fscore = (2 * beta) * (precision * recall) / (precision + beta * recall)
                if Fscore > max_Fscore:
                    max_Fscore = Fscore
                    bst_cl = cl
            if bst_cl is not None:
                domain_list.append(bst_cl)
            bst_cl = None
            if max_Fscore > f1score:
                f1score = max_Fscore
            else:
                break
        Best_dict[ct] = (f1score, tuple(domain_list), mtd)

    # Save_spt_from_BestDict(sptFile, Best_dict)
    # find markers for each region
    # for cell_type in Best_dict.keys():
    #     region = Best_dict[cell_type][1]
    #     mtd = Best_dict[cell_type][2]
    #     clu_mtd = mtd.split('+')[2]
    #     dcv_mtd = mtd.split('+')[1]
    #     dct_mtd = mtd.split('+')[0]
    #     adata = Load_spt_to_AnnData(sptFile, count=dct_mtd)
    #     sc.tl.filter_rank_genes_groups(adata, key=region, groupby=clu_mtd)
    return Best_dict, Best_Fscore


def Pipline_findgroups(sptFile, pipline, beta=1, greedy=1):
    sptinfo = sptInfo(sptFile)
    platform = sptinfo.configs['platform'][0]
    spmat = 'matrix'
    clu_mtds = sptinfo.get_clu_methods(spmat)
    dcv_mtds = sptinfo.get_dcv_methods(spmat)
    BST = None
    adata = Load_spt_to_AnnData(sptFile, spmat, platform=platform, loadDeconv=True)
    if pipline in ['CARD']:  # CARD uses the dominant cell-type as clusters
        cell_prop = adata.obsm[pipline]
        dominant_cell = cell_prop.idxmax(axis=1).astype('category')
        adata.obs[pipline] = dominant_cell
        clu_mtds = clu_mtds + ['CARD']
    for dcv in dcv_mtds:
        for clu in clu_mtds:
            if clu == pipline and dcv == pipline:
                perf = p_moransI_purity(adata, dcv_mtd=dcv, clu=clu)
                perf = perf[perf['moransI'] > 0.5]
                BST = suitable_location(adata, perf, clu, beta, greedy)
    return BST


def Fscore_Comparison(Best_dict: dict, Best_Fscore):
    # find f-score about ground_truth for each pipeline
    ground_truth_max = pd.DataFrame()
    Seurat_max = pd.DataFrame()
    Giotto_max = pd.DataFrame()
    # spTRS compare with Seurat and Giotto
    for mtd in Best_Fscore.keys():
        clu_mtd = mtd.split('+')[2]
        if clu_mtd == 'ground_truth':
            group = Best_Fscore[mtd]
            ground_truth_max = pd.concat([ground_truth_max, pd.DataFrame(group['F-score']).T])
        dcv_mtd = mtd.split('+')[1]
        dct_mtd = mtd.split('+')[0]
        if dct_mtd == 'matrix' and clu_mtd == 'Seurat' and dcv_mtd == 'Seurat':
            group = Best_Fscore[mtd]
            Seurat_max = pd.concat([Seurat_max, pd.DataFrame(group['F-score']).T])
        if dct_mtd == 'matrix' and clu_mtd == 'Giotto' and dcv_mtd == 'Giotto':
            group = Best_Fscore[mtd]
            Giotto_max = pd.concat([Giotto_max, pd.DataFrame(group['F-score']).T])
    ground_truth_max = pd.DataFrame(ground_truth_max.max()).T
    ground_truth_max.index = ['ground_truth']
    Seurat_max = pd.DataFrame(Seurat_max.max()).T
    Seurat_max.index = ['Seurat']
    ground_truth_max = pd.concat([ground_truth_max, Seurat_max])
    Giotto_max = pd.DataFrame(Giotto_max.max()).T
    Giotto_max.index = ['Giotto']
    ground_truth_max = pd.concat([ground_truth_max, Giotto_max])
    cell_type = []
    fscore = []
    for key in Best_dict.keys():
        cell_type.append(key)
        fscore.append(Best_dict[key][0])
    spTRS_max = pd.DataFrame(data=fscore, index=cell_type, columns=['Despot']).T
    ground_truth_max = pd.concat([ground_truth_max, spTRS_max])

    # find f-score about ground_truth for each pipeline
    undct_max = pd.DataFrame()
    dct_max = pd.DataFrame()
    for mtd in Best_Fscore.keys():
        dct_mtd = mtd.split('+')[0]
        if dct_mtd == 'matrix':
            group = Best_Fscore[mtd]
            undct_max = pd.concat([undct_max, pd.DataFrame(group['F-score']).T])
        else:
            group = Best_Fscore[mtd]
            dct_max = pd.concat([dct_max, pd.DataFrame(group['F-score']).T])
    undct_max = pd.DataFrame(undct_max.max()).T
    undct_max.index = ['matrix']
    ground_truth_max = pd.concat([ground_truth_max, undct_max])
    dct_max = pd.DataFrame(dct_max.max()).T
    dct_max.index = ['decontaminated']
    ground_truth_max = pd.concat([ground_truth_max, dct_max])
    # comparison among vae, es and na
    vae_max = pd.DataFrame()
    es_max = pd.DataFrame()
    na_max = pd.DataFrame()
    for mtd in Best_Fscore.keys():
        dcv_mtd = mtd.split('+')[1]
        if dcv_mtd in ["StereoScope", "spacexr", "SPOTlight_vae"]:
            group = Best_Fscore[mtd]
            vae_max = pd.concat([vae_max, pd.DataFrame(group['F-score']).T])
        elif dcv_mtd in ["StereoScope_es", "spacexr_es", "SPOTlight_es"]:
            group = Best_Fscore[mtd]
            es_max = pd.concat([es_max, pd.DataFrame(group['F-score']).T])
        elif dcv_mtd in ['StereoScope_na', 'SPOTlight']:
            group = Best_Fscore[mtd]
            na_max = pd.concat([na_max, pd.DataFrame(group['F-score']).T])
    vae_max = pd.DataFrame(vae_max.max()).T
    vae_max.index = ['vae']
    ground_truth_max = pd.concat([ground_truth_max, vae_max])

    es_max = pd.DataFrame(es_max.max()).T
    es_max.index = ['es']
    ground_truth_max = pd.concat([ground_truth_max, es_max])

    na_max = pd.DataFrame(na_max.max()).T
    na_max.index = ['na']
    pip_res = pd.concat([ground_truth_max, na_max])
    return pip_res


def Fscore_Comparison_pip(sptFile, has_ground_truth=False):
    sptinfo = sptInfo(sptFile)
    map_chains = sptinfo.get_map_chains()
    pip_names = ['Seurat', 'Giotto', 'CARD']
    pip_res = pd.DataFrame(map_chains.iloc[:, 0])
    pip_res.columns = ['Despot']
    for pip in pip_names:
        pip_perf = Pipline_findgroups(sptFile, pip, beta=1, greedy=1)
        pip_res[pip] = pip_perf.iloc[:, 0]
    if has_ground_truth:
        pip_gt, _ = spTRS_Find_bestGroup(sptFile, beta=1, greedy=1, clu_subset=['ground_truth'])
        col = ["cell-type", "F1-score", "domain", "dct", "dcv", "clu"]
        best_df = pd.DataFrame(columns=col)
        cell_types = pip_gt.keys()
        for ct in cell_types:
            Best_item = pip_gt[ct]
            mtds = Best_item[2].split('+')
            dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
            df_item = pd.DataFrame(data=[ct, Best_item[0], Best_item[1], dct_mtd, dcv_mtd, clu_mtd],
                                   index=col, columns=[ct]).T
            best_df = pd.concat([best_df, df_item])
        pip_res['ground_truth'] = best_df['F1-score']
    return pip_res.T

def Show_Comparison(sptFile,folder, figsize=(3.5,3), compare='platform',cell_type=None, title='F1-score', has_ground_truth=False):
    from vsl.palette import Set_palette
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams['font.size'] = 14
    pip_res = Fscore_Comparison_pip(sptFile, has_ground_truth)
    comp_mtds = ['Giotto', 'Seurat', 'CARD', 'Despot']
    if has_ground_truth:
        comp_mtds= ['Giotto', 'Seurat', 'CARD', 'ground_truth', 'Despot']

    full_ct = list(pip_res.columns)
    full_ct.sort()
    palette, cmap = Set_palette(len(full_ct))
    if cell_type is None:
        comp_ct = list(pip_res.columns)
    else:
        comp_ct = cell_type
    comp_ct.sort()
    comp_res = pip_res.loc[comp_mtds, comp_ct]
    x = np.arange(len(comp_mtds)) * 3
    all_width = 2.5  # the width of the bars
    # fig = None
    # if ax is None:
    fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)
    for m in range(len(comp_mtds)):
        comp_mtd = comp_mtds[m]
        measurement = comp_res.loc[comp_mtd, comp_ct]
        measurement = pd.DataFrame([np.nan_to_num(item) for item in measurement], index=comp_res.columns)
        measurement = measurement.sort_index()
        bar_width = all_width / sum(list(measurement.loc[:, 0] > 0))
        multiplier = 0
        for ct in measurement.index:
            if measurement.loc[ct, 0] > 0:
                offset = bar_width * multiplier
                color = palette[full_ct.index(ct)]
                rects = ax.bar(x[m] + offset, measurement.loc[ct, 0], bar_width, color=color)
                ax.bar_label(rects, padding=3, fmt='%.2f', rotation=90, fontsize=10)
                multiplier += 1
    ax.set_title(title, fontsize=16)
    ax.set_xticks(x+all_width/3, comp_mtds)
    ax.set_yticks([.0, .5, 1.0])
    ax.set_yticklabels([.0, .5, 1.0], fontsize=10)
    ax.set_ylim(-0.03, 1.1)
    ax.set_xlim(-0.5, max(x)+all_width)
    ax.spines.top.set_visible(False)
    ax.spines.right.set_visible(False)
    fig.savefig(folder+'/platform_comparison.svg', dpi=400)
    return fig

    # for i in range(len(comp_mtd)):
    #     attribute = comp_mtd[i]
    #     measurement = comp_res.iloc[i, :]
    #     measurement = [np.nan_to_num(item) for item in measurement]
    #     offset = width * multiplier
    #     rects = ax.bar(x + offset, measurement, width, label=attribute, color=cmap[i])
    #     ax.bar_label(rects, padding=3, fmt='%.2f', rotation=90)
    #     multiplier += 1
    #
    # # Add some text for labels, title and custom x-axis tick labels, etc.
    # ax.set_ylabel('F1-score')
    # ax.set_xticks(x + width, comp_ct, rotation=20)
    # ax.legend(loc='upper left', ncol=3)
    # ax.set_ylim(0, 1.1)


# 3D landscape for SCSPs
def Show_3D_landscape(sptFile, folder=None, cell_types=None, sf=None, pipline=None,imgPath=None, alpha=0.1):
    from mpl_toolkits.mplot3d import Axes3D
    from vsl.boundary import boundary_extract, show_edge
    from PIL import Image
    from pylab import ogrid
    from vsl.palette import Set_palette
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    sptinfo = sptInfo(sptFile)
    platform = sptinfo.get_platform()
    print("Platform: {0}".format(platform))
    # handle images
    if sptinfo.get_platform() != 'ST':
        img = Image.open(sptinfo.get_imgPath('low'))
        img0 = np.array(img) / 255
        img0[img0 > 1] = 1
        w, h = img.height, img.width

        # handle coordinates
        coords = sptinfo.get_coords()
        if sf is None:
            sf = sptinfo.get_sf(solution='low')
        x, y = coords.loc[:, 'image_row'] * sf, coords.loc[:, 'image_col'] * sf
        xmin, xmax = int(np.min(x) - 60), int(np.max(x) + 60)
        ymin, ymax = int(np.min(y) - 60), int(np.max(y) + 60)
        x -= xmin
        y -= ymin
        imgX, imgY = ogrid[0:(xmax - xmin), 0:(ymax - ymin)]
    else:
        img = Image.open(imgPath)
        img0 = np.array(img) / 255
        coords = sptinfo.get_coords()
        w, h = img.width, img.height
        sf = img.height / (np.max(coords.loc[:, 'image_col'])+1.5 - np.min(coords.loc[:, 'image_col']))
        x, y = (np.max(coords.loc[:, 'image_col'])+1.5-coords.loc[:, 'image_col'])*sf, (coords.loc[:, 'image_row']-0.2)*sf
        # xmin, xmax = int(np.min(x)), int(np.max(x))
        # ymin, ymax = int(np.min(y)), int(np.max(y))
        # x -= xmin
        # y -= ymin
        imgX, imgY = ogrid[0:img.height, 0:img.width]

    # handle the compared pipline or Despot
    # the full cell-types for SCSPs
    full_ct = list(sptinfo.get_map_chains().index)
    full_ct.sort()
    if pipline is None:
        Best_df = sptinfo.get_map_chains()
        figname = folder + '/landscape.svg'
    else:
        Best_df = Pipline_findgroups(sptFile, pipline, beta=1, greedy=1)
        Best_df['dct'] = 'matrix'
        Best_df['dcv'] = pipline
        Best_df['clu'] = pipline
        figname = folder + '/{0}.svg'.format(pipline)
    if cell_types is None:
        if len(Best_df) > 0:
            cell_types = list(Best_df.index)
            print(cell_types)
        else:
            cell_types = []
    else:
        # select the overlaps of provided cell_types and existing cell_types
        cell_types = list(set(cell_types).intersection(list(Best_df.index)))

    # arrange the cell_types
    spot3Ds = {}
    centroid = []
    cell_types.sort()  # first sort the cell-type by alphabets
    for cell_type in cell_types:
        Best_item = Best_df.loc[cell_type, :]
        domain = Best_item['domain']
        if type(domain) == int:
            domain = (domain,)
        dct_mtd, dcv_mtd, clu_mtd = Best_item['dct'], Best_item['dcv'], Best_item['clu']
        if clu_mtd in ['CARD']:   # if use the dominant cell-type as cluster, we directly generate it
            adata = Load_spt_to_AnnData(sptFile, 'matrix', platform=platform, loadDeconv=True)
            cell_prop = adata.obsm[clu_mtd]
            idents = cell_prop.idxmax(axis=1).astype('category')
        else:
            idents = sptinfo.get_idents(dct_mtd)[clu_mtd]

        spot3D = pd.concat([x, y, idents], axis=1)
        spot3D.columns = ['X', 'Y', 'level']
        # select domain
        selection = np.zeros_like(idents)
        for d in domain:
            selection += (np.array(idents) == d)
        selection = np.array(selection, dtype=bool)
        spot3D = spot3D[selection]
        spot3Ds[cell_type] = spot3D
        # the closer cell-type has a lower layer
        center = (- spot3D['X'].min() + spot3D['Y'].max()) + 2 * (- spot3D['X'].mean() + spot3D['Y'].mean())
        centroid.append(center)
    centroid = pd.DataFrame(centroid, index=cell_types)

    # plot the figure with Axes3D
    fig = plt.figure(figsize=(15, 10))
    ax1 = plt.axes(projection='3d')
    arranged_types = centroid.sort_values(by=0)
    palette, cmap = Set_palette(len(full_ct))
    if len(full_ct) > 20:
        palette = palette + palette
    z = 0  # the Z index
    prev_spot3D = None
    spot3D_layer = {}
    for cell_type in arranged_types.index:
        # print(full_ct.index(cell_type))
        color = palette[full_ct.index(cell_type)]  # get the related cell-type color
        spot3D = spot3Ds[cell_type]
        if len(spot3D_layer) == 0:
            z = 1
            spot3D_layer[z] = [cell_type]
        else:
            need_new_layer = True   # if need new layer
            for layer in spot3D_layer.keys():
                ct_in_layer = spot3D_layer[layer]  # all the cell types in a layer
                put_in_layer = True
                for ct in ct_in_layer: # if new cell type has overlap of existing layer, z=z+1.
                    if len(pd.merge(spot3Ds[ct][['X', 'Y']], spot3D[['X', 'Y']], how='inner')) > 0:
                        put_in_layer = False
                        break
                if put_in_layer:
                    z = layer
                    spot3D_layer[layer].append(cell_type)
                    need_new_layer = False   # cell_type has been put in existing layer, so we don't need a new layer
                    break
            if need_new_layer:
                z = max(list(spot3D_layer.keys())) + 1
                spot3D_layer[z] = [cell_type]
        pts = np.array(spot3D[['X', 'Y']])
        edges, centers = boundary_extract(pts, alpha, err=10e-5)
        spot3D['Z'] = np.ones_like(spot3D.index) * z
        if sptinfo.get_platform() == 'ST':
            s = 30
        else:
            s=10
        ax1.scatter3D(np.array(spot3D['X'], dtype=float),
                      np.array(spot3D['Y'], dtype=float),
                      np.array(spot3D['Z'], dtype=float),
                      color=color, label=cell_type, s=s, alpha=1)
        if sptinfo.get_platform() == 'ST':
            sel_line = np.random.randint(len(spot3D.index), size=max(round(len(spot3D.index) / 10), 5))
        else:
            sel_line = np.random.randint(len(spot3D.index), size=max(round(len(spot3D.index) / 20), 5))
        for item in sel_line:
            x = spot3D['X'][item]
            y = spot3D['Y'][item]
            zz = spot3D['Z'][item]
            ax1.plot([x, x], [y, y], [zz, 0], color=color, alpha=0.5, linewidth=1.25)
        show_edge(edges, ax1, z=0, color=color, linewidth=2, alpha=1, label=None)
    if sptinfo.get_platform() != 'ST':
        ax1.plot_surface(imgX, imgY, np.atleast_2d(0), rstride=10, cstride=10, facecolors=img0[xmin:xmax, ymin:ymax],
                     alpha=0.5, linewidth=0)
    else:
        ax1.plot_surface(imgX, imgY, np.atleast_2d(0), rstride=5, cstride=5, facecolors=img0,
                         alpha=0.5, linewidth=0)
    # ax1.scatter3D(spot3Dbg['X'], spot3Dbg['Y'], spot3Dbg['Z'], c='grey', label='background')
    ax1.set_xticks(ticks=range(w))
    ax1.set_yticks(ticks=range(h))
    ax1.set_zticks(ticks=range(int(z*1.5)))
    ax1.axis('off')
    ax1.legend(loc='lower left', ncol=3, frameon=False,
               bbox_to_anchor=(0.1, -0.05),
               handletextpad=0.3,
               borderpad=0.5,
               borderaxespad=1.05,
               columnspacing=0.7,
               handlelength=0.7,
               fontsize=14)
    fig.savefig(figname, dpi=400)


def Show_best_group(sptFile, cell_type, fscore, domain, mtd_chain: str, figname=None, save=True):
    info = sptInfo(sptFile)
    platform = info.configs['platform'][0]
    if platform != 'ST':
        platform = '10X_Visium'
    mtds = mtd_chain.split('+')
    dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
    abd_name = dct_mtd + '_' + dcv_mtd + '_' + clu_mtd + '_' + cell_type
    adata = Load_spt_to_AnnData(sptFile, h5data=dct_mtd, loadDeconv=True, platform=platform)
    adata.obs[abd_name] = adata.obsm[dcv_mtd][cell_type]
    fig = Show_matched_domains(adata, clu_mtd, dcv_mtd, domain, cell_type, fscore, platform=platform)
    if save:
        fig.savefig(figname, dpi=400)
    print("find {0} locate in {1} using methods chain:"
          "\n{2}, F-score={3}".format(cell_type, str(domain), abd_name, fscore))
    return fig


def Show_bash_best_group(sptFile, folder, Best_dict=None):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    sptinfo = sptInfo(sptFile)
    Best_df = sptinfo.get_map_chains()
    cell_types = []
    if Best_dict is not None:
        cell_types = list(Best_dict.keys())
    else:
        if len(Best_df) > 0:
            cell_types = list(Best_df.index)
            print(cell_types)

    for cell_type in cell_types:
        if Best_dict is not None:
            Best_item = Best_dict[cell_type]
            fscore = Best_item[0]
            domain = Best_item[1]
            mtds = Best_item[2].split('+')
            dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
        else:
            Best_item = Best_df.loc[cell_type, :]
            domain = Best_item['domain']
            fscore = Best_item['F1-score']
            dct_mtd, dcv_mtd, clu_mtd = Best_item['dct'], Best_item['dcv'], Best_item['clu']
        mtd_chain = dct_mtd + '+' + dcv_mtd + '+' + clu_mtd
        fig_cell_type = cell_type.replace('/', '.')  # illegal `/` in file path
        figname = folder + "/" + fig_cell_type + ".svg"
        fig = Show_best_group(sptFile, cell_type, fscore, domain, mtd_chain, figname, save=True)


def Show_row_best_group(sptFile, folder, cell_types, file_name='domains1.svg', imgPath=None, titles=None, cmap='cividis', ap=100):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    sptinfo = sptInfo(sptFile)
    platform = sptinfo.get_platform()
    if platform != 'ST':
        platform = '10X_Visium'
    Best_df = sptinfo.get_map_chains()
    nct = len(cell_types)
    if platform == '10X_Visium':
        fig, axs = plt.subplots(1, nct, figsize=(3 * nct, 4))
    else:
        fig, axs = plt.subplots(1, nct, figsize=(3.5 * nct, 3.2))
    plt.subplots_adjust(wspace=0.05, hspace=0)
    for i in range(nct):
        cell_type = cell_types[i]
        if titles:
            title = titles[i]
        else:
            title = cell_type
        Best_item = Best_df.loc[cell_type, :]
        domain = Best_item['domain']
        f1score = Best_item['F1-score']
        dct_mtd, dcv_mtd, clu_mtd = Best_item['dct'], Best_item['dcv'], Best_item['clu']
        adata = Load_spt_to_AnnData(sptFile, h5data=dct_mtd, loadDeconv=True, platform=platform)
        Ax_prop_domains(adata, clu_mtd, dcv_mtd, domain, cell_type, axs[i],
                        title=title, f1score=f1score, platform=platform, imgPath=imgPath, cmap=cmap, ap=ap)
    figname = folder + "/"+file_name
    fig.savefig(figname, dpi=400)


def Show_row_marker_genes(sptFile, folder, genes, cell_types=None, spmatrix=None,figname='gene.svg', cmap='coolwarm'):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    sptinfo = sptInfo(sptFile)
    Best_df = sptinfo.get_map_chains()
    platform = sptinfo.get_platform()
    print(platform)
    if platform != 'ST':
        platform = '10X_Visium'
    ngenes = len(genes)
    if platform == 'ST':
        fig, axs = plt.subplots(1, ngenes, figsize=(3.5 * ngenes, 3.2))
    else:
        fig, axs = plt.subplots(1, ngenes, figsize=(3 * ngenes, 4))
    plt.subplots_adjust(wspace=0.05, hspace=0)
    for i in range(ngenes):
        gene = genes[i]
        if spmatrix is None:
            cell_type = cell_types[i]
            Best_item = Best_df.loc[cell_type, :]
            dct_mtd, dcv_mtd, clu_mtd = Best_item['dct'], Best_item['dcv'], Best_item['clu']
        else:
            dct_mtd = spmatrix
        adata_sp = Load_spt_to_AnnData(sptFile, h5data=dct_mtd, platform=platform)
        sc.pp.normalize_total(adata_sp)
        sc.pp.log1p(adata_sp)
        Ax_expr_spgenes(adata_sp, gene, ax=axs[i], platform=platform, title=gene, cmap=cmap)
    figname = folder + "/" + figname
    fig.savefig(figname, dpi=400)


def Show_row_scmarkers(sptFile, folder, genes, figname='sc_gene.svg'):
    if not os.path.isdir(folder):
        os.makedirs(folder)
    ngenes = len(genes)
    fig, axs = plt.subplots(1, ngenes, figsize=(3 * ngenes, 2.6))
    plt.subplots_adjust(wspace=0.1, hspace=0)
    adata_sc = Load_spt_sc_to_AnnData(sptFile, h5data='scRNA_seq')
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata_sc = Spatial_Cluster_Analysis(adata_sc)
    for i in range(ngenes):
        gene = genes[i]
        sc.pl.umap(adata_sc, ax=axs[i], color=gene, show=False, size=12, frameon=False)
        axs[i].set_xlabel('')
        axs[i].set_ylabel('')
    figname = folder + "/" + figname
    fig.savefig(figname, dpi=400)


def spTRS_self_correlation(sptFile, alpha=1, method='pearsonr'):
    # using scipy to calculate pearson's r or spearman r
    from scipy.stats import pearsonr, spearmanr
    if method == 'spearmanr':
        corr_mtd = spearmanr
    else:
        corr_mtd = pearsonr
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
                    corr[x, y], p_val[x, y] = corr_mtd(comp.iloc[:, x], comp.iloc[:, y])
            coef = pd.DataFrame(coef, index=comp.columns, columns=comp.columns)

            # Save_spt_from_corrcoef(sptFile, coef, dcv, h5data=h5data)


def spTRS_combine_best_group(sptFile, Best_dict=None):
    sptinfo = sptInfo(sptFile)
    Best_df = sptinfo.get_map_chains()
    cell_types = []
    if Best_dict is not None:
        cell_types = list(Best_dict.keys())
    else:
        if len(Best_df) > 0:
            cell_types = list(Best_df.index)
            print(cell_types)
    abundance = pd.DataFrame()
    sptinfo = sptInfo(sptFile)
    dct_mtds = sptinfo.get_spmatrix()
    adatas = list(map(lambda dct: Load_spt_to_AnnData(sptFile, h5data=dct, loadDeconv=True), dct_mtds))
    adatas = dict(zip(dct_mtds, adatas))
    for cell_type in cell_types:
        if Best_dict is not None:
            Best_item = Best_dict[cell_type]
            domain = Best_item[1]
            mtds = Best_item[2].split('+')
            dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
        else:
            Best_item = Best_df.loc[cell_type, :]
            dct_mtd, dcv_mtd, clu_mtd = Best_item['dct'], Best_item['dcv'], Best_item['clu']
        abundance = pd.concat([abundance, adatas[dct_mtd].obsm[dcv_mtd][cell_type]], axis=1)
    geometry = [Point(xy) for xy in zip(adatas['matrix'].obs['array_row'], adatas['matrix'].obs['array_col'])]
    ad_gdf = GeoDataFrame(abundance, geometry=geometry, copy=True)
    wq = lps.weights.Queen.from_dataframe(ad_gdf)
    wq.transform = 'r'
    comp = abundance
    y_lag = pd.DataFrame(list(map(lambda i: lps.weights.lag_spatial(wq, comp.loc[:, i]), comp.columns))).T
    y_lag.columns, y_lag.index = comp.columns, comp.index
    return comp, y_lag


def spTRS_group_correlation(sptFile, Best_dict=None, alpha=1):
    comp, y_lag = spTRS_combine_best_group(sptFile, Best_dict)
    # set alpha=0 to cancel the efficiency of neighbor expressions
    comp = comp + alpha * y_lag
    corr, p_val = do_corr(comp, "pearsonr")
    return comp, corr, p_val



ct0 = ['Acinar cells', 'Cancer clone A', 'Cancer clone B', 'Endocrine cells', 'RBCs','Tuft cells',
'Fibroblasts','T cells & NK cells', 'mDCs B','Macrophages B','Monocytes','Mast cells', 'mDCs B',
       'Ductal - APOL1 high/hypoxic', 'Ductal - CRISP3 high/centroacinar like','Ductal - MHC Class II','Ductal - terminal ductal like',
       ]

ct1 = [
'CAFs myCAF-like','Cycling_Myeloid','Macrophage','Monocyte','DCs', 'Endothelial Lymphatic LYVE1','Plasmablasts', 'B cells Memory',
'Endothelial CXCL12','Endothelial RGS5','Endothelial ACKR1','PVL Immature','Cancer LumB SC','Cancer Cycling', 'Cancer Basal SC',
    'Mature Luminal','Luminal Progenitors','Myoepithelial',
'B cells Naive','CAFs MSC iCAF-like',
]

ct2 = ['PT(S1-S2)', 'PT(S3)', 'LH(AL)', 'LH(DL)', 'CD-PC', 'CNT', 'DCT']

ct3 = ['Normal Epithelial', 'Cancer Epithelial', 'CAFs', 'T-cells', 'Myeloid']

ct4 =['Cancer Epithelial', 'T-cells', 'CAFs', 'Plasmablasts',
       'Myeloid', 'Endothelial', 'PVL', 'Normal Epithelial','B-cells']

def Despot_pair_correlation(sptFile1, sptFile2, alpha=1, ct=None):
    comp1, y_lag1 = spTRS_combine_best_group(sptFile1)
    comp1 = comp1 + alpha * y_lag1
    comp2, y_lag2 = spTRS_combine_best_group(sptFile2)
    comp2 = comp2 + alpha * y_lag2
    corr, p_val = do_corr(comp1, comp2)
    idx = corr.columns
    diag_val = pd.DataFrame([corr.loc[i,i] for i in idx], index=idx)
    diag_val = diag_val.sort_values(by=0,ascending=False)
    if ct:
        corr = corr[ct]
        corr = corr.loc[np.flip(ct), :]  # index are sn, columns are sc
    else:
        corr = corr[diag_val.index]
        corr = corr.loc[np.flip(diag_val.index), :]  # index are sn, columns are sc
    return corr, p_val


def do_corr(comp, comp0=None, method="pearsonr"):
    # normalize the comp
    from scipy.stats import pearsonr, spearmanr
    if method == "pearsonr":
        f = pearsonr
    else:
        f = spearmanr
    ct_mean = np.mean(comp, axis=0)
    ct_std = np.std(comp, axis=0)
    comp = (comp - ct_mean) / ct_std
    if comp0 is None:
        comp0 = comp
        ct_mean = np.mean(comp0, axis=0)
        ct_std = np.std(comp0, axis=0)
        comp0 = (comp0 - ct_mean) / ct_std
    corr = np.zeros((comp.shape[1], comp0.shape[1]))
    p_val = np.zeros((comp.shape[1], comp0.shape[1]))
    for x in range(p_val.shape[0]):
        for y in range(p_val.shape[1]):
            corr[x, y], p_val[x, y] = f(comp.iloc[:, x], comp0.iloc[:, y])
    corr = pd.DataFrame(corr, index=comp.columns, columns=comp0.columns)
    p_val = pd.DataFrame(p_val, index=comp.columns, columns=comp0.columns)
    return corr, p_val


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    # cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=row_labels)
    if len(col_labels) > 0:
        ax.set_yticks(np.arange(data.shape[0]), labels=col_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True, color='white')

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    # ax.grid(which="minor", color="grey", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """
    from matplotlib import ticker
    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = 0.3

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    def func(x, pos):
        if x < 0:
            return ''
        else:
            return f"{x:.2f}".replace("0.", ".").replace("1.00", "1.0")

    valfmt = ticker.FuncFormatter(func)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(data[i, j] < threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def Show_paired_correlation(corr, folder):
    fig, axs = plt.subplots(1, 1, figsize=(6.5, 5.5), constrained_layout=True)
    im = heatmap(corr, corr.columns, corr.index, ax=axs, vmin=-1,
                    cmap="cividis", cbarlabel="Pearson's Correlation")

    from matplotlib import ticker
    annotate_heatmap(im, size=11)
    axs.set_ylabel('single-nucleus', fontsize=14)
    axs.set_xlabel('single-cell', fontsize=14)
    axs.set_title('Co-occurrence of sc & sn cell-types')
    fig.savefig(folder + '/corr.svg', dpi=400)


def Show_self_correlation(corr, folder, figsize=(7.5, 7.5)):
    fig, axs = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)
    im = heatmap(corr, corr.columns, corr.index, ax=axs, vmin=0,
                    cmap="cividis", cbarlabel="Pearson's Correlation")

    from matplotlib import ticker
    annotate_heatmap(im, size=11)
    axs.set_title('Self correlation of cell-types')
    fig.savefig(folder + '/corr.svg', dpi=400)


def Generate_New_idents(sptFile, Best_dict):
    cell_types = list(Best_dict.keys())
    sptinfo = sptInfo(sptFile)
    dct_mtds = sptinfo.get_spmatrix()
    adatas = list(map(lambda dct: Load_spt_to_AnnData(sptFile, h5data=dct), dct_mtds))
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


# generate venn graph for multimodal integration analysis
def Gen_venn(sptFile, folder, sptFile2=None, Best_dict=None, show_genes=2, cell_filter=None, pval=10e-6, effect=0.5):
    sptinfo = sptInfo(sptFile)
    platform = sptinfo.configs['platform'][0]
    if platform != 'ST':
        platform = '10X_Visium'
    Best_df = sptinfo.get_map_chains()
    cell_types = []
    if Best_dict is not None:
        cell_types = list(Best_dict.keys())
    else:
        if len(Best_df) > 0:
            cell_types = list(Best_df.index)

    dct_mtds = sptinfo.get_spmatrix()
    adatas = list(map(lambda dct: Load_spt_to_AnnData(sptFile, h5data=dct, platform=platform), dct_mtds))
    for i in range(len(adatas)):
        sc.pp.log1p(adatas[i])
    adatas = dict(zip(dct_mtds, adatas))

    # load the single cell data
    adata_sc = Load_spt_sc_to_AnnData(sptFile)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata_sc = Spatial_Cluster_Analysis(adata_sc)
    adata_sc0 = adata_sc.copy()
    if sptFile2 is not None:
        adata_sc1 = Load_spt_sc_to_AnnData(sptFile2)
        adata_sc1 = Spatial_Cluster_Analysis(adata_sc1)
        adata_sc2 = adata_sc1.copy()
    inter_genes = {}
    for cell_type in cell_types:
        if cell_filter is not None and cell_type in cell_filter:
            continue
        if Best_dict is not None:
            Best_item = Best_dict[cell_type]
            domain = Best_item[1]
            mtds = Best_item[2].split('+')
            dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
        else:
            Best_item = Best_df.loc[cell_type, :]
            domain = Best_item['domain']
            dct_mtd, dcv_mtd, clu_mtd = Best_item['dct'], Best_item['dcv'], Best_item['clu']
        adata_sp = adatas[dct_mtd]
        if sptFile2 is not None:
            inter_gene = Gen_venn3(adata_sp, adata_sc, adata_sc2, domain, clu_mtd, cell_type, cell_type, folder)
        else:
            inter_gene = Gen_venn2(adata_sp, adata_sc, domain, clu_mtd, cell_type, folder, pval, effect)
        marker_genes = list(inter_gene.iloc[0:show_genes, 0])
        # for marker in marker_genes:
        #     figname = folder + '/' + marker + "_" + cell_type.replace('/', '|') + "_expr.eps"
        #     Show_matched_genes(adata_sp=adata_sp,
        #                        adata_sc=adata_sc0,
        #                        genes=marker,
        #                        platform=platform,
        #                        figname=figname)
        inter_genes[cell_type] = inter_gene
    return inter_genes


def Gen_venn2(adata_sp, adata_sc, domain, clu_mtd, cell_type, folder, pval=10e-6, effect=0.5):
    # generating venn graph for intersection genes from scRNA-seq and ST
    from matplotlib_venn import venn2
    adata_sp.obs[clu_mtd] = adata_sp.obs[clu_mtd].astype('str')

    if type(domain) == tuple:
        groups = [str(i) for i in domain]
    else:
        groups = [str(domain)]
    sp_gene = Gen_DEgenes(adata_sp, clu_mtd, groups, pval=pval, effect=effect, only_positive=True)
    sc_gene = Gen_DEgenes(adata_sc, groupby="annotation", groups=[cell_type], pval=pval, effect=effect,only_positive=True)

    if len(sp_gene) > 0 and len(sc_gene) > 0:
        inter_genes = (set(sp_gene.iloc[:, 0]).intersection(set(sc_gene.iloc[:, 0])))
        if len(inter_genes) > 0:
            inter_genes = pd.DataFrame(list(inter_genes), columns=['names'])
            inter_genes['scores'] = np.mean(np.array([list(sp_gene.loc[list(inter_genes['names']), 'scores']),
                                                      list(sc_gene.loc[inter_genes['names'], 'scores'])]), axis=0)
            inter_genes = inter_genes.sort_values(by=['scores'], ascending=False)
            bg_genes = Gen_Background_genes(adata_sp, adata_sc)
            mia = MIA_analysis(len(inter_genes), len(bg_genes), len(sc_gene), len(sp_gene))
            # fig, axs = plt.subplots(1, 1, figsize=(3.5, 3.5), constrained_layout=True)
            # f = venn2(
            #     subsets=[set(sp_gene.iloc[:, 0]), set(sc_gene.iloc[:, 0])],
            #     set_labels=('domain:' + str(domain), cell_type),
            #     alpha=0.8,  # Transparency
            #     normalize_to=1.0, ax=axs)
            # f.get_patch_by_id('11').set_linestyle('--')  # Middle ring outer frame line type
            # f.get_patch_by_id('11').set_edgecolor('white')  # Middle ring outer frame line color
            # axs.annotate('MIA: %.3f' % mia,
            #              color='black',
            #              xy=f.get_label_by_id('11').get_position() + np.array([0, 0.05]),
            #              xytext=(20, 100),
            #              ha='center', textcoords='offset points',
            #              bbox=dict(boxstyle='round,pad=0.5', fc='grey', alpha=0.7),
            #              arrowprops=dict(arrowstyle='-|>', connectionstyle='arc3,rad=-0.5', color='black')
            #              )
            # figname = folder + "/" + cell_type.replace("/", ".") + "_.svg"
            # fig.savefig(figname, dpi=400)
        else:
            print('No overlapped genes detected, return empty DataFrame.')
            inter_genes = pd.DataFrame(columns=['names'])
    else:
        print('No DE genes detected, return empty DataFrame.')
        inter_genes=pd.DataFrame(columns=['names'])
    return inter_genes


def Gen_venn3(adata_sp, adata_sc1, adata_sc2, domain, clu_mtd, cell_type1, cell_type2, folder, pval=10e-6, effect=0.2):
    from matplotlib_venn import venn3
    adata_sp.obs[clu_mtd] = adata_sp.obs[clu_mtd].astype('str')
    if type(domain) == tuple:
        groups = [str(i) for i in domain]
    else:
        groups = [str(domain)]

    sp_gene = Gen_DEgenes(adata_sp, clu_mtd, groups, pval=pval, effect=0.8)
    sc_gene1 = Gen_DEgenes(adata_sc1, "annotation", [cell_type1], pval=pval, effect=effect)
    sc_gene2 = Gen_DEgenes(adata_sc2, "annotation", [cell_type2], pval=pval, effect=effect)

    fig, axs = plt.subplots(1, 1, figsize=(3.5, 3.5), constrained_layout=True)
    f = venn3(
        subsets=[set(sp_gene.iloc[:, 0]), set(sc_gene1.iloc[:, 0]), set(sc_gene2.iloc[:, 0])],
        set_labels=('domain:' + str(domain), cell_type1, cell_type2),
        alpha=0.5,  # Transparency
        normalize_to=1.0, ax=axs)

    f.get_patch_by_id('111').set_linestyle('--')  # line style of intersection genes part
    f.get_patch_by_id('111').set_edgecolor('white')  # 中间圈外框颜色
    inter_genes = (set(sp_gene.iloc[:, 0]).intersection(set(
        sc_gene1.iloc[:, 0])).intersection(set(sc_gene2.iloc[:, 0])))
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
    fig.savefig(folder + "/" + cell_type1 + "_" + cell_type2 + ".svg", dpi=400)
    return inter_genes


def spTRS_DE(sptFile, Best_dict, cell_type1, cell_type2):
    # use sptinfo to load decontamination methods
    sptinfo = sptInfo(sptFile)
    platform = sptinfo.configs['platform'][0]
    dct_mtds = sptinfo.get_spmatrix()
    adatas = list(map(lambda dct: Load_spt_to_AnnData(sptFile, h5data=dct, platform=platform), dct_mtds))
    adatas = dict(zip(dct_mtds, adatas))
    Best_item1 = Best_dict[cell_type1]
    Best_item2 = Best_dict[cell_type2]
    d1 = Best_item1[1]
    d2 = Best_item2[1]
    mtds1 = Best_item1[2].split('+')
    mtds2 = Best_item2[2].split('+')
    dct_mtd1, dcv_mtd1, clu_mtd1 = mtds1[0], mtds1[1], mtds1[2]
    dct_mtd2, dcv_mtd2, clu_mtd2 = mtds2[0], mtds2[1], mtds2[2]
    adata_sp = adatas[dct_mtd1]
    domain1 = adata_sp.obs[clu_mtd1] == d1
    adata_sp = adatas[dct_mtd2]
    domain2 = adata_sp.obs[clu_mtd2] == d2
    de = DE_analysis(adata_sp, domain1, domain2, cell_type1, cell_type2)
    return de


# generate differential expressed genes (adjusted p value < 10e-5, t-test_overestim_var; effect > 0.2, cohen's d)
def Gen_DEgenes(adata, groupby, groups, pval=10e-5, effect=0.2, only_positive=True):
    sc.tl.rank_genes_groups(adata, groupby=groupby, groups=groups)
    sp_d = cohensD_eff(adata, groupby, groups)
    DE_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    DE_genes = pd.DataFrame({"names": DE_genes.iloc[:, 0].apply(str.upper),
                             "scores": pd.DataFrame(adata.uns['rank_genes_groups']['scores']).iloc[:, 0],
                             "pvals": pd.DataFrame(adata.uns['rank_genes_groups']['pvals']).iloc[:, 0],
                             "pvals_adj": pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj']).iloc[:, 0],
                             "logfoldchanges": pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges']).iloc[:, 0]
                             })
    DE_genes.index = DE_genes['names']
    DE_genes = DE_genes[DE_genes['pvals_adj'] < pval]
    sp_d = sp_d[sp_d["effect"] > effect]
    for idx in DE_genes.index:
        if idx in sp_d.index:
            DE_genes.loc[idx, "effect"] = sp_d.loc[idx, 'effect']
    if only_positive:
        DE_genes = DE_genes[DE_genes['scores'] > 0]
    return DE_genes


def Gen_Background_genes(adata_sp: ad.AnnData, adata_sc1: ad.AnnData, adata_sc2=None):
    sp_genes = list(map(lambda x: str.upper(x), list(adata_sp.var_names)))
    sc1_genes = list(map(lambda x: str.upper(x), list(adata_sc1.var_names)))
    if adata_sc2 is not None:
        sc2_genes = pd.DataFrame(adata_sc2.var_names).apply(str.upper)
        sc_genes = set(sc1_genes).intersection(set(sc2_genes.iloc[:, 0]))
    else:
        sc_genes = set(sc1_genes)
    bg_genes = list(sc_genes.intersection(set(sp_genes)))
    return bg_genes


# cohen's d effective size
def cohensD_eff(adata, groupby, groups):
    groups = [type(adata.obs[groupby][0])(i) for i in groups]
    idx = adata.obs[groupby].isin(groups)
    adata.obs_names_make_unique()
    nidx = ~idx  # select the remaining index
    adata1 = adata[idx].copy()
    adata2 = adata[nidx].copy()
    x1 = np.array(adata1.X.todense())
    x2 = np.array(adata2.X.todense())
    nx1 = len(x1)
    nx2 = len(x2)
    dof = nx1 + nx2 - 2
    cohens_d = abs(np.mean(x1, axis=0) - np.mean(x2, axis=0)) / np.sqrt(
        ((nx1 - 1) * np.std(x1, axis=0, ddof=1) ** 2 + (nx2 - 1) * np.std(x2, axis=0, ddof=1) ** 2) / dof)
    cohens_ds = pd.DataFrame(cohens_d, index=[str.upper(i) for i in adata.var_names], columns=['effect'])
    return cohens_ds


# multimodal intersection analysis, using hyper-geometric cumulative distribution
def MIA_analysis(inter_hvg, total_genes, sc_hvg, sp_hvg):
    return -np.log10(hypergeom.sf(inter_hvg, total_genes, sc_hvg, sp_hvg))


def DE_analysis(adata: ad.AnnData, domain1, domain2, cell_type1, cell_type2):
    sc.pp.log1p(adata)
    adata.obs["domain"] = "domain0"
    df = adata.obs["domain"]
    df[domain1 == True] = "domain1"
    df[domain2 == True] = "domain2"
    adata.obs["domain"] = df
    sc.tl.rank_genes_groups(adata, groupby="domain", groups=['domain1'], reference='domain2')

    gene_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    data = pd.DataFrame({"names": gene_df.iloc[:, 0].apply(str.upper),
                         "scores": pd.DataFrame(adata.uns['rank_genes_groups']['scores']).iloc[:, 0],
                         "pvals": pd.DataFrame(adata.uns['rank_genes_groups']['pvals']).iloc[:, 0],
                         "padj": pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj']).iloc[:, 0],
                         "log2FoldChange": pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges']).iloc[:, 0], })
    data.loc[(data.log2FoldChange > 1) & (data.pvals < 0.05), 'type'] = 'up'
    data.loc[(data.log2FoldChange < -1) & (data.pvals < 0.05), 'type'] = 'down'
    data.loc[(abs(data.log2FoldChange) <= 1) | (data.padj >= 0.05), 'type'] = 'nosig'
    data['-logpadj'] = -data.padj.apply(math.log10)
    data[['log2FoldChange', 'padj', 'type', '-logpadj']].head()
    colors = ["#01c5c4", "#686d76", "#ff414d"]
    sns.set_palette(sns.color_palette(colors))
    x_threshold = 1
    y_threshold = 5

    x_min = -6
    x_max = 6
    y_min = 0
    y_max = 80

    # draw volcano
    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    sns.scatterplot(x='log2FoldChange', y='-logpadj', data=data,
                    hue='type',  # 颜色映射
                    edgecolor=None,
                    s=6,  # 点大小
                    ax=ax)

    # names = ['PHLDA2', 'MEST', 'KRT5', 'KRT19', 'TGFA', 'EPCAM', 'CALML5', 'MGP', 'CCND1', 'KRT16',
    #          "PTGDS", 'CCL19', 'CD3E', 'CXCR4', "LTB", "ACAP1", "TSC22D3", "CD79A"]
    data.index = data['names']
    data_up = data.iloc[:20, :]
    data_down = data.iloc[-20:, :]
    names = data_up.index
    ax.scatter(data_up['log2FoldChange'], data_up['-logpadj'],
               edgecolors='black',
               facecolors='none',
               alpha=0.7,
               s=10,
               linewidths=0.5)
    for name in names:
        ax.annotate(name,  # context
                    xy=(data_up['log2FoldChange'][name], data_up['-logpadj'][name]),  # location of annotation
                    xytext=(2, 1),
                    textcoords='offset points',
                    fontsize=4)

    names = data_down.index
    ax.scatter(data_down['log2FoldChange'], data_down['-logpadj'],
               edgecolors='black',
               facecolors='none',
               alpha=0.7,
               s=10,
               linewidths=0.5)
    for name in names:
        ax.annotate(name,  # context
                    xy=(data_down['log2FoldChange'][name], data_down['-logpadj'][name]),  # location of annotation
                    xytext=(-15, 1),
                    textcoords='offset points',
                    fontsize=4)
    # labels
    ax.set(xlim=(x_min, x_max), ylim=(y_min, y_max), title=cell_type1 + " vs. " + cell_type2)
    ax.set_xlabel("log2FC")
    ax.set_ylabel("-log10(padj)")
    ax.vlines(-x_threshold, y_min, y_max, color='dimgrey', linestyles='dashed', linewidth=1)
    ax.vlines(x_threshold, y_min, y_max, color='dimgrey', linestyles='dashed', linewidth=1)
    ax.hlines(y_threshold, x_min, x_max, color='dimgrey', linestyles='dashed', linewidth=1)
    # 移动图例位置
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    fig.savefig("CID4971_minor_1/" + cell_type1 + " vs. " + cell_type2 + ".eps", dpi=400)
    return data


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
