import pandas as pd

from utils.io import *
from utils.preprocess import *
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from vsl.palette import Set_palette
from vsl.boundary import boundary_extract, show_edge
from PIL import Image


# report the smdFile, including brief information, umap graphs, deconvolution composition
# PVL mapping right
# B-cell mapping wrong
def Despot_report(smdFile, count="matrix"):
    # report brief information
    adata = Load_smd_to_AnnData(smdFile, count)
    figd, axsd = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
    sc.pl.spatial(adata, color="ground_truth", show=False, ax=axsd)
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
    sns.heatmap(df0.T, ax=axsd, cmap='Purples', linewidths=0.2)
    figd.savefig("spacexr1.pdf")

    # report deconvolution composition
    W0 = W.copy()
    W0['annotation'] = adata.obs['ground_truth'].copy()
    df1 = pd.DataFrame(W0.groupby(['annotation']).mean())
    fige, axse = plt.subplots(1, 1, figsize=(6, 3.5), constrained_layout=True)
    width = 0.5
    x = df1.index
    bottom_y = np.array([0] * len(df1))
    for cell, color in zip(df1.columns, mcolors.XKCD_COLORS):
        y = df1.loc[:, cell]
        axse.bar(x, y, width, bottom=bottom_y, label=W.columns, color=color)
        bottom_y = y + bottom_y
    axse.set_title('The cellular composition of each region')
    axse.set_xticks(x, x)
    axse.legend(W.columns, loc=2, bbox_to_anchor=(1.0, 1.0), borderaxespad=0.)


def Show_Origin_histology(smdFile, figname, hires=True, key=None):
    adata = Load_smd_to_AnnData(smdFile, hires=hires)
    fig, axs = plt.subplots(1, 1, figsize=(6, 3), constrained_layout=True)
    sc.pl.spatial(adata, img_key='hires', ax=axs, color=key)
    fig.savefig(figname, dpi=400)


def Ax_expr_spgenes(adata_sp, genes, ax, platform='ST', title=None, cmap='coolwarm'):
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams['font.size'] = 16
    expr = adata_sp.to_df().loc[:, genes]

    if platform == 'ST':
        ax.scatter(x=adata_sp.obs['array_row'], y=adata_sp.obs['array_col'], c=expr, cmap=cmap)
        ax.set_xticks(ticks=[])
        ax.set_yticks(ticks=[])
        ax.spines.top.set_visible(False)
        ax.spines.bottom.set_visible(False)
        ax.spines.left.set_visible(False)
        ax.spines.right.set_visible(False)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title("{0}".format(genes))
    elif platform == '10X_Visium':
        # show img
        img_key = list(adata_sp.uns['spatial'].keys())[0]
        sf = adata_sp.uns['spatial'][img_key]['scalefactors']['tissue_lowres_scalef']
        x,y = adata_sp.obs['image_col']*sf, adata_sp.obs['image_row']*sf
        xmin, xmax = int(np.min(x) - 5), int(np.max(x) + 5)
        ymin, ymax = int(np.min(y) - 5), int(np.max(y) + 5)
        x -= xmin
        y -= ymin
        ax.imshow(adata_sp.uns['spatial'][img_key]['images']['lowres'][ymin:ymax, xmin:xmax], alpha=0.8)
        ax.scatter(x=x, y=y, c=expr,
                       cmap=cmap, marker='.', s=10)
        ax.spines.top.set_visible(False)
        ax.spines.bottom.set_visible(False)
        ax.spines.left.set_visible(False)
        ax.spines.right.set_visible(False)
        ax.set_xticks(ticks=[])
        ax.set_yticks(ticks=[])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(title)


def Ax_expr_scgenes(adata_sc, genes, ax, title=None):
    sc.pl.umap(adata_sc, color=genes, show=False, ax=ax, cmap='coolwarm')
    ax.spines.top.set_visible(False)
    ax.spines.bottom.set_visible(False)
    ax.spines.left.set_visible(False)
    ax.spines.right.set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(title)


def Show_matched_genes(adata_sp, adata_sc, genes, figname, platform='ST'):
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False

    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    Ax_expr_spgenes(adata_sp, genes, ax=axs[0], platform=platform, title="")
    Ax_expr_scgenes(adata_sc, genes, ax=axs[1], title="")
    fig.suptitle("{0}".format(genes))
    fig.savefig(figname, dpi=400)


def Show_row_spgenes(adata_sp, genes, figname, platform='ST'):
    ngenes = len(genes)
    fig, axs = plt.subplots(1, ngenes, figsize=(3 * ngenes, 3), constrained_layout=True)
    Ax_expr_spgenes(adata_sp, genes, ax=axs[0], platform=platform, title="")
    fig.savefig(figname, dpi=400)


def Ax_prop_domains(adata, clu_mtd, dcv_mtd, domains, cell_type, ax, title, f1score, platform='ST', imgPath=None, cmap='cividis', ap=100):
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams['font.size'] = 16
    # select domain
    selection = np.zeros_like(adata.obs[clu_mtd])
    for d in domains:
        selection += (np.array(adata.obs[clu_mtd], dtype=int) == d)
    selection = np.array(selection, dtype=bool)

    # property
    prop = adata.obsm[dcv_mtd][cell_type]

    if platform == '10X_Visium':
        adata.obs[cell_type] = prop
        img_key = list(adata.uns['spatial'].keys())[0]
        sf = adata.uns['spatial'][img_key]['scalefactors']['tissue_lowres_scalef']

        # show img
        x,y = (adata.obs['image_col'])*sf, adata.obs['image_row']*sf
        xmin, xmax = int(np.min(x) - 5), int(np.max(x) + 5)
        ymin, ymax = int(np.min(y) - 5), int(np.max(y) + 5)
        x -= xmin
        y -= ymin
        pts = np.array([x[selection], y[selection]]).T
        alpha = ap / np.max(pts)
        edges, centers = boundary_extract(pts, alpha, err=10e-5)
        ax.imshow(adata.uns['spatial'][img_key]['images']['lowres'][ymin:ymax, xmin:xmax], alpha=0.8)
        ax.scatter(x=x, y=y, c=prop,
                       cmap=cmap, marker='.', s=10)
        show_edge(edges, ax, label='F1score: %.3f' % f1score, linewidth=1.5)
        ax.set_xticks(ticks=[])
        ax.set_yticks(ticks=[])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.spines.top.set_visible(False)
        ax.spines.bottom.set_visible(False)
        ax.spines.left.set_visible(False)
        ax.spines.right.set_visible(False)
        ax.legend(loc='lower right',
                      handletextpad=0.3,
                      borderpad=0.5,
                      borderaxespad=1.05,
                      columnspacing=0.7,
                      handlelength=0.7,
                      fontsize=12)
        ax.set_title(title)
    elif platform == 'ST':
        adata.obs[cell_type] = prop
        if imgPath is not None:
            img = Image.open(imgPath)
            sf = img.height / (np.max(adata.obs['image_col']) + 1.5 - np.min(adata.obs['image_col']))
            print(sf)
        else:
            sf = 1
        x, y = (adata.obs['image_col'])*sf, (-adata.obs['image_row'])*sf
        pts = np.array([x[selection], y[selection]]).T
        alpha = 15 / np.max(pts)
        edges, centers = boundary_extract(pts, alpha, err=10e-5)
        # ax.imshow(np.array(img), alpha=0.8)
        ax.scatter(x=x, y=y, c=prop,
                       cmap='cividis', marker='o', s=25)
        show_edge(edges, ax, label='F1score: %.3f' % f1score, linewidth=1.5)
        ax.set_xticks(ticks=[])
        ax.set_yticks(ticks=[])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.spines.top.set_visible(False)
        ax.spines.bottom.set_visible(False)
        ax.spines.left.set_visible(False)
        ax.spines.right.set_visible(False)
        ax.legend(loc='lower right',
                      handletextpad=0.3,
                      borderpad=0.5,
                      borderaxespad=1.05,
                      columnspacing=0.7,
                      handlelength=0.7,
                      fontsize=12)
        ax.set_title(title)


def Show_matched_domains(adata, clu_mtd, dcv_mtd, domains, cell_type, f1score, platform='ST'):
    from vsl.boundary import boundary_extract, show_edge

    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams['font.size'] = 14
    palette, cmap = Set_palette(len(np.unique(adata.obs[clu_mtd])))
    ncols = 1
    if len(np.unique(adata.obs[clu_mtd])) > 10:
        ncols = 2
    # set min label=0
    origin_label = np.array(adata.obs[clu_mtd], dtype=int)
    label = np.array(adata.obs[clu_mtd], dtype=int)
    min_label = np.min(label)
    label = label - min_label
    # labels = list(map(lambda x: palette[x], label))

    # set vpad of legends
    vpad = 0.1 + 0.01 * len(np.unique(label))

    # select domain
    selection = np.zeros_like(adata.obs[clu_mtd])
    for d in domains:
        selection += (np.array(adata.obs[clu_mtd], dtype=int) == d)
    selection = np.array(selection, dtype=bool)
    adata0 = adata[selection]
    # label0 = np.array(adata0.obs[clu_mtd], dtype=int)
    # label0 = label0 - min_label
    # labels0 = list(map(lambda x: palette[x], label0))
    if len(domains) == 1:
        dm = domains[0] - min_label
    else:
        dm = tuple(domains - min_label)

    # property
    prop = adata.obsm[dcv_mtd][cell_type]

    # plot style1
    if platform == 'ST':
        fig, axs = plt.subplots(1, 2, figsize=(7.5, 3.2), constrained_layout=True)
        handles = axs[0].scatter(x=adata.obs['array_row'], y=adata.obs['array_col'], c=label, label=label,
                                 cmap=cmap, marker='o', edgecolors='#606060', alpha=0.8)
        # axs[0].scatter(x=adata.obs['array_row'], y=adata.obs['array_col'], c=[], edgecolors='grey', alpha=0.8, marker='o', linewidths=0.8)
        # axs[0].scatter(x=adata0.obs['array_row'], y=adata0.obs['array_col'], c=labels0, marker='o', alpha=0.7)
        axs[0].scatter(x=adata0.obs['array_row'], y=adata0.obs['array_col'], c=[], marker='o', edgecolors='r', linewidths=0.8, alpha=0.8)
        axs[0].set_xticks(ticks=[])
        axs[0].set_yticks(ticks=[])
        axs[0].spines.top.set_visible(False)
        axs[0].spines.bottom.set_visible(False)
        axs[0].spines.left.set_visible(False)
        axs[0].spines.right.set_visible(False)
        axs[0].set_title("{0}: {1}".format(cell_type, dm))
        axs[0].legend(handles=handles.legend_elements()[0],
                      labels=list(np.unique(label)),
                      ncol=ncols,
                      frameon=False,
                      bbox_to_anchor=(-0.23, vpad),
                      handletextpad=0.1,
                      borderpad=0.5,
                      loc='lower left',
                      borderaxespad=1.05,
                      columnspacing=0.7)

        axs[1].scatter(x=adata.obs['array_row'], y=adata.obs['array_col'], c=prop, marker='o', cmap='viridis')
        # axs[1].scatter(x=adata.obs['array_row'], y=adata.obs['array_col'], c=[],
        #                label='others', edgecolors='grey', alpha=0.8, marker='o', linewidths=0.8)
        axs[1].scatter(x=adata0.obs['array_row'], y=adata0.obs['array_col'], c=[],
                       label='SCSP', marker='o', edgecolors='r', linewidths=0.8, alpha=0.7)
        axs[1].set_xticks(ticks=[])
        axs[1].set_yticks(ticks=[])
        axs[1].spines.top.set_visible(False)
        axs[1].spines.bottom.set_visible(False)
        axs[1].spines.left.set_visible(False)
        axs[1].spines.right.set_visible(False)
        axs[1].legend(loc='best',
                      handletextpad=0.05,
                      borderpad=0.5,
                      borderaxespad=1.05,
                      columnspacing=0.7)
        plt.colorbar(axs[1].get_children()[1], ax=axs[1], pad=0)
        axs[1].set_title("F1-score of SCSP: {0}".format(f1score))
        return fig

    elif platform == '10X_Visium':
        adata.obs[cell_type] = prop
        adata0 = adata[selection]
        pts = np.array(adata0.obs[['image_col', 'image_row']])
        alpha = 50 / np.max(pts)
        edges, centers = boundary_extract(pts, alpha, err=10e-5)
        fig, axs = plt.subplots(1, 2, figsize=(7, 3.5), constrained_layout=True)
        sc.pl.spatial(adata,
                      color=[clu_mtd],
                      palette=palette,
                      alpha=0.7,
                      ax=axs[0])
        show_edge(edges, axs[0], label='SCSP')
        axs[0].legend(frameon=False,
                      labels=list(np.unique(label)),
                      ncol=ncols,
                      bbox_to_anchor=(-0.35, vpad),
                      handletextpad=0.1,
                      borderpad=0.5,
                      loc='lower left',
                      borderaxespad=1.05,
                      columnspacing=0.7)
        axs[0].spines.top.set_visible(False)
        axs[0].spines.bottom.set_visible(False)
        axs[0].spines.left.set_visible(False)
        axs[0].spines.right.set_visible(False)
        axs[0].set_xlabel('')
        axs[0].set_ylabel('')
        axs[0].set_title("{0}: {1}".format(cell_type, dm))

        adata.obs[cell_type] = prop
        sc.pl.spatial(adata,
                      img_key= None,
                      color=[cell_type],
                      cmap='viridis',
                      ax=axs[1])
        show_edge(edges, axs[1], label='SCSP')
        axs[1].set_xticks(ticks=[])
        axs[1].set_yticks(ticks=[])
        axs[1].set_xlabel('')
        axs[1].set_ylabel('')
        axs[1].spines.top.set_visible(False)
        axs[1].spines.bottom.set_visible(False)
        axs[1].spines.left.set_visible(False)
        axs[1].spines.right.set_visible(False)
        axs[1].legend(loc='lower right',
                      handletextpad=0.3,
                      borderpad=0.5,
                      borderaxespad=1.05,
                      columnspacing=0.7,
                      handlelength=0.7)
        axs[1].set_title("F1-score of SCSP: %.3f" % f1score)
        return fig
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False


def Show_prop_domains(adata, clu_mtd, dcv_mtd, domains, cell_type, f1score, platform='ST'):
    fig, axs = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
    Ax_prop_domains(adata,
                    clu_mtd,
                    dcv_mtd,
                    domains,
                    cell_type=cell_type,
                    ax=axs,
                    title=cell_type,
                    f1score=f1score,
                    platform=platform)
    return fig


def Show_Spatial(smdFile, figname, key, h5data='matrix', do_cluster=True, platform="10X_Visium", title=''):
    adata = Load_smd_to_AnnData(smdFile, h5data, platform=platform)
    fig, axs = plt.subplots(1, 1, figsize=(4, 3), constrained_layout=True)
    if do_cluster:
        adata = Remove_mito_genes(adata)
        adata = Filter_genes(adata)
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata = Spatial_Cluster_Analysis(adata)
    sc.pl.spatial(adata, ax=axs, color=key, spot_size=1)
    axs.set_title(title)
    axs.set_xlabel('')
    axs.set_ylabel('')
    fig.savefig(figname, dpi=400)


def Show_Violin(smdFile, figname, key, h5data='matrix'):
    adata = Load_smd_to_AnnData(smdFile, h5data)
    adata.var_names_make_unique()
    fig, axs = plt.subplots(1, 1, figsize=(4, 3), constrained_layout=True)
    adata = Remove_mito_genes(adata)
    adata = Filter_genes(adata)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata = Spatial_Cluster_Analysis(adata)
    sc.pl.violin(adata, ax=axs, groupby="BayesSpace", keys=key)
    fig.savefig(figname, dpi=400)


def Show_Umap(smdFile, figname, h5data='matrix', is_spatial=True, key='ground_truth', title=None, do_cluster=True, figsize=(6,5)):
    if is_spatial:
        adata = Load_smd_to_AnnData(smdFile, h5data)
    else:
        adata = Load_smd_sc_to_AnnData(smdFile, h5data)
    if do_cluster:
        # adata = Remove_mito_genes(adata)
        # adata = Filter_genes(adata)
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata = Spatial_Cluster_Analysis(adata)
    fig, axs = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)
    sc.pl.umap(adata, ax=axs, color=key, show=False, size=12)
    axs.set_title(title)
    axs.legend(loc='lower left', frameon=False, bbox_to_anchor=(1,0), fontsize=12,
                handletextpad=0.3,
              borderpad=0.5,
              borderaxespad=1.05,
              columnspacing=0.7,
              handlelength=0.7)
    fig.savefig(figname, dpi=400)

def Show_colored_Umap(smdFile, figname, h5data='SPCS_mat'):
    adata = Load_smd_to_AnnData(smdFile, h5data)
    adata = Remove_mito_genes(adata)
    adata = Filter_genes(adata)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor='seurat', inplace=True)
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata)
    # use umap and leiden for clustering
    sc.tl.umap(adata, n_components=3)
    umap_color = adata.obsm['X_umap']
    norm_color = (umap_color - np.nanmin(umap_color)) / (np.nanmax(umap_color) - np.nanmin(umap_color))
    fig, axs = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)

    # axs[0].scatter(x=adata.obs['array_row'], y=adata.obs['array_col'], c=[], edgecolors='grey', alpha=0.8, marker='o', linewidths=0.8)
    # axs[0].scatter(x=adata0.obs['array_row'], y=adata0.obs['array_col'], c=labels0, marker='o', alpha=0.7)
    axs.scatter(x=adata.obs['array_col'], y=-adata.obs['array_row'], c=norm_color, marker='.')
    axs.set_xticks(ticks=[])
    axs.set_yticks(ticks=[])
    axs.spines.top.set_visible(False)
    axs.spines.bottom.set_visible(False)
    axs.spines.left.set_visible(False)
    axs.spines.right.set_visible(False)
    fig.savefig(figname, dpi=400)


def Show_Spatial_and_Umap(smdFile, h5data='matrix', key='BayesSpace'):
    adata = Load_smd_to_AnnData(smdFile, h5data)
    adata = Remove_mito_genes(adata)
    adata = Filter_genes(adata)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata = Spatial_Cluster_Analysis(adata)
    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    sc.pl.spatial(adata, ax=axs[0], color=key, legend_loc=None, spot_size=1)
    axs[0].set_title("Tissue Plot")
    sc.pl.umap(adata, ax=axs[1], color=key, show=False, size=18)
    axs[1].set_title("Umap Plot")
    fig.show()


def Show_Domain_and_Composition(smdFile, cell_type):
    info = smdInfo(smdFile)
    # arrangement of map chain: ["cell-type", "F1-score", "domain", "dct", "dcv", "clu"]
    map_chain = info.get_map_chains(cell_type)
    adata = Load_smd_to_AnnData(smdFile, map_chain['dct'])
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


def Show_bar_composition(smdFile, h5data, method="Cell2Location", subtype="leiden"):
    info = smdInfo(smdFile)

    dcv_mtds = info.get_dcv_methods(h5data)
    adata = Load_smd_to_AnnData(smdFile, h5data, loadDeconv=True)
    adata = adata[adata.obs['BayesSpace'] == 6]
    adata = Remove_mito_genes(adata)
    adata = Filter_genes(adata)
    adata = Remove_RPMTgenes(adata)
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    adata = Spatial_Cluster_Analysis(adata)
    import scipy.stats as stats
    adata = adata[adata.obs['SpaGCN'] != 0].copy()
    adata.obs['SpaGCN'][adata.obs['image_col'] < 4000] = 1
    adata.obs['domains'] = adata.obs['SpaGCN'].apply(lambda x: 'Invasive Carcinoma' if x == 1 else 'DCIS')
    sc.tl.rank_genes_groups(adata, groupby="domains",
                            method='wilcoxon',
                            corr_method='benjamini-hochberg',
                            groups=['Invasive Carcinoma'], reference='DCIS')
    gene_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    data = pd.DataFrame({"names": gene_df.iloc[:, 0].apply(str.upper),
                         "scores": pd.DataFrame(adata.uns['rank_genes_groups']['scores']).iloc[:, 0],
                         "pvals": pd.DataFrame(adata.uns['rank_genes_groups']['pvals']).iloc[:, 0],
                         "padj": pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj']).iloc[:, 0],
                         "log2FoldChange": pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges']).iloc[:, 0], })
    data = data[data.log2FoldChange < 6]
    data.loc[(data.log2FoldChange > 1) & (data.pvals < 0.05), 'type'] = 'up'
    data.loc[(data.log2FoldChange < -1) & (data.pvals < 0.05), 'type'] = 'down'
    data.loc[(abs(data.log2FoldChange) <= 1) | (data.pvals >= 0.05), 'type'] = 'nosig'
    data['-logpadj'] = -data.padj.apply(math.log10)
    data.loc[data['-logpadj'] < 0.5, 'type'] = 'nosig'
    data[['log2FoldChange', 'padj', 'type', '-logpadj']].head()
    colors = ["#01c5c4", "#686d76", "#ff414d"]
    sns.set_palette(sns.color_palette(colors))
    x_threshold = 0.5
    y_threshold = 0.5

    x_min = -5
    x_max = 5
    y_min = 0
    y_max = 6

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
                    fontsize=10)

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
                    fontsize=10)
    # labels
    ax.set(xlim=(x_min, x_max), ylim=(y_min, y_max))
    ax.set_xlabel("log2FC")
    ax.set_ylabel("-log10(padj)")
    ax.vlines(-x_threshold, y_min, y_max, color='dimgrey', linestyles='dashed', linewidth=1)
    ax.vlines(x_threshold, y_min, y_max, color='dimgrey', linestyles='dashed', linewidth=1)
    ax.hlines(y_threshold, x_min, x_max, color='dimgrey', linestyles='dashed', linewidth=1)
    # 移动图例位置
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)


    df1 = pd.DataFrame(adata.obsm[method])
    df = pd.concat([df1, adata.obs[subtype]], axis=1)
    invasive = df[df["SpaGCN"] == 1]
    dcis = df[df["SpaGCN"] == 3]
    for ct in df1.columns:
        print(ct, stats.ttest_ind(invasive[ct], dcis[ct]))

    # important genes
    # focus_genes = ["ENTPD1", 'PDCD1', "CD200", "IL7R", "HAVCR2", "CXCL13", "TIGIT", "TCF7", "ZNF683", "GZMK", "GZMA",
    #                "LAG3"]
    # focus_genes = ['ENTPD1', 'PDCD1', 'CD200']
    # df = pd.concat([adata.obs["SpaGCN"], pd.DataFrame(adata[:, focus_genes].X.toarray(),
    #                                                   columns=focus_genes, index=adata.obs_names)], axis=1)
    # fig, axs = plt.subplots(1, 1, figsize=(5.3, 4))
    # sns.boxplot(x="SpaGCN", y='ENTPD1', data=df, ax=axs, width=0.7, inner=None,
    #                palette=['#82a7a9', '#c86c86', '#b3999b'])
    # sns.swarmplot(x="SpaGCN", y='ENTPD1', data=df, ax=axs, color='#404040')
    # axs.set_xticks(ticks=[0, 1], labels=['invasive cancer', 'DCIS'], fontname="sans-serif", fontsize=12)
    # axs.set_yticks(ticks=list(range(2)), fontname="sans-serif", fontsize=14)
    # axs.set_ylabel(ylabel="Enrichment score", fontname="sans-serif", fontsize=14)
    # axs.set_xlabel(xlabel="", fontname="sans-serif", fontsize=14)
    # fig.show()

    def get_cell_prop(method, cell_type, subtype):
        df = pd.DataFrame(adata.obsm[method], copy=True)
        if method == 'Cell2Location':
            df = df.div(df.sum(axis=1), axis='rows')
        cell_prop = df[cell_type].copy()
        df0 = pd.DataFrame()
        df0['cell_prop'] = cell_prop
        df0['cell_type'] = cell_type
        df0['domains'] = adata.obs[subtype].apply(lambda x:'Invasive Carcinoma' if x==1 else 'DCIS')
        return df0

    def get_gene_expr(focus_genes, subtype):
        df = pd.DataFrame(adata[:, focus_genes].X.toarray(),
                          columns=focus_genes, index=adata.obs_names)
        df_all = pd.DataFrame()
        for gene in focus_genes:
            df0 = pd.DataFrame()
            df0['cell_prop'] = df[gene]
            df0['cell_type'] = gene
            df0['domains'] = adata.obs[subtype].apply(lambda x: 'Invasive Carcinoma' if x == 1 else 'DCIS')
            df_all = pd.concat([df_all, df0])
        return df_all

    def comp_sig(df, cell_type):
        sig_type=[]
        invasive = df[df['domains'] == 'Invasive Carcinoma']
        dcis = df[df['domains'] == 'DCIS']
        if type(cell_type) == str:
            cell_type = [cell_type]
        from tqdm import tqdm
        for x,ct in tqdm(enumerate(cell_type)):
            data1 = invasive['cell_prop'][invasive['cell_type'] == ct]
            data2 = dcis['cell_prop'][dcis['cell_type'] == ct]
            stat, p = stats.ranksums(data1, data2)
            if p <= 0.05:
                sig_type.append(ct)
        return sig_type

    def comp_plot(df, axs, cell_type, label=None):
        sig_type=[]
        from matplotlib.ticker import FuncFormatter
        sns.boxplot(x='cell_type', y='cell_prop', hue='domains', data=df, ax=axs, width=0.7,
                    palette=['#82a7a9', '#c86c86', '#b3999b'])
        invasive = df[df['domains'] == 'Invasive Carcinoma']
        dcis = df[df['domains'] == 'DCIS']
        if type(cell_type) == str:
            cell_type = [cell_type]

        ymax = df['cell_prop'].max()+0.025 * df['cell_prop'].max()
        for x,ct in enumerate(cell_type):
            data1 = invasive['cell_prop'][invasive['cell_type'] == ct]
            data2 = dcis['cell_prop'][dcis['cell_type'] == ct]
            stat, p = stats.ranksums(data1, data2)
            axs.plot([x-0.2,x+0.2], [ymax,ymax], color="black", linewidth=1)
            text_y = ymax
            if p <= 0.05:
                sig_type.append(ct)
            if p > 0.05:
                pval = 'ns'
                text_y = ymax+0.03 * ymax
            elif p > 0.01:
                pval = '*'
            elif p > 0.001:
                pval = '**'
            elif p > 0.0001:
                pval = '***'
            else:
                pval = '****'
            axs.text(x, text_y, pval, verticalalignment='center', horizontalalignment='center')
        # sns.swarmplot(x='cell_type', y='cell_prop', hue='domains', data=df, ax=axs, color='#404040')
        # axs.set_yticks(ticks=np.arange(0.1, 0.25, 0.05), fontsize=10)
        if label is not None:
            axs.set_xticklabels(labels=label)
        else:
            axs.set_xticklabels(labels=cell_type)
        # axs.set_title("terminal Tex (p=0.0042**)", fontname="sans-serif", fontsize=14)
        axs.set_ylabel(ylabel="subtype proportion", fontname="sans-serif", fontsize=14)
        axs.set_xlabel(xlabel="", fontname="sans-serif", fontsize=14)
        axs.spines.top.set_visible(False)
        axs.spines.right.set_visible(False)
        axs.spines["left"].set_linewidth(1)
        axs.spines["bottom"].set_linewidth(1)
        def ks(x, pos):
            return '%.2f' % (x)
        axs.yaxis.set_major_formatter(FuncFormatter(ks))
        axs.tick_params(axis="y", which="major", length=3, width=1)
        axs.tick_params(axis="x", which="major", length=3, width=1)
        axs.get_legend().remove()
        return sig_type

    df_tex = pd.DataFrame()
    subtype = 'SpaGCN'

    # terminal Tex
    cell_type = 'CD200+ terminal Tex'
    df0 = get_cell_prop("Cell2Location", cell_type, subtype)
    df_tex = pd.concat([df_tex, df0])

    # GZMB+ Tex
    cell_type = 'GZMB+ Tex'
    df1 = get_cell_prop("Cell2Location", cell_type, subtype)
    df_tex = pd.concat([df_tex, df1])

    df_tn = pd.DataFrame()
    # TCF7+ Tn1
    cell_type = 'TCF7+ Tn1'
    df0 = get_cell_prop("Cell2Location", cell_type, subtype)
    df_tn = pd.concat([df_tn, df0])

    # ZNF683+ Tm
    cell_type = 'ZNF683+ Tm'
    df1 = get_cell_prop("Cell2Location", cell_type, subtype)
    df_tn = pd.concat([df_tn, df1])

    from matplotlib.gridspec import GridSpec
    # fig, axs = plt.subplots(2, 1, figsize=(4, 6), constrained_layout=True)
    fig = plt.figure(layout="constrained", figsize=(3.7, 6))
    gs = GridSpec(16, 1, figure=fig)

    axs = [fig.add_subplot(gs[:7, :]),fig.add_subplot(gs[7, :]),
           fig.add_subplot(gs[8:-1, :]), fig.add_subplot(gs[-1, :])]
    # Tex
    cell_types = ['CD200+ terminal Tex', 'GZMB+ Tex']
    labels = ['terminal Tex', 'GZMB+ Tex']
    comp_plot(df_tex, axs[0], cell_types, labels)
    axs[0].set_title("Exhausted T-cells",fontsize=14)
    axs[1].axis('off')
    # Tn
    cell_types = [ 'TCF7+ Tn1', 'ZNF683+ Tm']
    labels = ['TCF7+ Tn1', 'ZNF683+ Tm']
    comp_plot(df_tn, axs[2], cell_types, labels)
    elem, label = axs[2].get_legend_handles_labels()
    axs[2].set_title("non-Exhausted T-cells",fontsize=14)
    axs[3].axis('off')
    fig.legend(elem, label,loc='lower left',
               frameon=False,
               bbox_to_anchor=(0.1,-0.03),
               fontsize=14, ncol=2,
               handletextpad=0.3,
               borderpad=0.5,
               borderaxespad=1.05,
               columnspacing=0.7,
               handlelength=0.7)
    # comp_plot(df2, axs[2], 'GZMK+ Tex')
    fig.savefig("figures/Fig.5/subtype-diff.svg", dpi=300)

    fig, axs = plt.subplots(1, 1, figsize=(7, 3), constrained_layout=True)
    # focus_genes = ["ENTPD1", 'PDCD1', "CD200", "IL7R", "HAVCR2", "CXCL13", "TIGIT", "TCF7", "ZNF683", "GZMK", "GZMA",
    #                "LAG3"]
    # "IL7R", "CCL3"  sig. genes
    genes = ["IFI6", "IL7R", 'PDCD10', 'CCR5', 'CD52', 'LAG3']
    labels = genes
    df_all = get_gene_expr(genes, subtype)
    comp_plot(df_all, axs, genes, labels)
    axs.set_title('DEGs of T-cell domains')
    axs.set_ylabel('gene expression')
    fig.savefig("figures/Fig.5/subtype-DEGs2.svg", dpi=300)

    # GZMK+ Tex
    method = "StereoScope"
    df1 = pd.DataFrame(adata.obsm[method])
    df = pd.concat([df1, adata.obs[subtype]], axis=1)
    fig, axs = plt.subplots(1, 1, figsize=(5.3, 4))
    sns.violinplot(x="SpaGCN", y='GZMK+ Tex', data=df, ax=axs, width=0.7, inner=None,
                   palette=['#82a7a9', '#c86c86', '#b3999b'])
    sns.swarmplot(x="SpaGCN", y='GZMK+ Tex', data=df, ax=axs, color='#404040')
    axs.set_xticks(ticks=[0, 1], labels=['invasive cancer', 'DCIS'], fontname="sans-serif", fontsize=12)
    axs.set_title("GZMK+ Tex (p<.0001)", fontname="sans-serif", fontsize=14)
    axs.set_xlabel(xlabel="", fontname="sans-serif", fontsize=14)
    axs.set_ylabel(ylabel="cellular composition", fontname="sans-serif", fontsize=14)
    fig.savefig("subcluster/GZMK+ Tex.eps", dpi=300)

    # Naive T
    method = "Seurat"
    df1 = pd.DataFrame(adata.obsm[method])
    df = pd.concat([df1, adata.obs[subtype]], axis=1)
    fig, axs = plt.subplots(1, 1, figsize=(5.3, 4))
    sns.violinplot(x="SpaGCN", y='TCF7+ Tn1', data=df, ax=axs, width=0.7, inner=None,
                   palette=['#82a7a9', '#c86c86', '#b3999b'])
    sns.swarmplot(x="SpaGCN", y='TCF7+ Tn1', data=df, ax=axs, color='#404040')
    axs.set_xticks(ticks=[0, 1], labels=['invasive cancer', 'DCIS'], fontname="sans-serif", fontsize=12)
    axs.set_title("TCF7+ Tn1 (p<.0001)", fontname="sans-serif", fontsize=14)
    axs.set_ylabel(ylabel="cellular composition", fontname="sans-serif", fontsize=14)
    axs.set_xlabel(xlabel="", fontname="sans-serif", fontsize=14)
    fig.savefig("subcluster/TCF7+ Tn1.eps", dpi=300)

    # ZNF683+ Tm
    method = "SPOTlight"
    df1 = pd.DataFrame(adata.obsm[method])
    df = pd.concat([df1, adata.obs[subtype]], axis=1)
    fig, axs = plt.subplots(1, 1, figsize=(5.3, 4))
    sns.violinplot(x="SpaGCN", y='ZNF683+ Tm', data=df, ax=axs, width=0.7, inner=None,
                   palette=['#82a7a9', '#c86c86', '#b3999b'])
    sns.swarmplot(x="SpaGCN", y='ZNF683+ Tm', data=df, ax=axs, color='#404040')
    axs.set_xticks(ticks=[0, 1], labels=['invasive cancer', 'DCIS'], fontname="sans-serif", fontsize=12)
    axs.set_title("ZNF683+ Tm (p<.0001)", fontname="sans-serif", fontsize=14)
    axs.set_ylabel(ylabel="cellular composition", fontname="sans-serif", fontsize=14)
    axs.set_xlabel(xlabel="", fontname="sans-serif", fontsize=14)
    fig.savefig("subcluster/ZNF683+ Tm.eps", dpi=300)
    for ct in df1.columns:
        adata.obs[ct] = df1[ct]
        fig, axs = plt.subplots(1, 1, figsize=(4, 4))
        sc.pl.spatial(adata, color=ct, legend_loc=None, ax=axs)
        fig.savefig("subcluster/" + ct + ".eps", dpi=300)
    # Cell2Location, spacexr, stereoScope, Seurat
    # CD200+ Terminal Tex, pval=0.004 **, 0.047 *, 0.452, nan,
    # GZMB+ Tex, pval=0.019 *, 0.045 *, 0.109, 0.005 **,
    # TCF7+ Tn1, pval=0.003 **, 0.019 *, 0.0008 ***, 1.102e-14 ****
    # ZNF683+ Tm, pval=4.633e-08(SPOTlight), 5.854e-09 ****

    df = pd.DataFrame(df.groupby(subtype).mean())
    fig, axs = plt.subplots(1, 1, figsize=(7, 3.5), constrained_layout=True)
    width = 0.5
    x = df1.index
    bottom_y = np.array([0] * len(df1))
    for cell, color in zip(df1.columns, mcolors.XKCD_COLORS):
        y = df1.loc[:, cell]
        axs.bar(x, y, width, bottom=bottom_y, label=df1.columns, color=color)
        bottom_y = y + bottom_y
    axs.set_title('The cellular composition of each region')
    axs.set_xticks(x, x)
    fig.show()


def Show_intersection_genes(adata_sp, adata_sc1, gene, folder, adata_sc2=None, cell_type=""):
    fig, axs = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
    adata_sp.var_names = pd.Series(adata_sp.var_names).apply(str.upper)
    sc.pl.spatial(adata_sp, ax=axs, color=gene, show=False, spot_size=1)
    figname = folder + '/' + gene + "_" + cell_type.replace('/', '|') + ".eps"
    fig.savefig(figname, dpi=400)

    fig, axs = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
    adata_sc1.var_names = pd.Series(adata_sc1.var_names).apply(str.upper)
    sc.pl.umap(adata_sc1, ax=axs, color=gene, show=False)
    figname = folder + '/' + gene + "_" + cell_type.replace('/', '|') + "_sc1.eps"
    fig.savefig(figname, dpi=400)

    if adata_sc2 is not None:
        fig, axs = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata_sc2 = Spatial_Cluster_Analysis(adata_sc2)
        adata_sc2.var_names = pd.Series(adata_sc2.var_names).apply(str.upper)
        sc.pl.umap(adata_sc2, ax=axs, color=gene, size=2, show=False)
        figname = folder + '/' + gene + "_" + cell_type.replace('/', '|') + "_sc2.eps"
        fig.savefig(figname, dpi=400)


def Show_Heatmap(smdFile, figname, h5data='matrix', title=None):
    adata = Load_smd_to_AnnData(smdFile, h5data)

    import seaborn as sns
    import scanpy as sc
    sc.pp.log1p(adata)
    X = pd.DataFrame(adata.X.todense(), index=adata.obs_names, columns=adata.var_names)
    X = X.iloc[0:200, 0:200]
    fig, axs = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
    sns.heatmap(X, cmap="viridis", ax=axs)
    axs.set_title(title)
    fig.savefig(figname, dpi=400)


def Show_marker_heatmap(smdFile, figname, h5data='scRNA_seq', title=None,is_spatial=False, ):
    if is_spatial:
        adata = Load_smd_to_AnnData(smdFile, h5data)
    else:
        adata = Load_smd_sc_to_AnnData(smdFile, h5data)
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata = Spatial_Cluster_Analysis(adata)
        sc.tl.rank_genes_groups(adata, 'annotation', method='wilcoxon')
        sc.tl.dendrogram(adata, 'annotation')
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show_gene_labels=True)


# downsampling comparison. this function draws a figure to compare cell expression patterns
# between downsampling and VAE.
def Show_DS_comparison(smdFile, figname,figsize=(3.5, 4.5)):
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams['font.size'] = 14
    from utils.stereoScope import Easy_Sample
    from clu.Cluster_leiden import Spatial_Cluster_Analysis
    from sklearn import metrics
    scdata0 = Load_smd_sc_to_AnnData(smdFile, 'sc-ref-es')
    scdata1 = Load_smd_sc_to_AnnData(smdFile, 'sc-ref')
    scdata0 = Spatial_Cluster_Analysis(scdata0)
    scdata1 = Spatial_Cluster_Analysis(scdata1)
    ARI0 = metrics.adjusted_rand_score(scdata0.obs['clusters'], scdata0.obs['annotation'])
    ARI1 = metrics.adjusted_rand_score(scdata1.obs['clusters'], scdata1.obs['annotation'])
    fig, axs = plt.subplots(2, 1, figsize=figsize, constrained_layout=True)
    sc.pl.umap(scdata0, color='clusters', palette=sc.pl.palettes.default_20, size=14, ax=axs[0], show=False)
    axs[0].set_title("DS-derived ARI={:.2f}".format(ARI0))
    sc.pl.umap(scdata1, color='clusters', palette=sc.pl.palettes.default_20, size=14, ax=axs[1], show=False)
    axs[1].set_title("VAE-derived ARI={:.2f}".format(ARI1))
    fig.savefig(figname, dpi=400)
