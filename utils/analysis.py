from utils.common import *


# 对聚类结果的可视化
def Cluster_Visualization(adata, type='all'):
    print(adata)
    if type=='umap':
        sc.pl.umap(adata, color='clusters', palette=sc.pl.palettes.default_20)
    elif type=='spatial':
        sc.pl.spatial(adata, img_key="lowres", color="clusters", size=1.5)
    else:
        fig, axs = plt.subplots(2, 2, figsize=(16, 16), constrained_layout=True)
        sc.pl.spatial(adata, img_key="lowres", color=['total_counts'], ax=axs[0, 0], show=False)
        sc.pl.spatial(adata, img_key="lowres", color=['n_genes'], ax=axs[0, 1], show=False)
        sc.pl.umap(adata, color='clusters', palette=sc.pl.palettes.default_20, ax=axs[1, 0], show=False)
        sc.pl.spatial(adata, img_key="lowres", color="clusters", size=1.5, ax=axs[1, 1], show=False)
        plt.show()

# 差异表达分析
def Differential_Expression_Analysis(adata):
    sc.tl.rank_genes_groups(adata, "clusters", inplace=True)
    sc.pl.rank_genes_groups_heatmap(adata, groups="4", groupby="clusters", show=False)
    plt.show()
    return adata


# 将得到的前四个高可变基因输出
def HDG_Visualization(adata, gene_list):
    fig, axs = plt.subplots(2, 2, figsize=(16, 16), constrained_layout=True)
    sc.pl.spatial(adata, img_key="hires", color=gene_list[0], show=False, ax=axs[0, 0])
    sc.pl.spatial(adata, img_key="hires", color=gene_list[1], show=False, ax=axs[0, 1])
    sc.pl.spatial(adata, img_key="hires", color=gene_list[2], show=False, ax=axs[1, 0])
    sc.pl.spatial(adata, img_key="hires", color=gene_list[3], show=False, ax=axs[1, 1])
    plt.show()
