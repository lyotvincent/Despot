from SEDR.src.graph_func import graph_construction
from SEDR.src.utils_func import adata_preprocess
from SEDR.src.SEDR_train import SEDR_Train
import torch
import argparse
from utils.common import *
import pip
import pip._internal
import subprocess
import sys


def sedr_run(adata, n_clusters=20):
    def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
        '''
            arg1(adata)[AnnData matrix]
            arg2(fixed_clus_count)[int]

            return:
                resolution[int]
        '''
        for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
            sc.tl.leiden(adata, random_state=0, resolution=res)
            count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
            if count_unique_leiden == fixed_clus_count:
                break
        return res

    device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
    print('===== Using device: ' + device)
    # ################ Parameter setting
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', type=int, default=10, help='parameter k in spatial graph')
    parser.add_argument('--knn_distanceType', type=str, default='euclidean',
                        help='graph distance type: euclidean/cosine/correlation')
    parser.add_argument('--epochs', type=int, default=1000, help='Number of epochs to train.')
    parser.add_argument('--cell_feat_dim', type=int, default=200, help='Dim of PCA')
    parser.add_argument('--feat_hidden1', type=int, default=100, help='Dim of DNN hidden 1-layer.')
    parser.add_argument('--feat_hidden2', type=int, default=20, help='Dim of DNN hidden 2-layer.')
    parser.add_argument('--gcn_hidden1', type=int, default=32, help='Dim of GCN hidden 1-layer.')
    parser.add_argument('--gcn_hidden2', type=int, default=8, help='Dim of GCN hidden 2-layer.')
    parser.add_argument('--p_drop', type=float, default=0.02, help='Dropout rate.')
    parser.add_argument('--using_dec', type=bool, default=True, help='Using DEC loss.')
    parser.add_argument('--using_mask', type=bool, default=False, help='Using mask for multi-dataset.')
    parser.add_argument('--feat_w', type=float, default=10, help='Weight of DNN loss.')
    parser.add_argument('--gcn_w', type=float, default=0.1, help='Weight of GCN loss.')
    parser.add_argument('--dec_kl_w', type=float, default=10, help='Weight of DEC loss.')
    parser.add_argument('--gcn_lr', type=float, default=0.01, help='Initial GNN learning rate.')
    parser.add_argument('--gcn_decay', type=float, default=0.01, help='Initial decay rate.')
    parser.add_argument('--dec_cluster_n', type=int, default=10, help='DEC cluster number.')
    parser.add_argument('--dec_interval', type=int, default=20, help='DEC interval nnumber.')
    parser.add_argument('--dec_tol', type=float, default=0.00, help='DEC tol.')
    # ______________ Eval clustering Setting _________
    parser.add_argument('--eval_resolution', type=int, default=1, help='Eval cluster number.')
    parser.add_argument('--eval_graph_n', type=int, default=20, help='Eval graph kN tol.')
    parser.add_argument('--device', type=str, default=device, help='Using GPU or not.')
    params = parser.parse_args(args=[])
    adata.var_names_make_unique()
    adata_X = adata_preprocess(adata, min_cells=5, pca_n_comps=params.cell_feat_dim)
    graph_dict = graph_construction(adata.obsm['spatial'], adata.shape[0], params)
    params.cell_num = adata.shape[0]
    print('==== Graph Construction Finished')
    # ################## Model training
    sedr_net = SEDR_Train(adata_X, graph_dict, params)
    if params.using_dec:
        sedr_net.train_with_dec()
    else:
        sedr_net.train_without_dec()
    sedr_feat, _, _, _ = sedr_net.process()
    adata_sedr = ad.AnnData(sedr_feat)
    sc.pp.neighbors(adata_sedr, n_neighbors=params.eval_graph_n)
    sc.tl.umap(adata_sedr)
    eval_resolution = res_search_fixed_clus(adata_sedr, n_clusters)
    sc.tl.leiden(adata_sedr, key_added='clusters', resolution=eval_resolution)
    adata_sedr.obs_names = adata.obs_names
    adata.obs['clusters'] = adata_sedr.obs['clusters']

    return adata