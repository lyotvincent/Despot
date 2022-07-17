from utils.common import *


# the params are saved in a .json config
def Load_json_configs(path: str, encoding: str = 'utf8'):
    with open(path, 'r', encoding=encoding) as fp:
        json_data = json.load(fp)
    return json_data


def Load_spt_to_AnnData(sptFile: str,
                        count: str = "matrix",
                        hires: bool = False,
                        dtype='float32'):
    # use h5py to load origin info
    h5_obj = h5.File(sptFile, 'a')
    h5mat = h5_obj['matrix']
    h5img = h5_obj["sptimages"]
    features = h5mat['features']
    coord = h5img["coordinates"]
    sf = h5img["scalefactors"]

    # select the active data, default 'matrix'
    h5dat = h5_obj[count]

    # create X in AnnData, including data, indices, indptr, shape
    data = h5dat['data']
    indices = h5dat['indices']
    indptr = h5dat['indptr']
    shape = h5dat['shape']
    X = csc_matrix((data, indices, indptr), shape=shape, dtype=dtype).T

    # create obs and obsm
    obs = pd.DataFrame({'in_tissue': np.array(coord["tissue"][:], dtype='int32'),
                        'array_row': np.array(coord["row"][:], dtype='int32'),
                        'array_col': np.array(coord["col"][:], dtype='int32'),
                        'image_row': np.array(coord["imagerow"][:], dtype='int32'),
                        'image_col': np.array(coord["imagecol"][:], dtype='int32')},
                       index=np.array(coord["index"][:], dtype='str'))
    obs = obs.sort_index()
    obs_names = np.array(h5dat['barcodes'][:], dtype='str')
    obs = obs.loc[obs_names, :]

    # load idents
    idents = h5dat['idents']
    for ident in idents.keys():
        # change idents' datatype if they are numeric
        idf = np.array(idents[ident][:], dtype='str')
        if str.isnumeric(idf[0]):
            idf = np.array(idents[ident][:], dtype=int)
        # change idents to category
        idf = pd.Series(idf, index=np.array(h5mat['barcodes'][:], dtype='str'), dtype='category')
        obs[ident] = idf.loc[obs_names]
    obsm = np.array([obs['image_col'], obs['image_row']])
    obsm = obsm.T

    # create var
    if count != "matrix":  # which means Decont data exists
        h5fea = np.array(h5dat['features']['name'][:], dtype='str')
        var = pd.DataFrame({'gene_name': h5fea}, index=h5fea)
    else:
        # var = pd.DataFrame({'gene_ids': np.array(features["id"][:], dtype='str'),
        #                     'feature_types': np.array(features['feature_type'][:], dtype='str'),
        #                     'genome': np.array(features['genome'][:], dtype='str'),
        #                     'gene_name': np.array(features['name'][:], dtype='str')})
        # var.index = var['gene_ids']
        h5fea = np.array(h5dat['features']['name'][:], dtype='str')
        var = pd.DataFrame({'gene_name': h5fea}, index=h5fea)
    name = str(h5_obj["name"][0], encoding='utf8')

    adata = ad.AnnData(X, obs=obs, var=var)
    adata.obs_names = obs_names
    adata.obsm['spatial'] = obsm

    # create Images
    adata.uns['spatial'] = {name: {
        'images': {'lowres': h5img['lowres'][:].T},
        'scalefactors': {'spot_diameter_fullres': sf['spot'][0],
                         'tissue_hires_scalef': sf['hires'][0],
                         'fiducial_diameter_fullres': sf['fiducial'][0],
                         'tissue_lowres_scalef': sf['lowres'][0], },
        'metadata': {}
    }}
    if hires:
        adata.uns['spatial'][name]['images']['hires'] = h5img['hires'][:].T

    # create svgs
    if 'is_HVG' not in h5dat['features']:
        h5_obj.create_group(count + 'features/is_HVG')
    else:
        h5svg = h5dat['features']['is_HVG']
        HVGS = {}
        if len(h5svg) > 0:
            for svg in h5svg.keys():
                HVG = pd.DataFrame()
                if svg == "BayesSpace" or svg == "leiden":
                    pass
                else:
                    for elem in h5svg[svg].keys():
                        HVG = pd.concat([HVG, pd.Series(h5svg[svg][elem], name=elem)], axis=1)
                HVGS[svg] = HVG
        adata.uns["HVGs"] = HVGS

    # create deconv results
    if 'deconv' not in h5dat:
        h5_obj.create_group(count + '/deconv')
    else:
        h5dcv = h5dat['deconv']
        if len(h5dcv) > 0:
            for dcv in h5dcv.keys():
                shape = h5dcv[dcv]['shape']
                weights = np.array(h5dcv[dcv]['weights']).reshape(shape[1], shape[0]).T
                barcodes = np.array(h5dcv[dcv]['barcodes'], dtype='str')
                cell_type = np.array(h5dcv[dcv]['cell_type'], dtype='str')
                w = pd.DataFrame(weights, index=barcodes, columns=cell_type)
                adata.obsm[dcv] = w
    h5_obj.close()
    return adata


def Load_spt_sc_to_AnnData(sptFile):
    h5_obj = h5.File(sptFile, 'r')
    # the scRNA-seq Data is saved in 'scRNA_seq', using Sparse Matrix
    scRNA_seq = h5_obj['scRNA_seq']
    data = scRNA_seq['data']
    indices = scRNA_seq['indices']
    indptr = scRNA_seq['indptr']
    shape = scRNA_seq['shape']
    idents = scRNA_seq['idents']

    # create X in AnnData
    X = csc_matrix((data, indices, indptr), shape=shape, dtype='int32').T
    # create obs in AnnData
    obs = pd.DataFrame({'annotation': np.array(idents['annotation'][:], dtype='str')},
                       index=np.array(scRNA_seq['barcodes'][:], dtype='str'))
    # create var in AnnData
    var = pd.DataFrame(index=np.array(scRNA_seq['features']['name'][:], dtype='str'))
    # create AnnData
    adata = ad.AnnData(X, obs=obs, var=var)
    adata.uns['HVGs'] = np.array(scRNA_seq['features']['HVGs'][:], dtype='str')
    return adata


def Load_weight_to_AnnData(adata: ad.AnnData, Wfile):
    W = pd.read_csv(Wfile, sep='\t', index_col=0)
    adata.obsm['StereoScope'] = W.copy()
    # obs = adata.obs
    # adata.obs = pd.concat([obs, W], axis=1)
    # W0 = W.copy()
    # W0['annotation'] = adata.obs['BayesSpace'].copy()
    # df = W0.groupby(['annotation']).sum()
    # from sklearn.preprocessing import MinMaxScaler
    # scaler = MinMaxScaler()
    # tran = scaler.fit_transform(df)
    # df0 = pd.DataFrame(tran, columns=df.columns, index=df.index)
    #
    # figa, axsa = plt.subplots(4, 4, figsize=(12, 12), constrained_layout=True)
    # for i in range(4):
    #     for j in range(4):
    #         if 4 * i + j >= len(W.columns):
    #             break
    #         sc.pl.spatial(adata, img_key="lowres", color=W.columns[4 * i + j], size=1, ax=axsa[i, j], show=False)
    # sc.pl.spatial(adata, img_key="lowres", color="BayesSpace", size=1, ax=axsa[3, 2])
    return adata


def Load_spt_to_Benchmark(sptFile, h5data, mode: str = "cluster"):
    if mode == "cluster":
        with h5.File(sptFile, 'r') as f:
            bmc = f[h5data]['benchmark']['cluster']
            shape = bmc['shape']
            value = np.array(bmc['value']).reshape(shape[1], shape[0]).T

            methods = np.array(bmc['methods'], dtype="str")
            targets = np.array(bmc['targets'], dtype="str")
            bm = pd.DataFrame(value, index=methods, columns=targets)
        return bm
    elif mode == "estimate":
        bm = {}
        with h5.File(sptFile, 'r') as f:
            bme = f[h5data]['benchmark']['estimate']
            for method in bme.keys():
                bme_met = bme[method]
                # name = np.array(bme_met['name'], dtype='str')
                moransI = np.array(bme_met['moransI'], dtype="float32")
                pval = np.array(bme_met['pval'], dtype="float32")
                bm[method] = pd.DataFrame({"moransI": moransI, "pval": pval})
        return bm


def Save_tsv_from_scData(out_dir, scdata):
    if not osp.exists(out_dir):
        os.mkdir(out_dir)
    mat = scdata.X.todense()
    cnt = pd.DataFrame(mat, index=scdata.obs_names,
                       columns=scdata.var_names)
    opth_cnt = osp.join(out_dir, 'cnt_data.tsv')
    cnt.to_csv(opth_cnt, sep='\t', header=True, index=True, index_label='cell')

    mta = pd.DataFrame({'bio_celltype': scdata.obs['annotation']}, index=scdata.obs_names)
    opth_mta = osp.join(out_dir, 'mta_data.tsv')
    mta.to_csv(opth_mta, sep='\t', header=True, index=True, index_label='cell')


def Save_tsv_from_spData(out_dir, spdata):
    if not osp.exists(out_dir):
        os.mkdir(out_dir)
    mat = spdata.X.todense()
    spdata.var_names = spdata.var['gene_name']
    spdata.var_names_make_unique()
    cnt = pd.DataFrame(mat, index=spdata.obs_names,
                       columns=spdata.var_names)
    opth_cnt = osp.join(out_dir, 'spt_data.tsv')
    cnt.to_csv(opth_cnt, sep='\t', header=True, index=True, index_label='cell')


def Load_spt_to_stData(sptFile, count="matrix", hires=False, dataPath=None):
    import stlearn as st
    adata = Load_spt_to_AnnData(sptFile, count, hires, dataPath)
    adata = st.convert_scanpy(adata)
    return adata


def Save_spt_from_leiden(sptFile, adata, h5data='matrix'):
    h5_obj = h5.File(sptFile, 'r')
    fullgenes = np.array(h5_obj[h5data]['features']['name'][:], dtype='str')
    h5_obj.close()
    fg = pd.Series(np.zeros(len(fullgenes), dtype='bool'), index=fullgenes)
    with h5.File(sptFile, 'a') as f:
        # save idents
        ident = list(np.array(adata.obs['clusters'], dtype='int32'))
        if h5data + '/idents/leiden' not in f:
            f.create_dataset(h5data + '/idents/leiden', data=ident, dtype='int32')
        else:
            del f[h5data + '/idents/leiden']
            f.create_dataset(h5data + '/idents/leiden', data=ident, dtype='int32')

        # # save highly variable genes
        # hvg = adata.var['highly_variable']
        # hvgt = hvg[hvg is True]
        # for gene in hvgt.index:
        #     fg.loc[gene] = True
        # f.create_dataset(h5data + '/features/is_HVG/leiden', data=list(fg), dtype='bool')


def Save_spt_from_Squidpy_clu(sptFile, adata, h5data='matrix'):
    # h5_obj = h5.File(sptFile, 'r')
    # fullgenes = np.array(h5_obj[h5data]['features']['name'][:], dtype='str')
    # h5_obj.close()
    # fg = pd.Series(np.zeros(len(fullgenes), dtype='bool'), index=fullgenes)
    with h5.File(sptFile, 'a') as f:
        # save Squidpy feature idents, which clustered by images
        ident = list(np.array(adata.obs['features_cluster'], dtype='int32'))
        if h5data + '/idents/Squidpy_img' not in f:
            f.create_dataset(h5data + '/idents/Squidpy', data=ident, dtype='int32')
        else:
            del f[h5data + '/idents/Squidpy_img']
            f.create_dataset(h5data + '/idents/Squidpy', data=ident, dtype='int32')


def Save_spt_from_SpaGCN_clu(sptFile, adata, h5data='matrix'):
    with h5.File(sptFile, 'a') as f:
        # save idents
        ident = list(np.array(adata.obs['clusters'], dtype='int32'))
        if h5data + '/idents/SpaGCN' not in f:
            f.create_dataset(h5data + '/idents/SpaGCN', data=ident, dtype='int32')
        else:
            del f[h5data + '/idents/SpaGCN']
            f.create_dataset(h5data + '/idents/SpaGCN', data=ident, dtype='int32')
    print("Clustering with `SpaGCN` finished, idents saved in /" + h5data + '/idents/SpaGCN')


def Save_spt_from_stlearn(sptFile, stdata, h5data='matrix'):
    with h5.File(sptFile, 'a') as f:
        # save idents
        ident = list(np.array(stdata.obs['louvain'], dtype='int32'))
        f.create_dataset(h5data + '/idents/stlearn', data=ident, dtype='int32')
    print("Clustering with `stlearn` finished, idents saved in /" + h5data + '/idents/stlearn')


def Save_spt_from_SpatialDE(sptFile, adata, h5data='matrix'):
    # FSV: Fraction of variance explained by spatial variation 由空间变异解释的方差比例
    # P value: SpatialDE-logP Value,  表示空间变异的显著性
    with h5.File(sptFile, 'a') as f:
        # create group
        if h5data + '/fearures/is_HVG/SpatialDE' not in f:
            f.create_group(h5data + '/fearures/is_HVG/SpatialDE')
        else:
            del f[h5data + '/features/is_HVG/SpatialDE/FSV']
            del f[h5data + '/features/is_HVG/SpatialDE/pval']
            del f[h5data + '/features/is_HVG/SpatialDE/qval']
            del f[h5data + '/features/is_HVG/SpatialDE/name']
            del f[h5data + '/fearures/is_HVG/SpatialDE']
            f.create_group(h5data + '/fearures/is_HVG/SpatialDE')

        # save Gene Name
        name = list(np.array(adata.var.index, dtype='S'))
        f.create_dataset(h5data + '/features/is_HVG/SpatialDE/name', data=name)
        # save FSV
        FSV = list(np.array(adata.var['FSV'], dtype='float32'))
        f.create_dataset(h5data + '/features/is_HVG/SpatialDE/FSV', data=FSV, dtype='float32')
        # save P-Value
        pval = list(np.array(adata.var['pval'], dtype='float32'))
        f.create_dataset(h5data + '/features/is_HVG/SpatialDE/pval', data=pval, dtype='float32')
        # save Q-Value
        qval = list(np.array(adata.var['qval'], dtype='float32'))
        f.create_dataset(h5data + '/features/is_HVG/SpatialDE/qval', data=qval, dtype='float32')
        # FDR
    print("Estimation with `SpatialDE` finished, HVGs saved in /" + h5data + '/features/is_HVG/SpatialDE')


def Save_spt_from_SpaGCN_est(sptFile, adata, h5data='matrix'):
    with h5.File(sptFile, 'a') as f:
        # create group
        if h5data + '/fearures/is_HVG/SpaGCN' not in f:
            f.create_group(h5data + '/fearures/is_HVG/SpaGCN')
        else:
            del f[h5data + '/fearures/is_HVG/SpaGCN']
            del f[h5data + '/features/is_HVG/SpaGCN/name']
            del f[h5data + '/features/is_HVG/SpaGCN/pval']
            del f[h5data + '/features/is_HVG/SpaGCN/cluster']
            del f[h5data + '/features/is_HVG/SpaGCN/fold_change']
            f.create_group(h5data + '/fearures/is_HVG/SpaGCN')

        # save Gene Name
        name = list(np.array(adata.uns['filtered_info']['genes'], dtype='S'))
        f.create_dataset(h5data + '/features/is_HVG/SpaGCN/name', data=name)

        # save P-Value
        pval = list(np.array(adata.uns['filtered_info']['pvals_adj'], dtype='float32'))
        f.create_dataset(h5data + '/features/is_HVG/SpaGCN/pval', data=pval, dtype='float32')

        # save cluster
        cluster = list(np.array(adata.uns['filtered_info']['target_dmain'], dtype='int32'))
        f.create_dataset(h5data + '/features/is_HVG/SpaGCN/cluster', data=cluster, dtype='int32')

        # save fold change
        fc = list(np.array(adata.uns['filtered_info']['fold_change'], dtype='float32'))
        f.create_dataset(h5data + '/features/is_HVG/SpaGCN/fold_change', data=fc, dtype='float32')

    print("Estimation with `SpaGCN` finished, HVGs saved in /" + h5data + '/features/is_HVG/SpaGCN')


def Save_spt_from_StereoScope(sptFile, Wfile, h5data='matrix'):
    W = pd.read_csv(Wfile, sep='\t', index_col=0).T
    shape = [W.shape[1], W.shape[0]]
    weights = np.array(W).flatten()
    cell_type = list(W.index)
    barcodes = list(W.columns)

    with h5.File(sptFile, 'a') as f:
        if h5data + '/deconv/StereoScope' not in f:
            f.create_group(h5data + '/deconv/StereoScope')
        else:
            del f[h5data + '/deconv/StereoScope/weights']
            del f[h5data + '/deconv/StereoScope/shape']
            del f[h5data + '/deconv/StereoScope/cell_type']
            del f[h5data + '/deconv/StereoScope/barcodes']
            del f[h5data + '/deconv/StereoScope']
            f.create_group(h5data + '/deconv/StereoScope')
        f.create_dataset(h5data + '/deconv/StereoScope/weights', data=weights)
        f.create_dataset(h5data + '/deconv/StereoScope/cell_type', data=cell_type)
        f.create_dataset(h5data + '/deconv/StereoScope/barcodes', data=barcodes)
        f.create_dataset(h5data + '/deconv/StereoScope/shape', data=shape)
    print("Deconvolution with `StereoScope` finished, Weight saved in /" + h5data + '/deconv/StereoScope')


def Save_spt_from_Cell2Location(sptFile, adata, h5data='matrix'):
    W = adata.obsm['Cell2Location'].T
    shape = [W.shape[1], W.shape[0]]
    weights = np.array(W).flatten()
    cell_type = list(W.index)
    barcodes = list(W.columns)

    with h5.File(sptFile, 'a') as f:
        if h5data + '/deconv/Cell2Location' not in f:
            f.create_group(h5data + '/deconv/Cell2Location')
        else:
            del f[h5data + '/deconv/Cell2Location/weights']
            del f[h5data + '/deconv/Cell2Location/shape']
            del f[h5data + '/deconv/Cell2Location/cell_type']
            del f[h5data + '/deconv/Cell2Location/barcodes']
            del f[h5data + '/deconv/Cell2Location']
            f.create_group(h5data + '/deconv/Cell2Location')
        f.create_dataset(h5data + '/deconv/Cell2Location/weights', data=weights)
        f.create_dataset(h5data + '/deconv/Cell2Location/cell_type', data=cell_type)
        f.create_dataset(h5data + '/deconv/Cell2Location/barcodes', data=barcodes)
        f.create_dataset(h5data + '/deconv/Cell2Location/shape', data=shape)
    print("Deconvolution with `Cell2Location` finished, Weight saved in /" + h5data + '/deconv/Cell2Location')
