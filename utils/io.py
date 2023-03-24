from utils.common import *

# summary of a sptFile
class sptInfo:
    def __init__(self, sptFile):
        self.sptFile = sptFile
        # add names
        self.name = None
        with h5.File(sptFile, 'r') as f:
            mat = list(f.keys())
            if 'name' in mat:
                self.name = bytes2str(f['name'][:])[0]

        # add origin matrix
        self.spmatrixs = {}
        if 'matrix' in mat:
            self.spmatrixs['matrix'] = {}
        # add decontaminated matrix
        for m in mat:
            if m[-4:] == '_mat':
                self.spmatrixs[m] = {}

        # add scRNA-seq matrix
        self.scmatrixs = []
        if 'scRNA_seq' in mat:
            self.scmatrixs.append('scRNA_seq')

        # for each spmatrix, add clu_methods and dcv_methods
        with h5.File(self.sptFile, 'r') as f:
            for keys in self.spmatrixs.keys():
                h5mat = f[keys]
                clu_mtds = list(h5mat['idents'].keys())
                dcv_mtds = list(h5mat['deconv'].keys())
                self.spmatrixs[keys]['cluster_methods'] = clu_mtds
                self.spmatrixs[keys]['deconvolution_methods'] = dcv_mtds

        # add map_chain if available
        self.map_chains = pd.DataFrame()
        if 'map_chains' in mat:
            with h5.File(self.sptFile, 'r') as f:
                for key in f['map_chains'].keys():
                    dat = f['map_chains'][key][:]
                    if key in ['cell-type', 'clu', 'dct', 'dcv']:
                        dat = bytes2str(dat)
                    elif key == 'domain':
                        dat = bytes2tuple(dat)
                    dat = pd.DataFrame(dat, columns=[key])
                    self.map_chains = pd.concat([self.map_chains, dat], axis=1)
            self.map_chains.index = self.map_chains['cell-type']

        # add configs
        self.configs = {}
        with h5.File(self.sptFile, 'r') as f:
            if "configs" in f:
                cfg = f["configs"]
                for k in cfg.keys():
                    item = np.array(cfg[k], dtype='str')
                    if item.shape == ():    # for only one item
                        item = item.item()
                    else :     # for multi items
                        item = list(item)
                    self.configs[k] = item


    def get_spmatrix(self):
        return list(self.spmatrixs.keys())

    def get_clu_methods(self, spmat: str):
        if spmat == 'matrix':
            return self.spmatrixs[spmat]['cluster_methods']
        else:
            if spmat[-4:] != '_mat':
                spmat = spmat + '_mat'
            if spmat in self.spmatrixs.keys():
                return self.spmatrixs[spmat]['cluster_methods']
            else:
                return None

    def get_dcv_methods(self, spmat: str):
        if spmat == 'matrix':
            return self.spmatrixs[spmat]['deconvolution_methods']
        else:
            if spmat[-4:] != '_mat':
                spmat = spmat + '_mat'
            if spmat in self.spmatrixs.keys():
                return self.spmatrixs[spmat]['deconvolution_methods']
            else:
                return None

    def get_map_chains(self, cell_type=None):
        if cell_type in self.map_chains.index:
            return self.map_chains.loc[cell_type, :]
        else:
            return self.map_chains


# save configs to sptFile
def Save_spt_from_configs(sptFile, cfg_file="params.json", items: dict = None, change_cfg=False):
    def get_item(items,k):
        item = items[k]
        if type(item) == str:
            if item == "":
                item = "empty"
            item = bytes(item, encoding='utf-8')
        elif type(item) == list and len(item) > 0:
            if type(item[0]) == str:
                item = str2bytes(item)
        return item

    with h5.File(sptFile, 'a') as f:
        if "configs" not in f:
            # save configs from params.json if no configs or needs to change configs
            f.create_group("configs")
            cfg = Load_json_configs(cfg_file)
            for k in cfg.keys():
                item = get_item(cfg, k)
                f.create_dataset(name="configs/" + k, data=item)
        else:
            # if item is none, means no changes in configs
            if items is None:
                if change_cfg:
                    items = Load_json_configs(cfg_file)
                else:
                    return 0
            for k in items.keys():
                if k in f["configs"]:
                    del f["configs/" + k]
                item = get_item(items, k)
                print(item)
                f.create_dataset(name="configs/" + k, data=item)
    return 0


def Spt_check_corrections(sptFile):
    def check_matrix(h5obj):
        basic_item = ['barcodes', 'features', 'data', 'indices', 'indptr', 'shape']
        for item in basic_item:
            if item not in h5obj:
                print("check matrix corrections failed in ", h5obj)
                return False
        return True

    def check_images(h5obj):
        basic_item = ['coordinates/row', 'coordinates/col']
        for item in basic_item:
            if item not in h5obj:
                print("check matrix corrections failed in ", h5obj)
                return False
        return True

    true_flag = True
    # check configs
    Save_spt_from_configs(sptFile)

    with h5.File(sptFile, 'r') as f:
        if 'sptimages' not in f:
            print("check matrix corrections failed in sptimages.")
            return False
        else:
            if not check_images(f['sptimages']):
                return False

        if 'matrix' not in f:
            print("check matrix corrections failed in matrix.")
            return False
        else:
            if not check_matrix(f['matrix']):
                return False
    return true_flag

# the params are saved in a .json config
def Load_json_configs(path: str, encoding: str = 'utf-8'):
    with open(path, 'r', encoding=encoding) as fp:
        json_data = json.load(fp)
    return json_data


# convert bytes to str in array
def bytes2str(arr):
    arr = np.array(list(map(lambda x: str(x, encoding='utf8'), arr)))
    return arr


# convert str to bytes in array
def str2bytes(arr):
    arr = np.array(list(map(lambda x: bytes(x, encoding='utf-8'), arr)))
    return arr


def bytes2tuple(arr):
    arr = np.array(list(map(lambda x: tuple(x), arr)), dtype=object)
    arr = np.array(list(map(lambda x: (0, ) if len(x) == 0 else x, arr)), dtype=object)
    return arr


def tuple2bytes(arr):
    arr = np.array(list(map(lambda x: bytes(x), arr)))
    return arr


def Load_spt_to_AnnData(sptFile: str,
                        h5data: str = "matrix",
                        hires: bool = False,
                        loadDeconv: bool = False,
                        loadImg: bool = True,
                        platform: str = '10X_Visium',
                        dtype='float32'):
    # use h5py to load origin info
    with h5.File(sptFile, 'r') as h5_obj:
        h5mat = h5_obj['matrix']
        h5img = h5_obj["sptimages"]
        coord = h5img["coordinates"]
        sf = h5img["scalefactors"]

        # select the active data, default 'matrix'
        datas = h5data.split('/')
        data0 = datas[0]
        h5dat = h5_obj[data0]

        # create X in AnnData, including data, indices, indptr, shape
        data = h5dat['data']
        indices = h5dat['indices']
        indptr = h5dat['indptr']
        shape = h5dat['shape']
        X = csc_matrix((data, indices, indptr), shape=shape, dtype=dtype).T

        # create obs and obsm
        if platform == "ST":    # rotate anticlockwise 90 degrees for ST
            obs = pd.DataFrame({'in_tissue': np.array(coord["tissue"][:], dtype='int32'),
                                'array_row': np.array(coord["row"][:], dtype='int32'),
                                'array_col': np.array(coord["col"][:], dtype='int32'),
                                'image_row': np.array(-coord["imagecol"][:], dtype='int32'),
                                'image_col': np.array(coord["imagerow"][:], dtype='int32')},
                               index=bytes2str(coord["index"][:]))
        else:
            obs = pd.DataFrame({'in_tissue': np.array(coord["tissue"][:], dtype='int32'),
                                'array_row': np.array(coord["row"][:], dtype='int32'),
                                'array_col': np.array(coord["col"][:], dtype='int32'),
                                'image_row': np.array(coord["imagerow"][:], dtype='int32'),
                                'image_col': np.array(coord["imagecol"][:], dtype='int32')},
                               index=bytes2str(coord["index"][:]))
        obs = obs.sort_index()
        obs_names = bytes2str(h5dat['barcodes'][:])
        obs = obs.loc[obs_names, :]

        # load idents
        idents = h5dat['idents']
        for ident in idents.keys():
            # change idents' datatype if they are numeric
            idf = np.array(idents[ident][:], dtype='str')
            if str.isnumeric(idf[0]):
                idf = np.array(idents[ident][:], dtype=int)
            # change idents to category
            idf = pd.Series(idf, index=bytes2str(h5mat['barcodes'][:]), dtype='category')
            obs[ident] = idf.loc[obs_names]
        obsm = np.array([obs['image_col'], obs['image_row']])
        obsm = obsm.T

        # create var
        if data0 != "matrix":  # which means Decont data exists
            h5fea = bytes2str(h5dat['features']['name'][:])
            var = pd.DataFrame({'gene_name': h5fea}, index=h5fea)
        else:
            # var = pd.DataFrame({'gene_ids': np.array(features["id"][:], dtype='str'),
            #                     'feature_types': np.array(features['feature_type'][:], dtype='str'),
            #                     'genome': np.array(features['genome'][:], dtype='str'),
            #                     'gene_name': np.array(features['name'][:], dtype='str')})
            # var.index = var['gene_ids']
            h5fea = bytes2str(h5dat['features']['name'][:])
            var = pd.DataFrame({'gene_name': h5fea}, index=h5fea)
        name = str(h5_obj["name"][0], encoding='utf8')

        adata = ad.AnnData(X, obs=obs, var=var)
        adata.obs_names = obs_names
        adata.obsm['spatial'] = obsm
        print("Loading Gene expression finished.")
        # create Images
        if loadImg and 'lowres' in h5img:
            adata.uns['spatial'] = {name: {
                'images': {'lowres': h5img['lowres'][:].T},
                'scalefactors': {'spot_diameter_fullres': sf['spot'][0],
                                 'tissue_hires_scalef': sf['hires'][0],
                                 'fiducial_diameter_fullres': sf['fiducial'][0],
                                 'tissue_lowres_scalef': sf['lowres'][0], },
                'metadata': {}
            }}
            if hires and 'hires' in h5img:
                adata.uns['spatial'][name]['images']['hires'] = h5img['hires'][:].T

            print("Loading spatial images finished.")
        # create svgs
        if 'is_HVG' in h5dat['features']:
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

        # whether you need to do subset
        if len(datas) > 1:
            subarcodes = bytes2str(h5_obj[h5data + '/barcodes'][:])
            adata = adata[subarcodes]

        # create deconv results
        if loadDeconv:
            h5dcv = h5_obj[h5data + '/deconv']
            if len(h5dcv.keys()) > 0:
                for dcv in h5dcv.keys():
                    # print(dcv)
                    shape = h5dcv[dcv]['shape']
                    weights = np.array(h5dcv[dcv]['weights']).reshape(shape[1], shape[0]).T
                    barcodes = bytes2str(h5dcv[dcv]['barcodes'][:])
                    cell_type = bytes2str(h5dcv[dcv]['cell_type'][:])
                    w = pd.DataFrame(weights, index=barcodes, columns=cell_type)
                    # print(w)
                    adata.obsm[dcv] = w
            print("Loading deconvolution data finished.")

        if 'abundance' in h5dat.keys():
            h5abd = h5dat['abundance']
            if len(h5abd.keys()) > 0:
                for abd in h5abd.keys():
                    # print(dcv)
                    shape = h5abd[abd]['shape']
                    weights = np.array(h5abd[abd]['weights']).reshape(shape[1], shape[0]).T
                    barcodes = bytes2str(h5abd[abd]['barcodes'][:])
                    cell_type = bytes2str(h5abd[abd]['cell_type'][:])
                    w = pd.DataFrame(weights, index=barcodes, columns=cell_type)
                    # print(w)
                    adata.obsm[abd] = w

    # making var and obs unique
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def Load_spt_sc_to_AnnData(sptFile, h5data='scRNA_seq'):
    h5_obj = h5.File(sptFile, 'r')
    # the scRNA-seq Data is saved in 'scRNA_seq', using Sparse Matrix
    scRNA_seq = h5_obj[h5data]
    data = scRNA_seq['data']
    indices = scRNA_seq['indices']
    indptr = scRNA_seq['indptr']
    shape = scRNA_seq['shape']
    idents = scRNA_seq['idents']

    # create X in AnnData
    X = csc_matrix((data, indices, indptr), shape=shape, dtype='int32').T
    # create obs in AnnData
    anno = bytes2str(idents['annotation'][:])
    bar = bytes2str(scRNA_seq['barcodes'][:])
    obs = pd.DataFrame({'annotation': anno},
                       index=bar)
    # create var in AnnData
    fea = bytes2str(scRNA_seq['features']['name'][:])
    var = pd.DataFrame(index=fea)
    # create AnnData
    adata = ad.AnnData(X, obs=obs, var=var)
    if 'HVGs' in scRNA_seq['features']:
        hvg = bytes2str(scRNA_seq['features']['HVGs'][:])
        adata.uns['HVGs'] = np.array(hvg, dtype='str')
    adata.var_names_make_unique()
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


# Save subset under existing h5datas
def Save_spt_from_Subset(sptFile:str, h5data:str, subset_name:str, refer_idx=None, force=False):
    # if it's required to save a reference to existing h5data, using ref_idx to select the subset.
    if refer_idx is not None:
        if type(refer_idx[0] == str):
            refer_idx = str2bytes(refer_idx)
        with h5.File(sptFile, 'a') as f:
            # create group
            if h5data + '/' + subset_name in f:
                if force:
                    for key in f[h5data + '/' + subset_name].keys():
                        f.__delitem__(h5data + '/' + subset_name + '/' + key)
                    f.__delitem__(h5data + '/' + subset_name)
                else:
                    print("Subset Reference {0} is existed.".format(h5data + '/' + subset_name))
                    return
            f.create_group(h5data + '/' + subset_name)
            f.create_dataset(h5data + '/' + subset_name + '/barcodes', data=refer_idx)
            f.create_dataset(h5data + '/' + subset_name + '/type', data=b"reference")
            f.create_group(h5data + '/' + subset_name + '/deconv')
            f.create_group(h5data + '/' + subset_name + '/idents')
        print("Subset Reference is saved in {0}".format(h5data + '/' + subset_name))
    return


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


def Save_tsv_from_spData(out_dir, spdata, filename='spt_data.tsv', transport=False, index_label='cell'):
    if not osp.exists(out_dir):
        os.mkdir(out_dir)
    mat = spdata.X.todense()
    spdata.var_names = spdata.var['gene_name']
    spdata.var_names_make_unique()
    cnt = pd.DataFrame(mat, index=spdata.obs_names,
                       columns=spdata.var_names)
    if transport:
        cnt = cnt.T
    opth_cnt = osp.join(out_dir, filename)
    cnt.to_csv(opth_cnt, sep='\t', header=True, index=True, index_label=index_label)


def Save_meta_from_spData(out_dir, spdata, sample_name, filename='spt_meta.tsv'):
    if not osp.exists(out_dir):
        os.mkdir(out_dir)
    meta = pd.DataFrame(spdata.obs[['array_row', 'array_col']])
    meta.columns = ['X', 'Y']
    if 'spatial' in spdata.uns:
        diameter = spdata.uns['spatial'][sample_name]['scalefactors']['spot_diameter_fullres']
    else:
        diameter = 100
    diameters = pd.DataFrame(np.repeat(diameter/2, len(meta)), columns=['Spot_radius'], index=meta.index)
    opth_meta = osp.join(out_dir, filename)
    meta = pd.concat([meta,diameters], axis=1)
    meta.to_csv(opth_meta, header=True, index=True, index_label="Barcodes")
    return


def Load_spt_to_stData(sptFile, count="matrix", hires=False, dtype=None):
    import stlearn as st
    adata = Load_spt_to_AnnData(sptFile, count, hires, dtype)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
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


def Save_spt_from_SEDR(sptFile, adata, h5data='matrix'):
    with h5.File(sptFile, 'a') as f:
        # save idents
        ident = list(np.array(adata.obs['clusters'], dtype='int32'))
        if h5data + '/idents/SEDR' not in f:
            f.create_dataset(h5data + '/idents/SEDR', data=ident, dtype='int32')
        else:
            del f[h5data + '/idents/SEDR']
            f.create_dataset(h5data + '/idents/SEDR', data=ident, dtype='int32')
    print("Clustering with `SEDR` finished, idents saved in /" + h5data + '/idents/SEDR')


def Save_spt_from_stlearn(sptFile, stdata, h5data='matrix'):
    with h5.File(sptFile, 'a') as f:
        # save idents
        ident = list(np.array(stdata.obs['louvain'], dtype='int32'))
        f.create_dataset(h5data + '/idents/stlearn', data=ident, dtype='int32')
    print("Clustering with `stlearn` finished, idents saved in /" + h5data + '/idents/stlearn')


def Save_spt_from_SpatialDE(sptFile, adata, h5data='matrix'):
    # FSV: Fraction of variance explained by spatial variation
    # P value: SpatialDE-logP Value
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
        name = list(str2bytes(adata.var.index))
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


def Save_spt_from_StereoScope(sptFile, Wfile, h5data='matrix', name='StereoScope'):
    W = pd.read_csv(Wfile, sep='\t', index_col=0).T
    shape = [W.shape[1], W.shape[0]]
    weights = np.array(W).flatten()
    cell_type = list(str2bytes(W.index))
    barcodes = list(str2bytes(W.columns))

    with h5.File(sptFile, 'a') as f:
        if h5data + '/deconv/' + name not in f:
            f.create_group(h5data + '/deconv/' + name)
        else:
            del f[h5data + '/deconv/' + name + '/weights']
            del f[h5data + '/deconv/' + name + '/shape']
            del f[h5data + '/deconv/' + name + '/cell_type']
            del f[h5data + '/deconv/' + name + '/barcodes']
            del f[h5data + '/deconv/' + name]
            f.create_group(h5data + '/deconv/' + name)
        f.create_dataset(h5data + '/deconv/' + name + '/weights', data=weights)
        f.create_dataset(h5data + '/deconv/' + name + '/cell_type', data=cell_type)
        f.create_dataset(h5data + '/deconv/' + name + '/barcodes', data=barcodes)
        f.create_dataset(h5data + '/deconv/' + name + '/shape', data=shape)
    print("Deconvolution with `StereoScope` finished, Weight saved in /" + h5data + '/deconv/' + name)


def Save_spt_from_Cell2Location(sptFile, adata, h5data='matrix'):
    W = adata.obsm['Cell2Location'].T
    shape = [W.shape[1], W.shape[0]]
    weights = np.array(W).flatten()
    cell_type = list(str2bytes(W.index))
    barcodes = list(str2bytes(W.columns))

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


def Save_spt_from_txt(sptFile, txt_file, name, transport=False):
    count = pd.read_csv(txt_file, sep='\t', index_col=0)
    barcodes = list(str2bytes(count.index))
    features = list(str2bytes(count.columns))
    if transport:
        count = csr_matrix(count).T
    else:
        count = csr_matrix(count)
    with h5.File(sptFile, 'a') as f:
        if name not in f:
            f.create_group(name)
        else:
            if name + '/idents/annotation' in f:
                del f[name + '/idents/annotation']
            if name + '/idents' in f:
                del f[name + '/idents']
            if name + '/features/name' in f:
                del f[name + '/features/name']
            if name + '/features' in f:
                del f[name + '/features']
            if name + '/barcodes' in f:
                del f[name + '/barcodes']
            if name + '/shape' in f:
                del f[name + '/shape']
            if name + '/indptr' in f:
                del f[name + '/indptr']
            if name + '/indices' in f:
                del f[name + '/indices']
            if name + '/data' in f:
                del f[name + '/data']
            if name in f:
                del f[name]
            f.create_group(name)
        f.create_dataset(name + '/data', data=count.data, dtype='float32')
        f.create_dataset(name + '/indices', data=count.indices, dtype='int')
        f.create_dataset(name + '/indptr', data=count.indptr, dtype='int')
        shape = count.shape
        f.create_dataset(name + '/shape', data=shape, dtype='int')
        f.create_dataset(name + '/barcodes', data=barcodes)
        f.create_group(name + '/features')
        f.create_dataset(name + '/features/name', data=features)
        f.create_group(name + '/idents')
        f.create_group(name + '/deconv')
    print("Save txt finished, matrix data saved in " + name)


def Save_spt_from_ref(sptFile, reference: pd.DataFrame, label: pd.DataFrame, name='sc-ref'):
    with h5.File(sptFile, 'a') as f:
        if name not in f:
            f.create_group(name)
        else:
            if name + '/idents/annotation' in f:
                del f[name + '/idents/annotation']
            if name + '/idents' in f:
                del f[name + '/idents']
            if name + '/features/name' in f:
                del f[name + '/features/name']
            if name + '/features' in f:
                del f[name + '/features']
            if name + '/barcodes' in f:
                del f[name + '/barcodes']
            if name + '/shape' in f:
                del f[name + '/shape']
            if name + '/indptr' in f:
                del f[name + '/indptr']
            if name + '/indices' in f:
                del f[name + '/indices']
            if name + '/data' in f:
                del f[name + '/data']
            if name in f:
                del f[name]
            f.create_group(name)
        ref = csr_matrix(reference).T
        f.create_dataset(name + '/data', data=ref.data, dtype='float32')
        f.create_dataset(name + '/indices', data=ref.indices, dtype='int')
        f.create_dataset(name + '/indptr', data=ref.indptr, dtype='int')
        shape = ref.shape
        f.create_dataset(name + '/shape', data=shape, dtype='int')
        f.create_dataset(name + '/barcodes', data=reference.index, dtype='int')
        f.create_group(name + '/features')
        refname = list(str2bytes(reference.columns))
        annotation = list(str2bytes(label['annotation']))
        f.create_dataset(name + '/features/name', data=refname)
        f.create_group(name + '/idents')
        f.create_dataset(name + '/idents/annotation', data=annotation)
    print("references finished, reference scRNA-seq data saved in " + name)


def Save_spt_from_corrcoef(sptFile, coef: pd.DataFrame, dcv_mtd, h5data='matrix'):
    weights = np.array(coef).flatten()
    cell_type = list(str2bytes(coef.index))
    shape = [coef.shape[1], coef.shape[0]]
    with h5.File(sptFile, 'a') as f:
        if h5data + "/coef" not in f:
            f.create_group(h5data + "/coef")
        if h5data + "/coef/" + dcv_mtd not in f:
            f.create_group(h5data + "/coef/" + dcv_mtd)
        else:
            del f[h5data + '/coef/' + dcv_mtd + '/weights']
            del f[h5data + '/coef/' + dcv_mtd + '/shape']
            del f[h5data + '/coef/' + dcv_mtd + '/cell_type']
            del f[h5data + '/coef/' + dcv_mtd]
            f.create_group(h5data + '/coef/' + dcv_mtd)
        f.create_dataset(h5data + '/coef/' + dcv_mtd + '/weights', data=weights)
        f.create_dataset(h5data + '/coef/' + dcv_mtd + '/cell_type', data=cell_type)
        f.create_dataset(h5data + '/coef/' + dcv_mtd + '/shape', data=shape)
    print("pearson's R with " + dcv_mtd + " finished, index saved in /" + h5data + '/coef/' + dcv_mtd)


def Save_spt_from_BestDict(sptFile, Best_dict: dict):
    col = ["cell-type", "F1-score", "domain", "dct", "dcv", "clu"]
    best_df = pd.DataFrame(columns=col)
    cell_types = Best_dict.keys()
    for ct in cell_types:
        Best_item = Best_dict[ct]
        mtds = Best_item[2].split('+')
        dct_mtd, dcv_mtd, clu_mtd = mtds[0], mtds[1], mtds[2]
        df_item = pd.DataFrame(data=[ct, Best_item[0], Best_item[1], dct_mtd, dcv_mtd, clu_mtd],
                               index=col, columns=[ct]).T
        best_df = pd.concat([best_df, df_item])
    with h5.File(sptFile, 'a') as f:
        if "/map_chains" not in f:
            f.create_group("/map_chains")
        else:
            for key in f['map_chains']:
                del f['map_chains/' + key]
        for cols in best_df.columns:
            dat = list(best_df[cols])
            if type(dat[0]) == str:
                f.create_dataset("map_chains/" + cols, data=str2bytes(dat))
            elif type(dat[0]) == tuple:
                f.create_dataset("map_chains/" + cols, data=tuple2bytes(dat))
            else:
                f.create_dataset("map_chains/" + cols, data=dat)
