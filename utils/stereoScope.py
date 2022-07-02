from utils.common import *
from utils.io import *
from utils.scdg import Make_References, Save_tsv_from_ref

def Preprocess_hippo(pth, odir):
    ds = lp.connect(pth)

    celltype_1 = np.array([ '_'.join([str(x).replace(',','_'),str(y)]) for \
                           x,y in zip(ds.ca['Class'],ds.ca['Clusters']) ])

    celltype_2 = np.array([ '_'.join([str(x).replace(',','_'),str(y)]) for \
                           x,y in zip(ds.ca['Subclass'],ds.ca['Clusters']) ])


    new_ds_ca = dict(CellID = ds.ca['CellID'],
                     Clusters = ds.ca['Clusters'],
                     Subclass = ds.ca['Subclass'],
                     Celltype_1 = celltype_1,
                     Celltype_2 = celltype_2)

    new_ds_ra = dict(Gene = ds.ra['Gene'])
    new_basename = '_'.join(['mod',osp.basename(pth)])
    filename = os.path.join(odir,new_basename)
    lp.create(filename, ds[:,:], new_ds_ra, new_ds_ca)
    print(f"successfully created modifed loom-file >> {new_basename}")


def Preprocess_subsample(loom_pth, out_dir, label="Celltype_1",
                         lower_bound=25, upper_bound=250, add_filter=None):
    if not osp.exists(out_dir):
        os.mkdir(out_dir)
    tag = str(datetime.datetime.now())
    tag = re.sub('-| |:|\\.', '', tag)

    print(f"Unique Identifier for set >> {tag}")

    try:
        ds = lp.connect(loom_pth)
    except OSError:
        print(f'The file {loom_pth} either does not exist or',
              f" you might already have established a connection",
              f" to the file.")
        sys.exit(-1)

    if add_filter is not None:
        add_filter = ds.ca[add_filter[0]] == add_filter[1]
    else:
        add_filter = np.ones(ds.ca[label].shape[0]).astype(np.bool)

    uni_types = np.unique(ds.ca[label].flatten())

    use_cells = []
    for type in uni_types:
        zidx = np.where((ds.ca[label].flatten() == type) & add_filter)[0]
        nz = zidx.shape[0]

        if nz < lower_bound:
            print(f"{type} | was discarded due to insufficient number of cells")
            continue

        elif nz >= upper_bound:
            zidx = np.random.choice(zidx,
                                    size=upper_bound,
                                    replace=False,
                                    )

        fromz = zidx.shape[0]
        use_cells += zidx.tolist()

        print(f"{type} | Used {fromz} cells ")

    use_cells = np.array(use_cells)
    use_cells = np.sort(use_cells)

    new_cnt = ds[:, use_cells].T
    new_lbl = ds.ca[label].flatten()[use_cells]
    new_barcodes = ds.ca['CellID'][use_cells]
    _, idx = np.unique(new_barcodes, return_index=True)
    new_cnt = new_cnt[idx, :]
    new_lbl = new_lbl[idx]
    new_barcodes = new_barcodes[idx]
    genes = ds.ra['Gene']

    new_cnt = pd.DataFrame(new_cnt,
                           index=pd.Index(new_barcodes),
                           columns=pd.Index(genes),
                           )

    new_mta = pd.DataFrame(new_lbl,
                           index=pd.Index(new_barcodes),
                           columns=pd.Index(['bio_celltype']),
                           )

    name, counts = np.unique(new_lbl, return_counts=True)

    stats = pd.DataFrame(counts,
                         index=pd.Index(name),
                         columns=pd.Index(['members']),
                         )

    print(stats)
    print(f"Assembled Count Matrix with dimension >> nrow : {new_cnt.shape[0]} | ncol :  {new_cnt.shape[1]}")
    print(f"Assembled Meta File with >> nrow : {new_mta.shape[0]}")

    opth_cnt = osp.join(out_dir,
                        '.'.join([tag,
                                  'cnt_data.tsv',
                                  ]
                                 )
                        )

    opth_mta = osp.join(out_dir,
                        '.'.join([tag,
                                  'mta_data.tsv',
                                  ]
                                 )
                        )

    opth_sts = osp.join(out_dir,
                        '.'.join([tag,
                                  'stats.tsv',
                                  ]
                                 )
                        )

    new_cnt.to_csv(opth_cnt,
                   sep='\t',
                   header=True,
                   index=True,
                   index_label='cell',
                   )

    new_mta.to_csv(opth_mta,
                   sep='\t',
                   header=True,
                   index=True,
                   index_label='cell',
                   )

    stats.to_csv(opth_sts,
                 sep='\t',
                 header=True,
                 index=True,
                 index_label='cell',
                 )
    ds.close()


def myPreprocess(loom_pth, out_dir, standard_size=50, label='Celltype_1',add_filter=None):
    tag = str(datetime.datetime.now())
    tag = re.sub('-| |:|\\.', '', tag)

    print(f"Unique Identifier for set >> {tag}")
    try:
        ds = lp.connect(loom_pth)
    except OSError:
        print(f'The file {loom_pth} either does not exist or',
              f" you might already have established a connection",
              f" to the file.")
        sys.exit(-1)

    if add_filter is not None:
        add_filter = ds.ca[add_filter[0]] == add_filter[1]
    else:
        add_filter = np.ones(ds.ca[label].shape[0]).astype(np.bool)
    uni_types = np.unique(ds.ca[label].flatten())

    use_cells = []
    for type in uni_types:
        zidx = np.where((ds.ca[label].flatten() == type) & add_filter)[0]
        nz = zidx.shape[0]

        if nz < standard_size:
            print(f"{type} | was discarded due to insufficient number of cells")
            continue

        else:
            zidx = np.random.choice(zidx,
                                    size=standard_size,
                                    replace=False,
                                    )

        fromz = zidx.shape[0]
        use_cells += zidx.tolist()

        print(f"{type} | Used {fromz} cells ")
    use_cells = np.array(use_cells)
    use_cells = np.sort(use_cells)
    new_cnt = ds[:, use_cells].T
    new_lbl = ds.ca[label].flatten()[use_cells]
    new_barcodes = ds.ca['CellID'][use_cells]
    _, idx = np.unique(new_barcodes, return_index=True)
    new_cnt = new_cnt[idx, :]
    new_lbl = new_lbl[idx]
    new_barcodes = new_barcodes[idx]
    genes = ds.ra['Gene']

    new_cnt = pd.DataFrame(new_cnt,
                           index=pd.Index(new_barcodes),
                           columns=pd.Index(genes),
                           )

    new_mta = pd.DataFrame(new_lbl,
                           index=pd.Index(new_barcodes),
                           columns=pd.Index(['bio_celltype']),
                           )

    name, counts = np.unique(new_lbl, return_counts=True)

    stats = pd.DataFrame(counts,
                         index=pd.Index(name),
                         columns=pd.Index(['members']),
                         )

    print(stats)
    print(f"Assembled Count Matrix with dimension >> nrow : {new_cnt.shape[0]} | ncol :  {new_cnt.shape[1]}")
    print(f"Assembled Meta File with >> nrow : {new_mta.shape[0]}")
    opth_cnt = osp.join(out_dir,
                        '.'.join([tag,
                                  'cnt_data.tsv',
                                  ]
                                 )
                        )

    opth_mta = osp.join(out_dir,
                        '.'.join([tag,
                                  'mta_data.tsv',
                                  ]
                                 )
                        )

    opth_sts = osp.join(out_dir,
                        '.'.join([tag,
                                  'stats.tsv',
                                  ]
                                 )
                        )

    new_cnt.to_csv(opth_cnt,
                   sep='\t',
                   header=True,
                   index=True,
                   index_label='cell',
                   )

    new_mta.to_csv(opth_mta,
                   sep='\t',
                   header=True,
                   index=True,
                   index_label='cell',
                   )

    stats.to_csv(opth_sts,
                 sep='\t',
                 header=True,
                 index=True,
                 index_label='cell',
                 )
    ds.close()

# just use easy down sampling to get scRNA-seq references
def Easy_Sample(scdata, standard_size=25, add_filter=None):
    # High Variable Genes, saved in scdata.uns['HVGs']
    HVGs = scdata.uns['HVGs']

    # filter genes for other conditions
    if add_filter is None:
        add_filter = np.ones(len(scdata)).astype(np.bool)

    # all the cell types
    uni_types = np.unique(scdata.obs['annotation'])
    use_cells = []
    for type in uni_types:
        zidx = np.where((scdata.obs['annotation'] == type) & add_filter)[0]
        nz = zidx.shape[0]

        if nz < standard_size:
            print(f"{type} | was discarded due to insufficient number of cells")
            continue
        else:
            zidx = np.random.choice(zidx,size=standard_size,
                                    replace=False)
        fromz = zidx.shape[0]
        use_cells += zidx.tolist()

        print(f"{type} | Used {fromz} cells ")

    use_cells = np.array(use_cells)
    use_cells = np.sort(use_cells)
    scdata0 = scdata[use_cells, HVGs]
    return scdata0

def StereoScope_pp_EasySample(sptFile, out_dir='h5ads', h5data='matrix', standard_size=25 ,add_filter=None):
    # load spatial data
    spdata = Load_spt_to_AnnData(sptFile, h5data)
    # load scRNA-seq data
    scdata = Load_spt_sc_to_AnnData(sptFile)
    scdata0 = Easy_Sample(scdata, standard_size, add_filter)

    # save scRNA-seq data to tsv
    Save_tsv_from_scData(out_dir, scdata0)
    # save ST data to tsv
    Save_tsv_from_spData(out_dir, spdata)

def StereoScope_pp_VAE(sptFile, tempdir='h5ads', h5data='matrix', standard_size=25, add_filter=None, name="temp"):
    # load spatial data
    spdata = Load_spt_to_AnnData(sptFile, h5data)
    out_dir = tempdir + "/references/" + name
    Save_tsv_from_spData(out_dir, spdata)

    # load scRNA-seq data
    ref, lbl = Make_References(sptFile, ref_size=standard_size)
    Save_tsv_from_ref(ref, lbl, tempdir, name)


def StereoScope_cmd(sc_cnt, sc_labels, st_cnt,
                    sce=75000, ste=75000, out='h5ads/res', num=5000,
                    gpu=True, scb=100, stb=100):
    cmd = "stereoscope run --sc_cnt " + sc_cnt + \
              " --sc_labels " + sc_labels + \
              " --st_cnt " + st_cnt + \
              " -sce " + str(sce) + \
              " -ste " + str(ste) + \
              " -o " + out + \
              " -n " + str(num) + \
              " -scb " + str(scb) + \
              " -stb " + str(stb) + \
              " --gpu "
    return cmd

def StereoScope_run(stereo_dir, out_dir='h5ads/res',
                    scdata_epoch=5000, stdata_epoch=5000,
                    scdata_batch=100, stdata_batch=100,
                    num=3000, gpu=True):
    sc_cnt  = osp.join(stereo_dir, 'cnt_data.tsv')
    sc_labels = osp.join(stereo_dir, 'mta_data.tsv')
    st_cnt = osp.join(stereo_dir, 'spt_data.tsv')

    cmd = StereoScope_cmd(sc_cnt, sc_labels, st_cnt, out = out_dir,
                          sce = scdata_epoch, ste = stdata_epoch,
                          scb = scdata_batch, stb = stdata_batch,
                          num=num, gpu=gpu)
    os.system(cmd)
