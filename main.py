from utils.common import *
from utils.io import *
from utils.preprocess import *
from utils.geo import *


def Spt_init(sptFile: str, force: bool = False):
    # whether the file exists
    if os.path.isfile(sptFile) and force is False:
        with h5.File(sptFile, 'r') as f:
            if 'matrix' in f:
                print("SptFile:" + sptFile + " has initialized, skip it.")
                return
    print("Initializing sptFile...")
    cmd = "Rscript utils/Init.R"
    os.system(cmd)
    print("done.")


# we provide 2 methods for decontamination
def Pip_decont(sptFile, method='none', force=False):
    # checking whether need to do decontamination
    method_list = ['none', 'SpotClean', 'SPCS']
    with h5.File(sptFile, 'r') as f:
        if method not in method_list:
            print("Decontamination method:{0} has not been supported yet, default `none`.".format(method))
            return
        if (method + "_mat" in f) and (force == False):
            print(method + " has done, skip it.")
            return

    # do decontamination
    if method == 'SpotClean':
        print("Decontamination method: SpotClean.")
        os.system("Rscript dct/Decont_SpotClean.R")
    elif method == 'SPCS':
        print("Decontamination method: SPCS.")
        os.system("Rscript dct/Decont_SPCS.R")
    elif method == 'none':
        print("Decontamination method: none.")
        return


# we provide 7 inline methods for clustering and supports Squidpy and Giotto
def Pip_cluster(sptFile, h5data, method='leiden', force=False, tif=None):
    method_list = ['leiden', 'SpaGCN', 'stlearn', 'Seurat', 'BayesSpace', 'Giotto', 'Squidpy']
    with h5.File(sptFile, 'r') as f:
        if method not in method_list:
            print("Clustering method:{0} has not been supported yet, default `leiden`.".format(method))
            method = 'leiden'
        if (h5data + '/idents/' + method in f) and (force == False):
            print("Clustering with " + method + " has done, skip it.")
            return
    if method == 'leiden':
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        adata = Load_spt_to_AnnData(sptFile, h5data)
        adata = Remove_mito_genes(adata)
        adata = Filter_genes(adata)
        adata = Spatial_Cluster_Analysis(adata)
        Save_spt_from_leiden(sptFile, adata, h5data)
        # adata = Differential_Expression_Analysis(adata)
    elif method == 'SpaGCN':
        from clu.Cluster_SpaGCN import Spatial_Cluster_Analysis_SpaGCN
        adata = Load_spt_to_AnnData(sptFile, h5data)
        adata = Spatial_Cluster_Analysis_SpaGCN(adata)
        Save_spt_from_SpaGCN_clu(sptFile, adata, h5data)
        # adata = Differential_Expression_Analysis(adata)
    elif method == 'BayesSpace':
        # Using R scripts
        os.system("Rscript clu/Cluster_BayesSpace.R")
    elif method == "Seurat":
        # Using R scripts
        os.system("Rscript clu/Cluster_Seurat.R")
    elif method == "stlearn":
        from clu.Cluster_stlearn import Spatial_Cluster_Analysis_stSME
        stdata = Load_spt_to_stData(sptFile, hires=hires, dataPath=dataPath)
        stdata = Preprocess_stLearn(stdata)
        stdata = Spatial_Cluster_Analysis_stSME(stdata)
        Save_spt_from_stlearn(sptFile, stdata, h5data)
    elif method == "Giotto":
        os.system("Rscript clu/Cluster_Giotto.R")
    elif method == "Squidpy":
        from clu.Cluster_leiden import Spatial_Cluster_Analysis
        from clu.Cluster_Squidpy import Spatial_Cluster_Analysis_squidpy
        adata = Load_spt_to_AnnData(sptFile, h5data)
        if tif is None:
            print("Squidpy needs input image for cluster.Default leiden.")
            adata = Remove_mito_genes(adata)
            adata = Filter_genes(adata)
            adata = Spatial_Cluster_Analysis(adata)
            Save_spt_from_leiden(sptFile, adata, h5data)
        else:
            adata = Spatial_Cluster_Analysis_squidpy(adata, tif)
            Save_spt_from_Squidpy_clu(sptFile, adata, h5data)


def Pip_estimate(sptFile, h5data='matrix', method='SpatialDE', force=False):
    method_list = ['SpatialDE', 'SpaGCN', 'SPARK', 'Seurat', 'BayesSpace', 'Giotto', 'Squidpy']
    with h5.File(sptFile, 'r') as f:
        if method not in method_list:
            print("Clustering method:{0} has not been supported yet, default `leiden`.".format(method))
            method = 'SpatialDE'
        if (h5data + '/features/is_HVG/' + method in f) and (force == False):
            print("Estimating with " + method + " has done, skip it.")
            return
    if method == 'SpatialDE':
        from est.Estim_SpatialDE import HDG_Analysis_by_SpatialDE
        adata = Load_spt_to_AnnData(sptFile, h5data)
        adata.uns['HDGs'] = HDG_Analysis_by_SpatialDE(adata)
        Save_spt_from_SpatialDE(sptFile, adata, h5data)
    elif method == 'SpaGCN':
        from est.Estim_SpaGCN import HDG_Analysis_by_SpaGCN
        adata = Load_spt_to_AnnData(sptFile, h5data)
        adata = HDG_Analysis_by_SpaGCN(adata)
        Save_spt_from_SpaGCN_est(sptFile, adata, h5data)
    elif method == 'SPARK':
        # Using R scripts
        os.system("Rscript est/Estim_SPARK.R")
    elif method == 'trendsceek':
        # Using R scripts
        os.system("Rscript est/Estim_trendsceek.R")
    elif method == 'Giotto':
        # Using R scripts
        os.system("Rscript est/Estim_Giotto.R")


def Pip_deconv(sptFile, h5data='matrix', method="stereoScope", name='temp', force=False):
    method_list = ['Cell2Location', 'SPOTlight', 'SPOTlight_vae', 'SPOTlight_es','spacexr', 'spacexr_es',
                   'stereoScope', 'Seurat', 'Seurat_es', 'Giotto', 'stereoScope_es', 'stereoScope_na']
    with h5.File(sptFile, 'r') as f:
        if method not in method_list:
            print("Deconvolution method:{0} has not been supported yet, default `stereoScope`.".format(method))
            method = 'stereoScope'
        if (h5data + '/deconv/' + method in f) and (force == False):
            print("Deconvolution with " + method + " has done, skip it.")
            return
    if method == 'spacexr':
        do_vae = True
        with h5.File(sptFile, 'r') as f:
            if ('sc-ref' in f) and (force == False):
                do_vae = False
        if do_vae:
            from dcv.spacexr import spacexr_pp_VAE
            spacexr_pp_VAE(sptFile)
        os.system("Rscript dcv/Deconv_spacexr.R")
    elif method == 'spacexr_es':
        from dcv.spacexr import spacexr_pp_EasySample
        spacexr_pp_EasySample(sptFile)
        os.system("Rscript dcv/Deconv_spacexr_es.R")
    elif method == 'SPOTlight':
        # Using R scripts
        os.system("Rscript dcv/Deconv_SPOTlight.R")
    elif method == 'SPOTlight_es':
        do_es = True
        with h5.File(sptFile, 'r') as f:
            if ('sc-ref-es' in f) and (force == False):
                do_es = False
        if do_es:
            from dcv.spacexr import spacexr_pp_EasySample
            spacexr_pp_EasySample(sptFile)
        # Using EasyEnsample for SPOTlight reference
        os.system("Rscript dcv/Deconv_SPOTlight_es.R")
    elif method == 'SPOTlight_vae':
        do_vae = True
        with h5.File(sptFile, 'r') as f:
            if ('sc-ref' in f) and (force == False):
                do_vae = False
        if do_vae:
            from dcv.spacexr import spacexr_pp_VAE
            spacexr_pp_VAE(sptFile)
        # Using VAE for SPOTLight reference
        os.system("Rscript dcv/Deconv_SPOTlight_vae.R")
    elif method == 'stereoScope':
        from dcv.stereoScope import StereoScope_pp_VAE, StereoScope_run, StereoScope_pp_EasySample
        StereoScope_pp_VAE(sptFile, tempdir="h5ads", h5data=h5data, name=name)
        # StereoScope_pp_EasySample(sptFile, h5data=h5data)
        StereoScope_run(stereo_dir="h5ads/references/" + name, out_dir="h5ads/" + name + "_res")
        Wfile = os.listdir("h5ads/" + name + "_res/spt_data")[0]
        Wfile = "h5ads/" + name + "_res/spt_data/" + Wfile
        Save_spt_from_StereoScope(sptFile, Wfile, h5data, name='StereoScope')
        os.remove(Wfile)
    elif method == 'stereoScope_es':
        # stereoScope_na represents using easy sample in single-cell data
        from dcv.stereoScope import StereoScope_pp_VAE, StereoScope_run, StereoScope_pp_EasySample
        StereoScope_pp_EasySample(sptFile, h5data=h5data)
        StereoScope_run(stereo_dir="h5ads/references/" + name, out_dir="h5ads/" + name + "_res")
        Wfile = os.listdir("h5ads/" + name + "_res/spt_data")[0]
        Wfile = "h5ads/" + name + "_res/spt_data/" + Wfile
        Save_spt_from_StereoScope(sptFile, Wfile, h5data, name='StereoScope_es')
    elif method == 'stereoScope_na':
        # stereoScope_na represents no preprocessing in single-cell data
        from dcv.stereoScope import StereoScope_pp_na, StereoScope_run
        StereoScope_pp_na(sptFile, h5data=h5data)
        StereoScope_run(stereo_dir="h5ads/references/" + name, out_dir="h5ads/" + name + "_res")
        Wfile = os.listdir("h5ads/" + name + "_res/spt_data")[0]
        Wfile = "h5ads/" + name + "_res/spt_data/" + Wfile
        Save_spt_from_StereoScope(sptFile, Wfile, h5data, name='StereoScope_na')
    elif method == 'Cell2Location':
        from dcv.cell2location import Cell2Location_run
        adata = Cell2Location_run(sptFile)
        Save_spt_from_Cell2Location(sptFile, adata, h5data)
    elif method == 'Seurat':
        os.system("Rscript dcv/Deconv_Seurat.R")
    elif method == 'Giotto':
        os.system("Rscript dcv/Deconv_Giotto.R")


def Pip_cluVisual(sptFile, h5data, imgPath):
    adata = Load_spt_to_AnnData(sptFile, h5data)
    figa, axsa = plt.subplots(2, 3, figsize=(9, 6), constrained_layout=True)
    sc.pl.spatial(adata, img_key="lowres", color="Seurat", size=1.5, ax=axsa[0, 0], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="Giotto", size=1.5, ax=axsa[0, 1], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="ground_truth", size=1.5, ax=axsa[0, 2], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="BayesSpace", size=1.5, ax=axsa[1, 0], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="Squidpy", size=1.5, ax=axsa[1, 1], show=False)
    sc.pl.spatial(adata, img_key="lowres", color="stlearn", size=1.5, ax=axsa[1, 2], show=False)
    figa.savefig(imgPath)


def Pip_benchmark(sptFile, h5data):
    adata = Load_spt_to_AnnData(sptFile, h5data)
    clu = pd.DataFrame(adata.obs[['Seurat', 'SpaGCN', 'stlearn', 'BayesSpace', 'leiden', 'Giotto']], dtype='category')
    from sklearn import metrics
    ground_truth = adata.obs['ground_truth'].astype('category')
    benchmarks = pd.DataFrame(columns=clu.columns)
    NMI = [metrics.normalized_mutual_info_score(clu[method], ground_truth) for method in clu.columns]
    ARI = [metrics.adjusted_rand_score(clu[method], ground_truth) for method in clu.columns]
    benchmarks.loc['NMI'] = NMI
    benchmarks.loc['ARI'] = ARI
    benchmarks.to_csv("h5ads/151673_" + h5data + ".csv", index=True, header=True)


cfg = Load_json_configs("params.json")
sptFile = cfg['sptFile']

# decode decontamination method
decont_method = cfg['Decontamination']

# decode clustering method
clust_method = cfg['Clustering']

# decode estimation method
estim_method = cfg['Estimation']

# decode deconvolution method
deconv_method = cfg['Deconvolution']

# decode dataPath
dataPath = cfg['dataPath']

# need hires?
hires = cfg['load_hires']


def spTRS_Decont(sptFile, cfg, force=False):
    methods = cfg['Decontamination']
    for method in methods:
        Pip_decont(sptFile, method=method, force=force)


def spTRS_Cluster(sptFile, cfg, force=False):
    h5datas = []
    # whether need Decontamination
    Decont = cfg['Decontamination']
    for dec in Decont:
        if dec == "SpotClean":
            h5data = 'SpotClean_mat'
            h5datas.append(h5data)
        elif dec == "SPCS":
            h5data = "SPCS_mat"
            h5datas.append(h5data)
        elif dec == "none":
            h5data = "matrix"
            h5datas.append(h5data)
        else:
            print("no matched Decontamination method, default `none`")

    h5datas.append("matrix")

    # decode dataPath
    dataPath = cfg['dataPath']
    methods = cfg['Clustering']
    tif = dataPath + '/' + cfg['imgPath'] + '/' + cfg['fullres_path']
    for h5data in h5datas:
        for method in methods:
            Pip_cluster(sptFile, h5data, method=method, tif=tif, force=force)


def spTRS_Estimate(sptFile, cfg, force=False):
    h5datas = []
    # whether need Decontamination
    Decont = cfg['Decontamination']
    for dec in Decont:
        if dec == "SpotClean":
            h5data = 'SpotClean_mat'
            h5datas.append(h5data)
        elif dec == "SPCS":
            h5data = "SPCS_mat"
            h5datas.append(h5data)
        elif dec == "none":
            h5data = "matrix"
            h5datas.append(h5data)
        else:
            print("no matched Decontamination method, default `none`")

    h5datas.append("matrix")

    methods = cfg['Estimation']
    for h5data in h5datas:
        for method in methods:
            Pip_estimate(sptFile, h5data, method=method, force=force)


def spTRS_Deconv(sptFile, cfg, name='temp', force=False):
    do_scRNA_seq = False
    with h5.File(sptFile, 'r') as f:
        if 'scRNA_seq' not in f:
            print("scRNA_seq data not in sptFile, Generating it...")
            do_scRNA_seq = True
        else:
            print("Found scRNA_seq data in sptFile.")
    if do_scRNA_seq:
        if cfg['TabulaType'] == '':
            print("Using scRNA-seq provided by users...")
            os.system("Rscript sc/Generate_scRNA-seq.R")
        else:
            print("Using Tabula Data...")
            os.system("Rscript sc/Generate_TabulaData.R")
    h5datas = []
    # whether need Decontamination
    Decont = cfg['Decontamination']
    for dec in Decont:
        if dec == "SpotClean":
            h5data = 'SpotClean_mat'
            h5datas.append(h5data)
        elif dec == "SPCS":
            h5data = "SPCS_mat"
            h5datas.append(h5data)
        elif dec == "none":
            h5data = "matrix"
            h5datas.append(h5data)
        else:
            print("no matched Decontamination method, default `none`")

    h5datas.append("matrix")

    methods = cfg['Deconvolution']
    for h5data in h5datas:
        for method in methods:
            Pip_deconv(sptFile, h5data, method=method, name=name, force=force)


def spTRS_Benchmark(sptFile, cfg, mode: str = "cluster", force: bool = False):
    mode_list = ['cluster', 'estimate', 'deconvolution']
    h5datas = []
    # whether need Decontamination
    Decont = cfg['Decontamination']
    for dec in Decont:
        if dec == "SpotClean":
            h5data = 'SpotClean_mat'
            h5datas.append(h5data)
        elif dec == "SPCS":
            h5data = "SPCS_mat"
            h5datas.append(h5data)
        elif dec == "none":
            h5data = "matrix"
            h5datas.append(h5data)
        else:
            print("no matched Decontamination method, default `none`")

    if len(h5datas) == 0:
        h5datas.append("matrix")

    bm_done = True
    with h5.File(sptFile, 'r') as f:
        if mode not in mode_list:
            print("Benchmark mode:{0} has not been supported yet".format(mode))
            return
        for h5data in h5datas:
            if (h5data + '/benchmark/' + mode not in f) or (force is True):
                bm_done = False
                break
    if bm_done is True:
        print("Benchmark with " + mode + " has done, skip it.")
        return
    if mode == "cluster":
        print("Computing all cluster benchmarks in sptFile.")
        os.system("Rscript clu/Benchmark_cluster.R")
    elif mode == "estimate":
        print("Computing all estimate benchmarks in sptFile.")
        os.system("Rscript est/Benchmark_estimate.R")


name = cfg['name']

Spt_init(sptFile=sptFile, force=True)
spTRS_Decont(sptFile, cfg)
spTRS_Cluster(sptFile, cfg)
spTRS_Estimate(sptFile, cfg, force=True)
spTRS_Deconv(sptFile, cfg)
# spTRS_Benchmark(sptFile, cfg, mode="cluster", force=False)
# Pip_cluVisual(sptFile, h5data, imgPath="151673_none_cluster.pdf")
# Cluster_Visualization(adata, type="spatial")
# sc.pl.spatial(adata, img_key="lowres", color="ground_truth", size=1.5)
