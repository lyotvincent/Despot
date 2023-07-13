import os

from utils.io import *

from utils.api import Create_h5datas, Despot_Deconv
from utils.geo import spTRS_Find_bestGroup, spTRS_group_correlation,\
    Show_self_correlation, Show_bash_best_group, Gen_venn


def Spt_init(sptFile: str, force: bool = False):
    # whether the file exists
    if os.path.isfile(sptFile) and force is False:
        # check the corrections of sptFile
        if Spt_check_corrections(sptFile):
            return

    print("Initializing sptFile...")
    cmd = "Rscript utils/Init.R"
    os.system(cmd)
    print("done.")


import shutil
cfg_list = os.listdir('configs')
for cfg_name in cfg_list:
    cfg_path = 'configs/' + cfg_name
    shutil.copy(cfg_path, dst="params.json")
    cfg = Load_json_configs("params.json")
    sptFile = cfg['sptFile']
    name = cfg['name']
    platform = cfg['platform']
    pythonPath = cfg["pythonPath"]
    Save_spt_from_configs(sptFile, items=cfg)
    # need hires?
    hires = cfg['load_hires']
    print("=========Despot Info============")
    print("sptFile:{0}".format(sptFile))
    print("dataset name:{0}".format(name))
    print("platform:{0}".format(platform))
    print("Using hires img: {0}".format(hires))
    print("=========Despot Start===========")
    # Spt_init(sptFile=sptFile)
    # Despot_Decont(sptFile, cfg)
    # Despot_Cluster(sptFile, cfg)
    Despot_Deconv(sptFile, cfg)
    print("=========Despot Finish==========")
# # Best_dict, Groups = Despot_Find_bestGroup(sptFile, greedy=3)
# folder = "FFPE_Mouse_K"
# comp, corr, p_val = Despot_group_correlation(sptFile, Best_dict, alpha=1)
# Show_self_correlation(corr, folder)
# Show_bash_best_group(sptFile, Best_dict=None, folder=folder)
# inter_genes = Gen_venn(sptFile, folder, Best_dict=None, show_genes=1, cell_filter=None)

# Despot_Benchmark(sptFile, cfg, mode="cluster", force=False)
# Pip_cluVisual(sptFile, h5data, imgPath="151673_none_cluster.pdf")
# Cluster_Visualization(adata, type="spatial")
# sc.pl.spatial(adata, img_key="lowres", color="ground_truth", size=1.5)
