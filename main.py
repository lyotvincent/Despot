import json
import shutil
from utils.io import *
from utils.api import Create_h5datas, Despot_Decont, Despot_Cluster, Despot_Deconv
from utils.geo import Despot_Find_bestGroup, Despot_group_correlation,\
    Show_self_correlation, Show_bash_best_group, Gen_venn
from utils.check import *
print("Importing requirements...")

Check_Environments()
Check_Requirements({"anndata", "h5py","matplotlib", "numpy", "pandas", "scanpy", "scipy", "torch"})

def SMD_init(smdFile: str, force: bool = False):
    # whether the file exists
    if os.path.isfile(smdFile) and force is False:
        # check the corrections of sptFile
        if SMD_check_corrections(smdFile):
            return

    print("Initializing smdFile...")
    cmd = "Rscript utils/Init.R"
    os.system(cmd)
    print("Done.")


if __name__ == "__main__":
    cfg_list = os.listdir('configs')
    for cfg_name in cfg_list:
        if cfg_name == 'PDAC-A.json':
            cfg_path = 'configs/' + cfg_name
            shutil.copy(cfg_path, dst="params.json")
            cfg = Load_json_configs("params.json")
            smdFile = cfg['smdFile']
            name = cfg['name']
            platform = cfg['platform']

            # set python path
            pythonPath = sys.executable
            cfg["pythonPath"] = pythonPath

            # set R lib path
            cfg["R_library_Path"] = 'sptranr/R'

            # set working dictionary
            cfg["working_dir"] = os.path.dirname(os.path.abspath(__file__))

            # set venv
            cfg["venv"] = sys.prefix.split('/')[-1]
            cfg_json = json.dumps(cfg)

            # set running configs
            Save_json_configs(cfg, "params.json")
            Save_smd_from_configs(smdFile, items=cfg)
            # need hires?
            hires = cfg['load_hires']
            print("=========Despot Info============")
            print("smdFile:{0}".format(smdFile))
            print("dataset name:{0}".format(name))
            print("platform:{0}".format(platform))
            print("Using hires img: {0}".format(hires))
            print("=========Despot Start===========")
            SMD_init(smdFile=smdFile, force=True)
            Despot_Decont(smdFile, cfg)
            Despot_Cluster(smdFile, cfg)
            Despot_Deconv(smdFile, cfg)
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
