import os

import pandas as pd

from utils.io import Load_smd_to_AnnData, Save_tsv_from_spData, Save_meta_from_spData, smdInfo, Save_smd_from_txt


def sprod_install():
    print("Dependencies will be installed when Using SPROD for the first time.")
    # handle python dependencies
    os.system("git clone https://github.com/yunguan-wang/SPROD \n"
              "pip install SPROD \n"
              "python SPROD/test_examples.py")


def sprod_pp(smdFile, tempdir='temps'):
    smdinfo = smdInfo(smdFile)
    adata = Load_smd_to_AnnData(smdFile)
    print("Saving Count data to {0}.".format(tempdir))
    Save_tsv_from_spData(tempdir, adata, filename="Counts.txt", transport=False, index_label="Barcodes")
    print("Saving metadata to {0}.".format(tempdir))
    Save_meta_from_spData(tempdir, adata, sample_name=smdinfo.name, filename="Spot_metadata.csv")


def sprod_run(sprod_dir, out_dir, pythonPath):
    cmd = pythonPath + " SPROD/sprod.py " + sprod_dir + " " + out_dir
    os.system("{0} SPROD/setup.py install --user \n"
              "{1}".format(pythonPath, cmd))


def sprod_save(smdFile, out_dir):
    txt_file = os.path.join(out_dir, "sprod_Denoised_matrix.txt")
    Save_smd_from_txt(smdFile, txt_file, name="SPROD_mat", transport=True)
    os.remove(txt_file)
