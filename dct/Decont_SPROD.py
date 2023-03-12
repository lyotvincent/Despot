import os

import pandas as pd

from utils.io import Load_spt_to_AnnData, Save_tsv_from_spData, Save_meta_from_spData, sptInfo, Save_spt_from_txt


def Handle_Dependencies():
    print("Dependencies will be installed when Using SPROD for the first time.")
    # handle python dependencies

    # handle R dependencies
    pass


def sprod_pp(sptFile, tempdir='h5ads/temps'):
    sptinfo = sptInfo(sptFile)
    adata = Load_spt_to_AnnData(sptFile)
    print("Saving Count data to {0}.".format(tempdir))
    Save_tsv_from_spData(tempdir, adata, filename="Counts.txt", transport=False, index_label="Barcodes")
    print("Saving metadata to {0}.".format(tempdir))
    Save_meta_from_spData(tempdir, adata, sample_name=sptinfo.name, filename="Spot_metadata.csv")


def sprod_run(sprod_dir, out_dir, pythonPath):
    cmd = pythonPath + " SPROD/sprod.py " + sprod_dir + " " + out_dir
    os.system("{0} SPROD/setup.py install --user \n"
              "{1}".format(pythonPath, cmd))


def sprod_save(sptFile, out_dir):
    txt_file = os.path.join(out_dir, "sprod_Denoised_matrix.txt")
    Save_spt_from_txt(sptFile, txt_file, name="SPROD_mat", transport=True)
    os.remove(txt_file)