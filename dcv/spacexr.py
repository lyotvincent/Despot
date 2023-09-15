import time

import pandas as pd
from utils.common import *
from utils.io import *
from utils.scdg import Make_References, Easy_Sample


def spacexr_pp_VAE(smdFile, standard_size=10, add_filter=None):
    # load scRNA-seq data
    start = time.time()
    ref, lbl = Make_References(smdFile, ref_size=standard_size)
    Save_smd_from_ref(smdFile, ref, lbl, name="sc-ref")
    end = time.time()
    print("pp_vae using: {0}s.".format(end-start))


def spacexr_pp_EasySample(smdFile, standard_size=25, add_filter=None):
    # load scRNA-seq data
    scdata = Load_smd_sc_to_AnnData(smdFile)
    scdata.var_names_make_unique()
    scdata0 = Easy_Sample(scdata, standard_size, add_filter)
    ref = pd.DataFrame(scdata0.X.todense(), columns=scdata0.var_names)
    lbl = pd.DataFrame(scdata0.obs['annotation'], copy=True)
    lbl.index = bytes2str(ref.index)
    lbl['annotation'] = bytes2str(lbl['annotation'])
    Save_smd_from_ref(smdFile, ref, lbl, name="sc-ref-es")

import weka.classifiers