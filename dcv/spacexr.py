import time

import pandas as pd
from utils.common import *
from utils.io import *
from utils.scdg import Make_References, Easy_Sample


def spacexr_pp_VAE(sptFile, standard_size=25, add_filter=None):
    # load scRNA-seq data
    start = time.time()
    ref, lbl = Make_References(sptFile, ref_size=standard_size)
    Save_spt_from_ref(sptFile, ref, lbl, name="sc-ref")
    end = time.time()
    print("pp_vae using: {0}s.".format(end-start))


def spacexr_pp_EasySample(sptFile, standard_size=25, add_filter=None):
    # load scRNA-seq data
    scdata = Load_spt_sc_to_AnnData(sptFile)
    scdata0 = Easy_Sample(scdata, standard_size, add_filter)
    ref = pd.DataFrame(scdata0.X.todense(), columns=scdata0.var_names)
    lbl = pd.DataFrame(scdata0.obs['annotation'], copy=True)
    lbl.index = ref.index
    Save_spt_from_ref(sptFile, ref, lbl, name="sc-ref-es")

import weka.classifiers