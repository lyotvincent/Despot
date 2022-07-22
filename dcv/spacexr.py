from utils.common import *
from utils.io import *
from utils.scdg import Make_References, Save_tsv_from_ref

def spacexr_pp_VAE(sptFile, tempdir='h5ads', h5data='matrix', standard_size=25, add_filter=None, name="temp"):
    # load spatial data
    spdata = Load_spt_to_AnnData(sptFile, h5data)
    out_dir = tempdir + "/references/" + name
    Save_tsv_from_spData(out_dir, spdata)

    # load scRNA-seq data
    ref, lbl = Make_References(sptFile, ref_size=standard_size)
    Save_tsv_from_ref(ref, lbl, tempdir, name)