import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
def colorbar(folder):
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(figsize=(0.8, 2), constrained_layout=True)
    norm = colors.Normalize(vmin=0, vmax=1)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap="cividis"),
                 orientation="vertical",
                 cax=ax,
                 drawedges=False,
                 ticks=[0.0, 1.0],
                 ticklocation='left')
    cb.set_label('cell-type proportion', loc='center')
    cb.set_ticklabels([0, 1], rotation=90)
    fig.savefig(folder + "/cell_type_cb2.svg")
