import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
def colorbar(folder, figsize=(2.2, 0.78)):
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    norm = colors.Normalize(vmin=0, vmax=1)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap="cividis"),
                 orientation="horizontal",
                 cax=ax,
                 drawedges=False,
                 ticks=[0, 1],

                 ticklocation='bottom')
    cb.set_label('Pearson\'s correlation', loc='center')
    cb.set_ticklabels(['<0', '1'], rotation=0)
    fig.savefig(folder + "/corr_cb1.svg")

def colorbar_v(folder, figsize=(0.78, 2.2), label='cell-type proportion', side='right', cmap='cividis', ticklabels=None):
    if ticklabels is None:
        ticklabels = ['0', '1']
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    norm = colors.Normalize(vmin=0, vmax=1)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                 orientation="vertical",
                 cax=ax,
                 drawedges=False,
                 ticks=[0, 1],
                 ticklocation=side)
    cb.set_label(label, loc='center', fontsize=16)
    cb.set_ticklabels(ticklabels, rotation=90)
    fig.savefig(folder + "/{0}_v{1}.svg".format(label, side))


def colorbar_h(folder, figsize=(2.2, 0.78), label='gene expression', side='bottom', cmap='viridis', ticklabels=None):
    if ticklabels is None:
        ticklabels = ['0', '1']
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["font.size"] = 14
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    norm = colors.Normalize(vmin=0, vmax=1)
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                 orientation="horizontal",
                 cax=ax,
                 drawedges=False,
                 ticks=[0, 1],
                 ticklocation=side)
    cb.set_label(label, loc='center', fontsize=16)
    cb.set_ticklabels(ticklabels, rotation=0)
    fig.savefig(folder + "/{0}_h{1}_{2}.svg".format(label, side, cmap))
