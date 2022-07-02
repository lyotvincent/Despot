import esda
import pandas as pd
import geopandas as gpd
from geopandas import GeoDataFrame
import libpysal as lps
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from utils.io import *

adata = Load_spt_to_AnnData("h5ads/FFPE_Mouse_Kidney.h5spt")
geometry = [Point(xy) for xy in zip(adata.obs['array_row'], adata.obs['array_col'])]
ad_gdf = GeoDataFrame(adata.obsm['spacexr'], geometry=geometry)
wq = lps.weights.Queen.from_dataframe(ad_gdf)
wq.transform = 'r'
MI = list(map(lambda x: esda.moran.Moran(ad_gdf.iloc[:, x], wq), list(range(14))))
I = list(map(lambda x: MI[x].I, list(range(14))))
print(I)

bcell = ad_gdf.iloc[:, 0]
bcell_lag = lps.weights.lag_spatial(wq, bcell)
b, a = np.polyfit(x=np.array(bcell), y=np.array(bcell_lag), deg=1)
plt.plot(bcell, bcell_lag, '.', color='firebrick')

 # dashed vert at mean of the price
plt.vlines(bcell.mean(), bcell_lag.min(), bcell_lag.max(), linestyle='--')
 # dashed horizontal at mean of lagged price
plt.hlines(bcell_lag.mean(), bcell.min(), bcell.max(), linestyle='--')

# red line of best fit using global I as slope
plt.plot(bcell, a + b*bcell, 'r')
plt.title('Moran Scatterplot')
plt.ylabel('Spatial Lag of Price')
plt.xlabel('Price')
plt.show()
# y = ad_gdf.iloc[:, 0]
# ylag = lps.weights.lag_spatial(wq, y)
# mi = esda.moran.Moran(y, wq)
li = esda.moran.Moran_Local(bcell, wq)
sig = 1 * (li.p_sim < 0.05)
hotspot = 1 * (sig * li.q==1)
coldspot = 3 * (sig * li.q==3)
doughnut = 2 * (sig * li.q==2)
diamond = 4 * (sig * li.q==4)
spots = hotspot + coldspot + doughnut + diamond
spot_labels = [ '0 ns', '1 hot spot', '2 doughnut', '3 cold spot', '4 diamond']
p_moransI = [spot_labels[i] for i in spots]
from matplotlib import colors
hmap = colors.ListedColormap([ 'lightgrey', 'red', 'lightblue', 'blue', 'pink'])
f, ax = plt.subplots(1, figsize=(9, 9))
ad_gdf.assign(cl=p_moransI).plot(column='cl', categorical=True,
        k=2, cmap=hmap, linewidth=0.1, ax=ax,
        edgecolor='white', legend=True)
ax.set_axis_off()
plt.show()


def p_moransI_purity(adata, dcv_mtd: str, clu: str, comp: str):
    geometry = [Point(xy) for xy in zip(adata.obs['array_row'], adata.obs['array_col'])]
    # first calculate global moransI
    ad_gdf = GeoDataFrame(adata.obsm[dcv_mtd], geometry=geometry)
    wq = lps.weights.Queen.from_dataframe(ad_gdf)
    wq.transform = 'r'
    y = ad_gdf.loc[comp]
    MI = esda.moran.Moran(y, wq)
    print("global Moran's I: {0}".format(MI.I))
    y_lag = lps.weights.lag_spatial(wq, y)
    # then calculate local moransI
    li = esda.moran.Moran_Local(bcell, wq)
    # using lisa cluster to find hotpoint, coldpoint, doughnut and diamond
    sig = 1 * (li.p_sim < 0.05)
    hotspot = 1 * (sig * li.q == 1)
    coldspot = 3 * (sig * li.q == 3)
    doughnut = 2 * (sig * li.q == 2)
    diamond = 4 * (sig * li.q == 4)
    spots = hotspot + coldspot + doughnut + diamond
    spot_labels = ['0 ns', '1 hot spot', '2 doughnut', '3 cold spot', '4 diamond']
    p_moransI = [spot_labels[i] for i in spots]
    # select the point in lisa cluster
    adata.obs['p_moransI_'+comp] = p_moransI
    adata0 = adata[adata.obs['p_moransI'+comp] > 0]
    # calculate the purity with the segmentation


    # scale the local moransI purity
    pass

