import numpy as np
from numpy import transpose
import xarray as xr

from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair, ALL_TIMES, destagger)
from scipy.stats import pearsonr, kendalltau
from xarray import corr
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy

wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")


# trad1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
# wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")


#####SELECT A 3D VARIABLE
VAR1 = getvar(wrf_file, "W", timeidx=ALL_TIMES)
K = getvar(wrf_file, "RAINC", timeidx=0)
VAR2 = getvar(wrf_file, "QCLOUD", timeidx=ALL_TIMES)


level = [5, 10, 15, 17, 18, 20, 22, 23, 24]
'''
LEVEL= [17]
for LEVEL in LEVEL:
    for i in range(0,299,1):
        for j in range(0,231,1):
            VAR1 = VAR1[:,LEVEL,i,j]*1000
            VAR2 = VAR2[:,LEVEL,i,j]
            #corr_pearson, _ = pearsonr(VAR1, VAR2)
            corr_kendall, _ = kendalltau(VAR1, VAR2)

'''


font = {'family': 'serif',
        'color': 'black',
        'weight': 'bold',
        'size': 12,
        }

for L in level:
    var1 = VAR1[:,L,:,:]
    var2 = VAR2[:,L,:,:]
    CORR = corr(var1, var2, dim="Time")
    lats, lons = latlon_coords(CORR)
    cart_proj = get_cartopy(K)
    ax = plt.axes(projection=cart_proj)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    # ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.add_feature(cartopy.feature.LAKES)
    # ax.add_feature(cartopy.feature.RIVERS)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='--')
    # ax.set_xlim(cartopy_xlim(RAINNC766_1))
    # ax.set_ylim(cartopy_ylim(RAINNC766_1))
    # ax.gridlines(color="black", linestyle="dotted")
    # v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),30) ##COLORBAR LIMITS
    v = np.linspace(-1, 1, 30)
    plt.contourf(lons, lats, CORR, v, cmap=plt.get_cmap("nipy_spectral"), antialiased=False,
                 transform=cartopy.crs.PlateCarree())
    plt.colorbar(ax=ax, shrink=.81)
    plt.title("PCC(QCLOUD-W) at level=" + str(L), weight='bold', fontdict=font)
    # plt.figure(figsize=(8,6))
    plt.savefig("PCC-QCLOUD-W AT LEVEL =" + str(L) + ".png", dpi=300)
    plt.show()


