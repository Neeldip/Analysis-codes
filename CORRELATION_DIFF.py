import numpy as np
from numpy import transpose
import xarray as xr

from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair, ALL_TIMES, destagger)
from scipy.stats import pearsonr, kendalltau
from xarray import corr,DataArray
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy

wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset("G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
TRAD = Dataset(r"G:\WRF_Chem_Output\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
TRAD1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

# trad1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00").

K = getvar(wrf_file, "RAINC", timeidx=0)

#####SELECT A 3D VARIABLE
###VARIABLE 1##############
#V11 = getvar(TRAD, "TEMPERATURE", timeidx=ALL_TIMES)
V11 = getvar(wrf_file, "QCLOUD", timeidx=ALL_TIMES)
V12 = getvar(wrf_file, "CLDFRA", timeidx=ALL_TIMES)

###VARIABLE 2##############
#V21 = getvar(TRAD1, "TEMPERATURE", timeidx=ALL_TIMES)
V21 = getvar(wrf_file1, "QCLOUD", timeidx=ALL_TIMES)
V22 = getvar(wrf_file1, "CLDFRA", timeidx=ALL_TIMES)

level = [5, 10,11,12,13,14, 15, 17, 18, 20, 22, 23, 24,25]
#level = [5]
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
lats, lons = latlon_coords(K)

for L in level:
    #####SELECT A 3D VARIABLE
    ###VARIABLE 1##############
    #VAR1_1 = wrf_file.variables['T'][:, L, :, :]
    #print(VAR1_1)
    #VAR1_2 = wrf_file.variables['P'][:, L, :, :]
    #print(VAR1_2)
    ###VARIABLE 2##############
    #VAR2_1 = wrf_file.variables['T'][:, L, :, :]
    #VAR2_2 = wrf_file.variables['P'][:, L, :, :]
    #print(VAR2_1)
    #print(VAR2_2)
    VAR1_1 = V11[:, L, :, :]
    VAR1_2 = V12[:, L, :, :]
    VAR2_1 = V21[:, L, :, :]
    VAR2_2 = V22[:, L, :, :]
    DIFF_1 = VAR1_1 - VAR2_1
    #print(DIFF_1)
    DIFF_2 = (VAR1_2 - VAR2_2)*100
    #print(DIFF_2)
    var1 = DataArray(DIFF_1)
    #print(var1)
    var2 = DataArray(DIFF_2)
    #print(var2)
    CORR = corr(var1, var2, dim="Time")
    #CORR = corr(var1, var2, dim="dim_0")
    #print(CORR)

    cart_proj = get_cartopy(K)
    ax = plt.axes(projection=cart_proj)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='--')
    v = np.linspace(-1, 1, 30)
    plt.contourf(lons, lats, CORR, v, cmap=plt.get_cmap("nipy_spectral"), antialiased=False,
                 transform=cartopy.crs.PlateCarree())
    plt.colorbar(ax=ax, shrink=.81)
    plt.title("PCC(QCLOUD-CLDFRA) at level=" + str(L), weight='bold', fontdict=font)
    # plt.figure(figsize=(8,6))
    plt.savefig("PCC-QCLOUD-CLDFRA AT LEVEL =" + str(L) + ".png", dpi=300)
    plt.show()
    


