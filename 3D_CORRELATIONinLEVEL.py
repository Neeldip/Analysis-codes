import numpy as np
from numpy import transpose
import xarray as xr
from scipy.stats import spearmanr
from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES,destagger)
from scipy.stats import pearsonr,kendalltau
from xarray import corr
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

#####SELECT A 3D VARIABLE
VAR1 = getvar(trad, "QVAPOR", timeidx=ALL_TIMES)
K = getvar(wrf_file, "RAINC", timeidx=0)
VAR2 = getvar(wrf_file, "CLDFRA", timeidx=ALL_TIMES)
VAR1=VAR1.astype("float64")


CORR= corr(VAR1,VAR2,dim="Time")
lats, lons = latlon_coords(CORR)
font = {'family': 'serif',
        'color': 'black',
        'weight': 'bold',
        'size': 12,
        }

for L in range(0,44,1):
    C=CORR[L,:,:]
    cart_proj = get_cartopy(K)
    ax = plt.axes(projection=cart_proj)
    states = NaturalEarthFeature(category="cultural", scale="50m",facecolor="none",name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    # ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.add_feature(cartopy.feature.LAKES)
    # ax.add_feature(cartopy.feature.RIVERS)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
    # ax.set_xlim(cartopy_xlim(RAINNC766_1))
    # ax.set_ylim(cartopy_ylim(RAINNC766_1))
    # ax.gridlines(color="black", linestyle="dotted")
    # v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),30) ##COLORBAR LIMITS
    v = np.linspace(-1.0, 1.0, 21)
    #plt.figure(figsize=(11, 7))
    plt.contourf(lons, lats, C, v, cmap=plt.get_cmap("jet"), transform=cartopy.crs.PlateCarree(),extend='both')
    plt.colorbar(ax=ax, shrink=.81)
    #plt.title("PCC(QCLOUD-W) at level=" + str(L), weight='bold', fontdict=font)

    plt.savefig(r"F:\WRF-CHEM ANALYSIS\WRF level analysis\CLDFRA-TEMP CORR\4xBC\CLDFRA-TEMP AT LEVEL =" + str(L) + ".png", dpi=600,bbox_inches='tight')
    plt.show()

VAR1=None
VAR2=None