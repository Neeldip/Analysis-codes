import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
from xarray import open_dataset
import pandas as pd

#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
ds = open_dataset(r"H:\WRF_Chem_Output\201\August\wrfout_d01_2018-08-10_00%3A00%3A00")
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\August\wrfout_d01_2018-08-10_00%3A00%3A00")

time = ALL_TIMES
#time = 60
#level = range(0,25,1)
HGT = ds.variables["HGT"][0,:,:]
C = getvar(wrf_file, "HGT", timeidx=1)

#ds = pd.DataFrame(CLD[24:768,156,136])
#ds.to_csv("pblh_ghy_hong.csv")
#print(ds)
#if time == ALL_TIMES:
#        CLD = CLD.mean("Time")
        #CLD1 = CLD1.mean("Time")
#print(CLD[156,136])
 #lats, lons = w.latlon_coords(CLD)

CLD_DIFF = (C)#- CLD1)*1000
lats,lons =w.latlon_coords(C)
#print(lats)
#print(lons)
#font = {'family': 'serif',
#     'color':  'black',
#     'weight': 'normal',
#     'size': 12,
#       }
cart_proj = get_cartopy(C)
ax = plt.axes(projection=cart_proj)
states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8)
#ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.COASTLINE)
#ax.add_feature(cartopy.feature.LAKES)
#ax.add_feature(cartopy.feature.RIVERS)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
#ax.set_xlim(cartopy_xlim(RAINNC766_1))
#ax.set_ylim(cartopy_ylim(RAINNC766_1))
#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),30) ##COLORBAR LIMITS
v = np.linspace(0,7000,10)
#plt.figure(figsize=(8,6))
plt.contourf(lons, lats, HGT,v, cmap=plt.get_cmap("nipy_spectral"),antialiased=False,
                 transform=cartopy.crs.PlateCarree(),extend="max")
cb=plt.colorbar(ax=ax, shrink=.93,ticks=v)
#cb.set_label("PBLH (m)",size=12)
#plt.title("(d) MYNN3-April",weight='normal',fontdict=font)
#plt.title("(b) ACM2")
#plt.savefig("PBLH_QNSE_April.png",dpi=300,bbox_inches='tight')
plt.show()