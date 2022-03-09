from netCDF4 import Dataset
from wrf import get_cartopy,getvar
import numpy as np
from xarray import open_mfdataset
from nco import *
import xarray as xr
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim
from pandas import DataFrame
import pandas as pd
#ds = open_mfdataset(r"I:\MERRA2_AOD\March\2018\5-10\*.nc",combine='by_coords',concat_dim='time')
#ds = open_mfdataset(r"I:\MERRA2_AOD\2001-2020\Monthwise\11.Nov\*.nc",combine='by_coords',concat_dim='time')
#ds = open_mfdataset(r"I:\MERRA2_AOD\2001-2020\Monthwise\3.Mar\*.nc",combine='by_coords',concat_dim='time')
ds = open_mfdataset(r"I:\MERRA2_AOD\2001-2020\Monthwise\12.Dec\*.nc",combine='by_coords',concat_dim='time')
#ds1 =ds.mean(dim="time")
lats = ds.variables["lat"][:]
lons = ds.variables["lon"][:]
list1=[23.883,26.1,27.483,24.83,26.670,23.7307,26.503,23.8103,26.0207,24.817,
       26.7509,25.6751,28.0619,25.5788]
list2=[91.25,91.74,95.016,92.850,92.820,92.717,90.5536,90.4125,89.9743,93.9368,
       94.2037,94.1086,95.326,91.8933]
list=['AGARTALA',"GUWAHATI","DIBRUGARH","SILCHAR","TEZPUR","AIZAWL","BONGAIGAON","DHAKA","DHUBURI","IMPHAL",
      "JORHAT","KOHIMA",'PASIGHAT','SHILLONG']
for sel_lat,sel_lon,loc in zip(list1,list2,list):
        #sel_lat =26.1
    #sel_lon = 91.74
    a = abs(lats-sel_lat)+abs(lons-sel_lon)
    i,j = np.unravel_index(a.argmin(), a.shape)   #i=lat j=lon
    print(i)
    print(j)
    precp = ds.variables["AODANA"][:,i,j]
    precp = DataFrame(precp)
    print(precp)
    precp.to_csv(loc+".csv")
'''
#
#wrf = Dataset(r"G:\WRF_Chem_Output\202new\March\wrfout_d01_2018-03-10_00%3A00%3A00")
#lat = getvar(wrf,"lat",timeidx=1)
#lon = getvar(wrf,"lon",timeidx=1)


#lat_max= np.max(lat)
lat_max= 40
lat_min= np.min(lat)
#lon_max= np.max(lon)
lon_max= 103.5
lon_min= 0



a = abs(lat-lat_min)+abs(lon-lon_min)
i,j = np.unravel_index(a.argmin(), a.shape)
b = abs(lat-lat_max)+abs(lon-lon_max)
k,l = np.unravel_index(b.argmin(), b.shape)


precp = precp[i:k,j:l] ###unit =mm/hr,change unit to mm, 720 for APRIL
lat = ds.variables["lat"][i:k]
lon = ds.variables["lon"][j:l]
#precp = diff.transpose("lat","lon")

font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
plt.figure(figsize=(10,6))
RAINNC766_1 = getvar(wrf, "RAINNC", timeidx=0)
cart_proj = get_cartopy(RAINNC766_1)
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
ax.set_xlim(cartopy_xlim(RAINNC766_1))
ax.set_ylim(cartopy_ylim(RAINNC766_1))
v = np.linspace(0,1,20)
plt.contourf(lon,lat,precp,v,cmap="nipy_spectral",extend='max',orientation="vertical",
             transform=cartopy.crs.PlateCarree())
plt.colorbar(ticks=v,shrink=0.96)
plt.title("AOD_MERRA2_April11-15")
plt.savefig(r'F:\AOD_ANALYSIS\MERRA2_April11-15.png',dpi=300,bbox_inches='tight')
plt.show()
'''
