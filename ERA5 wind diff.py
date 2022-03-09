from xarray import open_mfdataset
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from wrf import get_cartopy,getvar
import numpy as np
from xarray import open_mfdataset
from nco import *
import xarray as xr
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
import wrf as w
import math

file = open_mfdataset(r"I:\observation_data\APRIL\DATA\Wind_data\4.ERA5-2018APR_U_V_3hr.nc")
wrf = Dataset(r"G:\WRF_Outputs\Monsoon\ACM2\wrfout_d03.nc")
lat = getvar(wrf,"lat",timeidx=0)
lon = getvar(wrf,"lon",timeidx=0)
C = getvar(wrf,"RAINC",timeidx=0)
lat_max= np.max(lat)
lat_min= np.min(lat)
lon_max= np.max(lon)
lon_min= np.min(lon)
#lat_max= np.max(lat)
#lat_min= np.min(lat)
#lon_max= np.max(lon)
#lon_min= np.min(lon)

lat1 = file.variables["latitude"][:] ##lat starts from 35......21 (WRONG)
lat1 = lat1[::-1] ##reversed the index
#print(lat1)
lon1 = file.variables["longitude"][:]

a = abs(lat1-lat_min)+abs(lon1-lon_min)
i,j = np.unravel_index(a.argmin(), a.shape)
#print(a)
b = abs(lat1-lat_max)+abs(lon1-lon_max)
#print(b)
k,l = np.unravel_index(b.argmin(), b.shape)

print(i)
#print(j)
print(k)
#print(l)
i=22
k=56
lat2 = file.variables["latitude"][i:k]
lon2 = file.variables["longitude"][j:l]

u = file.variables["u"][:,88,i:k,j:l]
v = file.variables["v"][:,88,i:k,j:l]

u = u.mean(axis=0)
v = v.mean(axis=0)
DIR= 180+(180/3.14)*(np.arctan2(v,u))
print(DIR)
#lats,lons =w.latlon_coords(SPD)
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
ax.add_feature(cartopy.feature.BORDERS, linestyle='--')
#ax.set_xlim(cartopy_xlim(RAINNC766_1))
#ax.set_ylim(cartopy_ylim(RAINNC766_1))
#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),30) ##COLORBAR LIMITS
v = np.linspace(0,360,16)
#plt.figure(figsize=(8,6))
plt.contourf(lon2, lat2, DIR,v, cmap=plt.get_cmap("jet"),antialiased=False,
                 transform=cartopy.crs.PlateCarree())
plt.colorbar(ax=ax, shrink=.90,ticks=v,extend="max")
#plt.title("(d) MYNN3-April",weight='normal',fontdict=font)
#plt.title("(f) YSU-April")
#plt.savefig("QFX_YSU_April.png",dpi=300,bbox_inches='tight')
plt.show()










'''


#u10 = file.variables["u10"][:,:,:]
#v10 = file.variables["v10"][:,:,:]
u10 = file.variables["u"][:,88,:,:]
v10 = file.variables["v"][:,88,:,:]
print(u10)
#lats= file.variables["latitude"][:]
#lons= file.variables["longitude"][:]
lats= file.variables["latitude"][:]
lons= file.variables["longitude"][:]
sel_lat =26.1
sel_lon = 91.583
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)

#print(i)-->latitude
#print(j)-->longitude
u = u10[:,i,j]
v = v10[:,i,j]
#u = np.asarray(u)
#v = np.asarray(v)
u=pd.DataFrame(u)
v=pd.DataFrame(v)
print(u)
print(v)
u.to_csv("ERA5_Model_U10_apr.csv")
v.to_csv("ERA5_Model_V10_apr.csv")

import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
import pandas as pd
#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")

#wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset("G:\WRF_Outputs\Monsoon\ACM2\wrf_trad_d03.nc")
wrf_file1 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrfout_d03.nc")
wrf_file2 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYJ\wrfout_d03.nc")
wrf_file3 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrfout_d03.nc")
wrf_file4 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_YSU\wrfout_d03.nc")
wrf_file5 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_QNSE\wrfout_d03.nc")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
time = ALL_TIMES
#time = 60
#level = range(0,25,1)
C = getvar(wrf_file1, "RAINC", timeidx=1)
CLD = getvar(wrf_file, "DIR", timeidx=time)
#CLD1 = getvar(wrf_file, "QFX", timeidx=time)
#ds = pd.DataFrame(CLD[24:744,156,136])
#ds.to_csv("pblh_ghy_acm2.csv")
#print(ds)
if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        #CLD1 = CLD1.mean("Time")
#print(CLD[156,136])
 #lats, lons = w.latlon_coords(CLD)

CLD_DIFF = CLD[0,:,:]#- CLD1)*1000
lats,lons =w.latlon_coords(CLD)
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
ax.add_feature(cartopy.feature.BORDERS, linestyle='--')
#ax.set_xlim(cartopy_xlim(RAINNC766_1))
#ax.set_ylim(cartopy_ylim(RAINNC766_1))
#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),30) ##COLORBAR LIMITS
v = np.linspace(np.min(CLD),np.max(CLD),20)
#plt.figure(figsize=(8,6))
plt.contourf(lons, lats, CLD_DIFF,v, cmap=plt.get_cmap("jet"),antialiased=False,
                 transform=cartopy.crs.PlateCarree())
plt.colorbar(ax=ax, shrink=.90,ticks=v,extend="max")
#plt.title("(d) MYNN3-April",weight='normal',fontdict=font)
#plt.title("(f) YSU-April")
#plt.savefig("QFX_YSU_April.png",dpi=300,bbox_inches='tight')
plt.show()
'''
