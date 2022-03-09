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
#*****************concatenate************************
ds = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\*.nc4",combine='by_coords',concat_dim='time')
#precp = ds.variables["precipitation"][:,:,:]
ds.to_netcdf(r'I:\GPM\May_2001-2020_monthly\May_2001-2020_monthlyprecip.nc')

ds = open_mfdataset(r'I:\GPM\May_2001-2020_monthly\May_2001-2020_monthlyprecip.nc')
precp = ds.variables["precipitation"][:,:,:]
precp = precp.mean(axis=0)
print(precp)
#wrf = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
lat = getvar(wrf,"lat",timeidx=1)
lon = getvar(wrf,"lon",timeidx=1)
#lat_max= np.max(lat)
#lat_min= np.min(lat)
#lon_max= np.max(lon)
#lon_min= np.min(lon)

#lat_max= np.max(lat)
lat_max= 32
lat_min= np.min(lat)
#lon_max= np.max(lon)
lon_max= 120
lon_min= np.min(lon)

lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]

a = abs(lat-lat_min)+abs(lon-lon_min)
i,j = np.unravel_index(a.argmin(), a.shape)
b = abs(lat-lat_max)+abs(lon-lon_max)
k,l = np.unravel_index(b.argmin(), b.shape)


#precp = ds.variables["precipitation"][:]
#precp = precp.squeeze(dim='time')
#precp = precp.sum(dim='time')
precp = precp[j:l,i:k]*744 ###unit =mm/hr,change unit to mm, 720 for APRIL,744 for march and may
lat = ds.variables["lat"][i:k]
lon = ds.variables["lon"][j:l]
#precp = precp*720 ###change unit to mm
#lat = ds.variables["lat"][:]
#lon = ds.variables["lon"][:]
#print(precp)
#print(np.shape(precp))
#precp = np.squeeze(precp)



precp = precp.transpose("lat","lon")
#print(precp)
font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
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
v = np.linspace(np.min(precp),np.max(precp),20)
plt.contourf(lon,lat,precp,v,cmap="nipy_spectral",orientation="vertical",
             transform=cartopy.crs.PlateCarree())
plt.colorbar(ticks=v,shrink=0.93)
plt.title("2001-2020 mean")
plt.savefig('F:\GPM_ANALYSIS\May\GPM_MARCH_2001-2020_mean.png',dpi=300,bbox_inches='tight')
plt.show()