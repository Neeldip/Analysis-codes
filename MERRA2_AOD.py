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


#ds = open_mfdataset(r"I:\MERRA2_AOD\March\2018\5-10\*.nc",combine='by_coords',concat_dim='time')
#ds = open_mfdataset(r"I:\MERRA2_AOD\2001-2020\Monthwise\11.Nov\*.nc",combine='by_coords',concat_dim='time')
ds = open_mfdataset(r"I:\MERRA2_AOD\April\2018\10-19\*.nc",combine='by_coords',concat_dim='time')
ds1 =ds.mean(dim="time")

precp = ds1.variables["AODANA"][:,:]

wrf = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
lat = getvar(wrf,"lat",timeidx=1)
lon = getvar(wrf,"lon",timeidx=1)


lat_max= np.max(lat)
#lat_max= 40
lat_min= np.min(lat)
lon_max= np.max(lon)
#lon_max= 103.5
lon_min= 70

lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]

a = abs(lat-lat_min)+abs(lon-lon_min)
i,j = np.unravel_index(a.argmin(), a.shape)
b = abs(lat-lat_max)+abs(lon-lon_max)
k,l = np.unravel_index(b.argmin(), b.shape)


precp = precp[i-1:k+2,j:l+2] ###unit =mm/hr,change unit to mm, 720 for APRIL
lat = ds.variables["lat"][i-1:k+2]
lon = ds.variables["lon"][j:l+2]
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
plt.contourf(lon,lat,precp,v,cmap="jet",orientation="vertical",
             transform=cartopy.crs.PlateCarree())
plt.colorbar(ticks=v,shrink=0.96)
plt.title("(c) MERRA-2")
plt.savefig(r'F:\AOD_ANALYSIS\MERRA2_April_10-19.png',dpi=600,bbox_inches='tight')
plt.show()