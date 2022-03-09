from netCDF4 import Dataset
from wrf import get_cartopy,getvar
import numpy as np
from xarray import open_mfdataset
from nco import *
import xarray as xr
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim,ll_to_xy
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from mpl_toolkits.basemap import *
from xarray import open_dataset
#ds = open_mfdataset(r"I:\GPM\2018\April_4-9\*.nc4",combine='by_coords',concat_dim='time')
ds = open_mfdataset(r"F:\WRF-CHEM ANALYSIS\rainfall\LOCATION ANALYSIS\26.67129 89.4322\GPM\8.30-9.30\*.nc4",combine='by_coords',concat_dim='time')
#wrf = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
wrf = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
lat = getvar(wrf,"lat",timeidx=1)
lon = getvar(wrf,"lon",timeidx=1)
lat_max= np.max(lat)
lat_min= np.min(lat)
lon_max= np.max(lon)
lon_min= np.min(lon)
#lat_max= np.max(lat)
#lat_min= np.min(lat)
#lon_max= np.max(lon)
#lon_min= np.min(lon)

lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]

a = abs(lat-lat_min)+abs(lon-lon_min)
#print(a)
i,j = np.unravel_index(a.argmin(), a.shape)
b = abs(lat-lat_max)+abs(lon-lon_max)
k,l = np.unravel_index(b.argmin(), b.shape)
#print(b)

precp = ds.variables["precipitationCal"][:]
precp = (precp.sum(dim='time'))/2
precp1 = precp[j:l+1,i:k+1]
lat = ds.variables["lat"][i:k+1]
lon = ds.variables["lon"][j:l+1]

precp = precp1.transpose("lat","lon")
#lats,lons =latlon_coords(RAINNC766_1)
ds2 = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
RAINC231 = getvar(ds2, "RAINC", timeidx=0)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
cart_proj = get_cartopy(RAINC231)
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
ax.set_xlim(cartopy_xlim(RAINC231))
ax.set_ylim(cartopy_ylim(RAINC231))
#ax.gridlines(color="black", linestyle="dotted")
v = np.linspace(np.min(precp),np.max(precp),21) ##COLORBAR LIMITS
#v = np.linspace(-250,250,16) ##RAINNC LIMITS
#v = np.linspace(-100,100,16) ##RAINC LIMITS
#v = np.linspace(-100,100,21) ##TOTAL LIMITS
plt.contourf(lon, lat, precp,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
plt.colorbar(ax=ax, shrink=.81,ticks=v)
#plt.title("RAINCcontribution_change_aerfeedback - no_bc_absorbtion rainfall difference")
#plt.figure(figsize=(8,6))
plt.savefig(r'F:\WRF-CHEM ANALYSIS\rainfall\LOCATION ANALYSIS\26.67129 89.4322\GPM\GPM_8.30-9.30.png',dpi=600,bbox_inches='tight')
plt.show()





#v = np.linspace(np.min(precp),100,30)
#plt.contourf(lon,lat,precp,v,cmap="nipy_spectral",orientation="vertical")
#plt.colorbar()
#plt.savefig('GPM_10_19APRIL.png',dpi=500)
#plt.show()