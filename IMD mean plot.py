from xarray import open_mfdataset
import numpy as np
import warnings

warnings.filterwarnings("ignore")


import netCDF4 as nc

fn = r'F:\IMD\IMD_sum_10-19_Apr_2018.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')

ds.set_fill_on()

lat = ds.createDimension('lat', 129)
lon = ds.createDimension('lon', 135)

lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('Rainfall', 'f4', ('lat', 'lon',),fill_value=np.nan)
value.units = 'mm/day'
lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = np.arange(6.5, 38.75, 0.25)
lons[:] = np.arange(66.5, 100.25, 0.25)

ds1 = open_mfdataset(r"I:\observation_data\APRIL\DATA\Rainfall\IMD_gridded_datasets\RAW data\NETCDF\New folder\Working\apr\New folder\2018-Apr.nc") #<----------------

for i in range(0,129,1):
    for j in range(0,135,1):
        RAINFALL = ds1.variables["RAINFALL"][9:19, i, j]
        RAINFALL = np.asarray(RAINFALL)
        print(np.shape(RAINFALL))
        if np.isnan(RAINFALL).any() == True:
            z = np.nan
            #print(z)
        else:
            z=np.sum(RAINFALL)
            #z = np.mean(RAINFALL)
            #print(z)
        value[i, j] = z

ds.close()

###########################IMD plot.py############################
from netCDF4 import Dataset
from wrf import get_cartopy, getvar
from xarray import open_mfdataset

import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from wrf import get_cartopy, cartopy_xlim, cartopy_ylim

ds = open_mfdataset(r"F:\IMD\IMD_sum_10-19_Apr_2018.nc")  # <----------------

wrf =  nc.Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")

lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]
precp = ds.variables["Rainfall"][:]

font = {'family': 'serif',
        'color': 'black',
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
# ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.COASTLINE)
# ax.add_feature(cartopy.feature.LAKES)
# ax.add_feature(cartopy.feature.RIVERS)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
ax.set_xlim(cartopy_xlim(RAINNC766_1))
ax.set_ylim(cartopy_ylim(RAINNC766_1))
v = np.linspace(np.min(precp), np.max(precp), 11)
# plt.contourf(lon,lat,precp,v,cmap="nipy_spectral",orientation="vertical",
#             transform=cartopy.crs.PlateCarree(),extend='max')
plt.contourf(lon,lat,precp,v,cmap="jet",orientation="vertical",
             transform=cartopy.crs.PlateCarree())
#plt.pcolormesh(lon, lat, precp, vmin=0, vmax=20, transform=cartopy.crs.PlateCarree(), cmap="jet")
plt.colorbar(ticks=v, shrink=0.93)
#plt.title("December")  # <----------------
plt.savefig(r'F:\IMD\IMD_APR_10-19.png', dpi=600, bbox_inches='tight')  # <----------------
plt.show()
