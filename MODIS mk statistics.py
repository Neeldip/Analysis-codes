import numpy as np
from pyhdf.SD import SD, SDC
import netCDF4 as nc

#https://hdfeos.org/software/pyhdf.php
#https://hdfeos.org/examples/c_grid_lonlat.php

import pymannkendall as pk

fn = 'F:\AOD_ANALYSIS\statistics\mann kendall\MODIS_z_statistics_DEC.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')
ds.set_fill_on()

#lat = ds.createDimension('lat', 35)
#lon = ds.createDimension('lon', 44)
lat = ds.createDimension('lat', 18)
lon = ds.createDimension('lon', 23)
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('value', 'f4', ('lat', 'lon',),fill_value=np.nan)
value.units = 'Unknown'
lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = np.arange(34.5, 16.5, -1)
lons[:] = np.arange(81.5, 104.5, 1)

#FILE_NAME = r'I:\MYD08_M3\2018\DEC.hdf'

list=[r'I:\MOD08_M3\2001\DEC.hdf',
      r'I:\MOD08_M3\2002\DEC.hdf',
      r'I:\MOD08_M3\2003\DEC.hdf',
      r'I:\MOD08_M3\2004\DEC.hdf',
      r'I:\MOD08_M3\2005\DEC.hdf',
      r'I:\MOD08_M3\2006\DEC.hdf',
      r'I:\MOD08_M3\2007\DEC.hdf',
      r'I:\MOD08_M3\2008\DEC.hdf',
      r'I:\MOD08_M3\2009\DEC.hdf',
      r'I:\MOD08_M3\2010\DEC.hdf',
      r'I:\MOD08_M3\2011\DEC.hdf',
      r'I:\MOD08_M3\2012\DEC.hdf',
      r'I:\MOD08_M3\2013\DEC.hdf',
      r'I:\MOD08_M3\2014\DEC.hdf',
      r'I:\MOD08_M3\2015\DEC.hdf',
      r'I:\MOD08_M3\2016\DEC.hdf',
      r'I:\MOD08_M3\2017\DEC.hdf',
      r'I:\MOD08_M3\2018\DEC.hdf',
      r'I:\MOD08_M3\2019\DEC.hdf',
      r'I:\MOD08_M3\2020\DEC.hdf',]

#for i in range(49,84,1):
for i in range(55, 73, 1):
    #for j in range(240,284,1):
    for j in range(261, 284, 1):
        data_array = []
        for FILE_NAME in list:
            hdf = SD(FILE_NAME, SDC.READ)
            # Read dataset with PANDAS
            #DATAFIELD_NAME='Cloud_Optical_Thickness_Liquid_Mean_Mean'
            DATAFIELD_NAME='AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean'
            #DATAFIELD_NAME='Cloud_Fraction_Mean_Mean'
            #DATAFIELD_NAME='Cloud_Effective_Radius_Liquid_Mean_Mean'

            AOD = hdf.select(DATAFIELD_NAME)
            data = AOD[i,j]
            if data == -9999:
                data = np.nan
            else:
                data = data/1000
                # ds = ds*0.009999999776482582 #COD & droplet radius
                # ds = ds*9.999999747378752/100000 #cloud fraction

            #print(data)
            data_array.append(data)
            #print(data_array)
            if len(data_array) == 20:
                if np.count_nonzero(np.isnan(data_array)) > 10:  #count # of nan values
                    z = np.nan                                   #if nan values > 10,
                else:
                    #z=pk.hamed_rao_modification_test(data_array)
                    z = pk.original_test(data_array)
                    #z = np.asarray(z)
                    #z = z.astype('float')
                print(z)
                print(i)
                print(j)
                value[i-55, j-261] = z
                #data_array.clear()
ds.close()

########################IMD plot.py############################
from netCDF4 import Dataset
from wrf import getvar
import numpy as np
from xarray import open_mfdataset
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim

ds = open_mfdataset(r"F:\AOD_ANALYSIS\statistics\mann kendall\MODIS_z_statistics_DEC.nc") #<----------------

wrf = Dataset(r"G:\WRF_Outputs\Monsoon\ACM2\wrfout_d02_2018-07-01_00%3A00%3A00")

lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]
precp = ds.variables["value"][:]


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
#v = np.linspace(-5,5,20)
v = [-1.96,0,1.96]
#plt.contourf(lon,lat,precp,v,cmap="nipy_spectral",orientation="vertical",
#             transform=cartopy.crs.PlateCarree())
plt.pcolormesh(lon,lat,precp,vmin=-1.96,vmax=1.96,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
plt.colorbar(ticks=v,shrink=0.93)
plt.title("December")  #<----------------
plt.savefig(r'F:\AOD_ANALYSIS\statistics\mann kendall\MODIS_z_statistics_DEC.png',dpi=300,bbox_inches='tight')  #<----------------
plt.show()
