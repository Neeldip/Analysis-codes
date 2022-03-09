from xarray import open_mfdataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import sys
#np.set_printoptions(threshold=sys.maxsize)
#sys.stdout = open(r'F:\IMD\z_statistics_apr.txt', 'w')   #<----------------
import pymannkendall as pk

#*****************concatenate************************
#ncrcat * IMD.nc

ds = open_mfdataset(r"I:\observation_data\APRIL\DATA\Rainfall\IMD_gridded_datasets\RAW data\NETCDF\New folder\Working\apr_2001-2020.nc") #<----------------
'''
for i in range(0,129,1):
    for j in range(0,135,1):
        RAINFALL = ds.variables["RAINFALL"][:, i, j]
        RAINFALL = np.asarray(RAINFALL)

        if np.isnan(RAINFALL).any() == True:
            z = np.nan
            print(z)
        else:
            z=pk.original_test(RAINFALL,alpha=0.05)
            #z = pk.hamed_rao_modification_test(RAINFALL, alpha=0.05)
            print(z)
sys.stdout.close()
'''

######################crete netcdf.py########################
import netCDF4 as nc


fn = r'F:\IMD\Trend analysis\mann kendall\slope\2011-2020\z_apr_2001-2020.nc'   #<----------------
ds1 = nc.Dataset(fn, 'w', format='NETCDF4')
ds1.set_fill_on()

lat = ds1.createDimension('lat', 129)
lon = ds1.createDimension('lon', 135)

lats = ds1.createVariable('lat', 'f4', ('lat',))
lons = ds1.createVariable('lon', 'f4', ('lon',))
value = ds1.createVariable('value', 'f4', ('lat', 'lon',),fill_value=np.nan)
value.units = 'MK test slope'
lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = np.arange(6.5, 38.75, 0.25)
lons[:] = np.arange(66.5, 100.25, 0.25)

#df = pd.read_csv(r"F:\IMD\z_statistics_apr.txt",header=None)  #<----------------
#df = np.asarray(df)
#print(np.nanmax(df))

#df = np.reshape(df,(129,135))

#for i in range(0,129,1):
#    for j in range(0,135,1):
#        if np.isnan(df[i,j]).any() == True:
#           value[i,j] = np.nan
#
#        else:
#           value[i,j] = df[i,j]

for i in range(0,129,1):
    for j in range(0,135,1):
        Rain = ds.variables["RAINFALL"][:,i,j]
        if np.isnan(Rain).any() == True:
           value[i,j] = np.nan

        else:
            RAIN =[]
            #for timestep in range(0,620,31): ##for months with 31 days
            for timestep in range(0,600,30):  ##for months with 30 days
            #for timestep in range(0,560,28):  ##for aprruary
                #rain = ds.variables["RAINFALL"][timestep:timestep + 31, i, j]
                rain = ds.variables["RAINFALL"][timestep:timestep + 30, i, j]
                #rain = ds.variables["RAINFALL"][timestep:timestep + 28, i, j]
                IMD = np.sum(rain)
                RAIN.append(IMD)
                if len(RAIN) == 20:
                    rain = pk.original_test(RAIN)
                    #rain = np.std(RAIN)
                    print(i)
                    print(j)
                    value[i,j] = rain


ds1.close()

###########################IMD plot.py############################
from netCDF4 import Dataset
from wrf import get_cartopy,getvar
from xarray import open_mfdataset
from nco import *
import xarray as xr
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim
import h5py
ds = open_mfdataset(r"F:\IMD\Trend analysis\mann kendall\slope\2011-2020\z_apr_2001-2020.nc") #<----------------

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
#v = np.linspace(np.min(precp),np.max(precp),20)
v = [-1.96,0,1.96]
#plt.contourf(lon,lat,precp,v,cmap="nipy_spectral",orientation="vertical",
#             transform=cartopy.crs.PlateCarree())
#plt.pcolormesh(lon,lat,precp,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
plt.pcolormesh(lon,lat,precp,vmin=-1.96,vmax=1.96,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
plt.colorbar(ticks=v,shrink=0.93)
plt.title("April")  #<----------------
plt.savefig(r'F:\IMD\Trend analysis\mann kendall\slope\2011-2020\z_apr_2001-2020.png',dpi=300,bbox_inches='tight')  #<----------------
plt.show()