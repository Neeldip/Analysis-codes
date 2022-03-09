

##CALCULATES combined absolute anomaly for each year considering march,april,may
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
#ds = open_mfdataset(r"I:\GPM\2018\April_2001-2020_monthly\*.nc4",combine='by_coords',concat_dim='time')
#precp = ds.variables["precipitation"][:,:,:]
#ds.to_netcdf(r'I:\GPM\2018\April_2001-2020_monthly\April_2001-2020_monthlyprecip.nc')
mean1 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\March_2001-2020_monthlyprecip.nc")
mean2 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\April_2001-2020_monthlyprecip.nc")
mean3 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\May_2001-2020_monthlyprecip.nc")

###only change the year
mar = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20190301-S000000-E235959.03.V06B.HDF5.nc4")
apr = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20190401-S000000-E235959.04.V06B.HDF5.nc4")
may =   open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20200501-S000000-E235959.05.V06B.HDF5.nc4")

precp1 = mean1.variables["precipitation"][:,:,:]
precp1 = precp1.mean(axis=0)  ###calculate 2001-2020 mean,coverts 3d to 2d
precp2 = mean2.variables["precipitation"][:,:,:]
precp2 = precp2.mean(axis=0)
precp3 = mean3.variables["precipitation"][:,:,:]
precp3 = precp3.mean(axis=0)




mar = mar.variables["precipitation"][:,:,:]
mar = mar.squeeze(dim='time') ##change to 2D and remove 'time' dimension
apr = apr.variables["precipitation"][:,:,:]
apr = apr.squeeze(dim='time')
may = may.variables["precipitation"][:,:,:]
may = may.squeeze(dim='time')

wrf = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
lat = getvar(wrf,"lat",timeidx=1)
lon = getvar(wrf,"lon",timeidx=1)
#lat_max= np.max(lat)
lat_max= 32
lat_min= np.min(lat)
#lon_max= np.max(lon)
lon_max= 120
lon_min= np.min(lon)

lat = mean1.variables["lat"][:]
lon = mean1.variables["lon"][:]

a = abs(lat-lat_min)+abs(lon-lon_min)
i,j = np.unravel_index(a.argmin(), a.shape)
b = abs(lat-lat_max)+abs(lon-lon_max)
k,l = np.unravel_index(b.argmin(), b.shape)

##MAA=sum of mean absolute anomaly
SAA = (abs(mar[j:l,i:k] - precp1[j:l,i:k])*744+abs(apr[j:l,i:k] - precp2[j:l,i:k])*720)#+abs(may[j:l,i:k] - precp3[j:l,i:k])*744) ###unit =mm/hr,change unit to mm, 720 for APRIL

lat = mean1.variables["lat"][i:k]
lon = mean1.variables["lon"][j:l]

SAA = SAA.transpose("lat","lon")

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
v = np.linspace(np.min(SAA),np.max(SAA),20)
plt.contourf(lon,lat,SAA,v,cmap="nipy_spectral",orientation="vertical",
             transform=cartopy.crs.PlateCarree())
plt.colorbar(ticks=v,shrink=0.93)
plt.title("2019-sum of absolute anomaly of Mar,Apr")
plt.savefig('F:\GPM_ANALYSIS\GPM_2019_SAA_MAR_APR.png',dpi=300,bbox_inches='tight')
plt.show()