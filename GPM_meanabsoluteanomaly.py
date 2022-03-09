


##IT CALCULATES THE OVERALL MEAN ABSOLUTE ANOMALY FOR EACH MONTH
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
ds = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\April_2001-2020_monthlyprecip.nc")
ds1 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20010401-S000000-E235959.04.V06B.HDF5.nc4")
ds2 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20020401-S000000-E235959.04.V06B.HDF5.nc4")
ds3 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20030401-S000000-E235959.04.V06B.HDF5.nc4")
ds4 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20040401-S000000-E235959.04.V06B.HDF5.nc4")
ds5 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20050401-S000000-E235959.04.V06B.HDF5.nc4")
ds6 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20060401-S000000-E235959.04.V06B.HDF5.nc4")
ds7 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20070401-S000000-E235959.04.V06B.HDF5.nc4")
ds8 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20080401-S000000-E235959.04.V06B.HDF5.nc4")
ds9 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20090401-S000000-E235959.04.V06B.HDF5.nc4")
ds10 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20100401-S000000-E235959.04.V06B.HDF5.nc4")
ds11 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20110401-S000000-E235959.04.V06B.HDF5.nc4")
ds12 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20120401-S000000-E235959.04.V06B.HDF5.nc4")
ds13 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20130401-S000000-E235959.04.V06B.HDF5.nc4")
ds14 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20140401-S000000-E235959.04.V06B.HDF5.nc4")
ds15 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20150401-S000000-E235959.04.V06B.HDF5.nc4")
ds16 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20160401-S000000-E235959.04.V06B.HDF5.nc4")
ds17 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20170401-S000000-E235959.04.V06B.HDF5.nc4")
ds18 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20180401-S000000-E235959.04.V06B.HDF5.nc4")
ds19 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20190401-S000000-E235959.04.V06B.HDF5.nc4")
ds20 = open_mfdataset(r"I:\GPM\April_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20200401-S000000-E235959.04.V06B.HDF5.nc4")
'''
ds = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\March_2001-2020_monthlyprecip.nc")
ds1 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20010301-S000000-E235959.03.V06B.HDF5.nc4")
ds2 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20020301-S000000-E235959.03.V06B.HDF5.nc4")
ds3 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20030301-S000000-E235959.03.V06B.HDF5.nc4")
ds4 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20040301-S000000-E235959.03.V06B.HDF5.nc4")
ds5 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20050301-S000000-E235959.03.V06B.HDF5.nc4")
ds6 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20060301-S000000-E235959.03.V06B.HDF5.nc4")
ds7 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20070301-S000000-E235959.03.V06B.HDF5.nc4")
ds8 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20080301-S000000-E235959.03.V06B.HDF5.nc4")
ds9 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20090301-S000000-E235959.03.V06B.HDF5.nc4")
ds10 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20100301-S000000-E235959.03.V06B.HDF5.nc4")
ds11 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20110301-S000000-E235959.03.V06B.HDF5.nc4")
ds12 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20120301-S000000-E235959.03.V06B.HDF5.nc4")
ds13 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20130301-S000000-E235959.03.V06B.HDF5.nc4")
ds14 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20140301-S000000-E235959.03.V06B.HDF5.nc4")
ds15 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20150301-S000000-E235959.03.V06B.HDF5.nc4")
ds16 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20160301-S000000-E235959.03.V06B.HDF5.nc4")
ds17 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20170301-S000000-E235959.03.V06B.HDF5.nc4")
ds18 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20180301-S000000-E235959.03.V06B.HDF5.nc4")
ds19 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20190301-S000000-E235959.03.V06B.HDF5.nc4")
ds20 = open_mfdataset(r"I:\GPM\March_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20200301-S000000-E235959.03.V06B.HDF5.nc4")
'''
'''
ds = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\May_2001-2020_monthlyprecip.nc")
ds1 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20010501-S000000-E235959.05.V06B.HDF5.nc4")
ds2 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20020501-S000000-E235959.05.V06B.HDF5.nc4")
ds3 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20030501-S000000-E235959.05.V06B.HDF5.nc4")
ds4 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20040501-S000000-E235959.05.V06B.HDF5.nc4")
ds5 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20050501-S000000-E235959.05.V06B.HDF5.nc4")
ds6 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20060501-S000000-E235959.05.V06B.HDF5.nc4")
ds7 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20070501-S000000-E235959.05.V06B.HDF5.nc4")
ds8 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20080501-S000000-E235959.05.V06B.HDF5.nc4")
ds9 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20090501-S000000-E235959.05.V06B.HDF5.nc4")
ds10 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20100501-S000000-E235959.05.V06B.HDF5.nc4")
ds11 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20110501-S000000-E235959.05.V06B.HDF5.nc4")
ds12 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20120501-S000000-E235959.05.V06B.HDF5.nc4")
ds13 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20130501-S000000-E235959.05.V06B.HDF5.nc4")
ds14 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20140501-S000000-E235959.05.V06B.HDF5.nc4")
ds15 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20150501-S000000-E235959.05.V06B.HDF5.nc4")
ds16 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20160501-S000000-E235959.05.V06B.HDF5.nc4")
ds17 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20170501-S000000-E235959.05.V06B.HDF5.nc4")
ds18 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20180501-S000000-E235959.05.V06B.HDF5.nc4")
ds19 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20190501-S000000-E235959.05.V06B.HDF5.nc4")
ds20 = open_mfdataset(r"I:\GPM\May_2001-2020_monthly\3B-MO.MS.MRG.3IMERG.20200501-S000000-E235959.05.V06B.HDF5.nc4")
'''
precp = ds.variables["precipitation"][:,:,:]
precp = precp.mean(axis=0)  ###calculate 2001-2020 mean
precp1 = ds1.variables["precipitation"][:,:,:]
precp1 = precp1.squeeze(dim='time') ##change to 2D and remove 'time' dimension
precp2 = ds2.variables["precipitation"][:,:,:]
precp2 = precp2.squeeze(dim='time')
precp3 = ds3.variables["precipitation"][:,:,:]
precp3 = precp3.squeeze(dim='time')
precp4 = ds4.variables["precipitation"][:,:,:]
precp4 = precp4.squeeze(dim='time')
precp5 = ds5.variables["precipitation"][:,:,:]
precp5 = precp5.squeeze(dim='time')
precp6 = ds6.variables["precipitation"][:,:,:]
precp6 = precp6.squeeze(dim='time')
precp7 = ds7.variables["precipitation"][:,:,:]
precp7 = precp7.squeeze(dim='time')
precp8 = ds8.variables["precipitation"][:,:,:]
precp8 = precp8.squeeze(dim='time')
precp9 = ds9.variables["precipitation"][:,:,:]
precp9 = precp9.squeeze(dim='time')
precp10 = ds10.variables["precipitation"][:,:,:]
precp10 = precp10.squeeze(dim='time')
precp11 = ds11.variables["precipitation"][:,:,:]
precp11 = precp11.squeeze(dim='time')
precp12 = ds12.variables["precipitation"][:,:,:]
precp12 = precp12.squeeze(dim='time')
precp13 = ds13.variables["precipitation"][:,:,:]
precp13 = precp13.squeeze(dim='time')
precp14 = ds14.variables["precipitation"][:,:,:]
precp14 = precp14.squeeze(dim='time')
precp15 = ds15.variables["precipitation"][:,:,:]
precp15 = precp15.squeeze(dim='time')
precp16 = ds16.variables["precipitation"][:,:,:]
precp16 = precp16.squeeze(dim='time')
precp17 = ds17.variables["precipitation"][:,:,:]
precp17 = precp17.squeeze(dim='time')
precp18 = ds18.variables["precipitation"][:,:,:]
precp18 = precp18.squeeze(dim='time')
precp19 = ds19.variables["precipitation"][:,:,:]
precp19 = precp19.squeeze(dim='time')
precp20 = ds20.variables["precipitation"][:,:,:]
precp20 = precp20.squeeze(dim='time')

wrf = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
lat = getvar(wrf,"lat",timeidx=1)
lon = getvar(wrf,"lon",timeidx=1)


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


diff = ((abs(precp1[j:l,i:k] - precp[j:l,i:k])+abs(precp2[j:l,i:k] - precp[j:l,i:k])+abs(precp3[j:l,i:k] - precp[j:l,i:k])+abs(precp4[j:l,i:k] - precp[j:l,i:k])+abs(precp5[j:l,i:k] - precp[j:l,i:k])+
        abs(precp6[j:l,i:k] - precp[j:l,i:k])+abs(precp7[j:l,i:k] - precp[j:l,i:k])+abs(precp8[j:l,i:k] - precp[j:l,i:k])+abs(precp9[j:l,i:k] - precp[j:l,i:k])+
        abs(precp10[j:l,i:k] - precp[j:l,i:k])+abs(precp11[j:l,i:k] - precp[j:l,i:k])+abs(precp12[j:l,i:k] - precp[j:l,i:k])+abs(precp13[j:l,i:k] - precp[j:l,i:k])+
        abs(precp14[j:l,i:k] - precp[j:l,i:k])+abs(precp15[j:l,i:k] - precp[j:l,i:k])+abs(precp16[j:l,i:k] - precp[j:l,i:k])+abs(precp17[j:l,i:k] - precp[j:l,i:k])+
        abs(precp18[j:l,i:k] - precp[j:l,i:k])+abs(precp19[j:l,i:k] - precp[j:l,i:k])+abs(precp20[j:l,i:k] - precp[j:l,i:k]))*744)/20 ###unit =mm/hr,change unit to mm, 720 for APRIL
lat = ds.variables["lat"][i:k]
lon = ds.variables["lon"][j:l]
precp = diff.transpose("lat","lon")

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
plt.title("mean absolute rainfall anomaly from 2001-2020")
plt.savefig('F:\GPM_ANALYSIS\April\GPM_April_mean_absolute_anomaly.png',dpi=300,bbox_inches='tight')
plt.show()