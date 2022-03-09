
#1.FIRST DO REGRIDDING OF MERRA2 DATA USING MERRA2 REGRIDDING.PY
#2.CDO REMAPBIL,TEMPLATE.NC INFILE OUTFILE

from xarray import open_dataset,open_mfdataset
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
from wrf import get_cartopy,getvar
from xarray import open_mfdataset
#from nco import *
import xarray as xr
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter

ds = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
lat= ds.variables["XLAT"][0,:,0]
lon = ds.variables["XLONG"][0,0,:]

ds1 =open_mfdataset([r"I:\Merra 2 M2T1NXAER\regridded\apr\10.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\11.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\12.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\13.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\14.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\15.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\16.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\17.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\18.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\apr\19.nc4",])


#bcsmass = ds1.variables["BCSMASS"][:,:,:]
#bcsmass = np.nanmean(bcsmass,axis=0)*10**9 ##CONVERT FROM KG/M3 TO UG/M3

#dusmass = ds1.variables["DUSMASS"][:,:,:]
#dusmass = np.nanmean(dusmass,axis=0)*10**9
#ocsmass = ds1.variables["OCSMASS"][:,:,:]
#ocsmass = np.nanmean(ocsmass,axis=0)*10**9
so4smass = ds1.variables["SO4SMASS"][:,:,:]
so4smass = np.nanmean(so4smass,axis=0)*10**9


#ds2 = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
ds2 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
RAINC231 = getvar(ds2, "RAINC", timeidx=0)
cart_proj = get_cartopy(RAINC231)
#ax = plt.axes(projection=ccrs.PlateCarree())
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


##for gridlines
'''
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.1, linestyle='-')
gl.xlabels_top = False
gl.xlabels_bottom = True
gl.ylabels_left = True
gl.ylabels_right = False
gl.xlines = True
gl.ylines = True
gl.xlocator = mticker.FixedLocator([88, 90, 92, 94, 96, 98,100])
gl.ylocator = mticker.FixedLocator([20, 22, 24, 26, 28, 30])
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 15, 'color': 'gray'}
gl.xlabel_style = {'color': 'black', 'weight': 'normal'}
'''
#ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())

ax.set_xticks([75, 80, 85, 90, 95, 100], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())


#v = np.linspace(np.min(precp),np.max(precp),21)
#v = [-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0]
#v = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]
#v = [-1.96,0,1.96]
#v= np.linspace(-30.0,30,11)
#v= np.linspace(0,50,21)
#v= np.linspace(-50.0,50.0,21)
#v= np.linspace(-0.006,0.006,21)
v= np.linspace(0,10,21)
plt.contourf(lon,lat,so4smass,v,cmap="jet",transform=cartopy.crs.PlateCarree(),extend='max',antialiased='False')
#plt.pcolormesh(lon,lat,precp,vmin=-1.96,vmax=1.96,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
plt.colorbar(ax=ax,ticks=v,shrink=0.96)
#plt.title("(d) decrease 0-0.21 mm/hr")

plt.savefig(r'F:\WRF-CHEM ANALYSIS Chap 5\concentration and fraction\Merra2_so4.png',dpi=600,bbox_inches='tight')
plt.show()




