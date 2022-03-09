import pandas as pd
from netCDF4 import Dataset
from wrf import get_cartopy,getvar
import numpy as np
from xarray import open_mfdataset
#from nco import *
import xarray as xr
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter
import warnings
warnings.filterwarnings("ignore")

#ds = open_mfdataset(r"I:\Merra 2 M2T1NXAER\columnintegrated\MERRA2_columnintegrated_aerosol_mass_dec.nc")
ds = open_mfdataset(r"F:\WRF-CHEM ANALYSIS Chap 5\rainfall\apr_rainfall_evaluation_NOR.nc")
#ds = open_mfdataset(r"F:\WRF-CHEM ANALYSIS Chap 2\WRF-Chem rainfall evaluation\apr_rainfall_evaluation_nor_new.nc")
#ds = open_mfdataset(r"F:\WRF-CHEM ANALYSIS NEW\rainfall\New folder\apr_rainfall(RAINNC)_by_frequency-NOR-NO_BC_ABS.nc")
wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf = open_mfdataset('F:\IMD\z_modified_statistics_DEC.nc4')

lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]
#precp = ds.variables["AOD"][9:19,:,:]
#precp = ds.variables["BC_conc"][0,:,:]
#precp = ds.variables["NH4_FRACTION"][0,31:83,70:134]
precp = ds.variables["mean error"][:,:]
#print(np.size(precp))

precp = np.asarray(precp)
print(precp)
'''
for i in range(0,231,1):
    for j in range(0,182,1):
        precp[i,j]=np.nan

for i in range(0,123,1):
    for j in range(182,273,1):
        precp[i,j]=np.nan

for i in range(215,231,1):
    for j in range(182,273,1):
        precp[i,j]=np.nan

for i in range(0,231,1):
    for j in range(273,299,1):
        precp[i,j]=np.nan
'''
#print(precp)
#precp=np.reshape(precp,(1,69069))
#print(precp)
#precp= pd.DataFrame(precp)
#precp.to_csv('precp.csv')
#precp= precp.dropna(axis=0)
#precp=np.asarray(precp)
#print(np.size(precp))
#precp=np.isnan(precp)
#print(precp)
#print(np.size(precp))
#precp1 = precp[157:231, 182:299]
#precp1 = precp[157:204, 174:263]
#print(np.nanmean(precp1))
#'''
precp1 = precp[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = precp[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = precp[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = precp[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = precp[168:178, 174:216]
print(np.nanmean(precp1))
'''
#precp1 = precp[123:203, 174:252]
#print(np.nanmean(precp1))
precp1 = precp[:, :]
print(np.nanmean(precp1))
#precp1 = precp[157:186, 174:216]
#print(np.nansum(precp1))
#precp1 = precp[157:186, 216:252]
#print(np.nansum(precp1))
#precp1 = precp[186:203, 216:252]
#print(np.nansum(precp1))
#precp1 = precp[123:157, 200:237]
#print(np.nansum(precp1))
#precp1 = precp[168:178, 174:216]
#print(np.nansum(precp1))
#precp1 = precp[123:203, 174:252]
#print(np.nansum(precp1))

#z1 = precp1[np.where(precp1 > 0)]
#print(z1)
#z1 = np.size(z1)
#print(z1)

#z2 = precp1[np.where(precp1 == 0)]
#z2 = np.size(z2)
#print(z2)

#z3 = precp1[np.where(precp1 < 0)]
#z3 = np.size(z3)
#print(z3)

#print('total gridpoints=',z1+z2+z3)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }

ds2 = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
RAINC231 = getvar(ds2, "RAINC", timeidx=0)

#RAINNC766_1 = getvar(wrf, "RAINNC", timeidx=0)
#RAINNC766_1 = wrf.variables["value"][:,:]
#plt.figure(figsize=(10,6))
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

ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())


#ax.set_xticks([75, 80, 85, 90, 95, 100], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())
#v = np.linspace(np.min(precp),np.max(precp),21)
#v = [-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0]
#v = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]
#v = [-1.96,0,1.96]
#v= np.linspace(-30.0,30,11)
v= np.linspace(0,10,21)
#v= np.linspace(-50.0,50.0,21)
#print(lon)
#print(lat)

#plt.contourf(lon,lat,precp,v,cmap="jet",orientation="vertical",transform=cartopy.crs.PlateCarree(),extend='max',antialiased='False')
#plt.pcolormesh(lon,lat,precp,vmin=-1.96,vmax=1.96,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
#plt.colorbar(ax=ax,ticks=v,shrink=0.96)
#plt.title("(d) decrease 0-0.21 mm/hr")
'''
#plt.savefig(r'F:\WRF-CHEM ANALYSIS NEW\rainfall\New folder\RAIN_dec_NOR-NOBCABS(OLD).png',dpi=600,bbox_inches='tight')
#plt.savefig(r'D:\PHD\My PhD Reports\Manuscript3\Revision 1\mean_may.png',dpi=600,bbox_inches='tight')
#plt.savefig(r'F:\WRF-CHEM ANALYSIS\aerosol fraction and concentration\BC_frcation_april_level 15.png',dpi=600,bbox_inches='tight')
#plt.show()


