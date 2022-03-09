import pandas as pd
import wrf as w
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim,ll_to_xy
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from matplotlib import cbook
from netCDF4 import Dataset
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
#from mpl_toolkits.basemap import *
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter
from xarray import open_dataset


#ds = Dataset("G:\WRF_Chem_Output\April\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset("G:\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_dust_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"G:\WRF_Chem_Output\202\NO_BC_ABS\wrfout_d01_2018-04-10_00%3A00%3A00")
ds = Dataset(r"G:\WRF_Chem_Output\202\NO_BC_ABS\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"G:\WRF_Chem_Output\202\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")

RAINC23 = w.getvar(ds, "RAINC", timeidx=0)
RAINNC23 = w.getvar(ds, "RAINNC", timeidx=0)
RAINC766 = w.getvar(ds, "RAINC", timeidx=240)
RAINNC766 = w.getvar(ds, "RAINNC", timeidx=240)
RAINC= RAINC766 - RAINC23
RAINNC =RAINNC766 - RAINNC23
RAIN=RAINC + RAINNC

#RAIN= (RAINC/RAIN)*100

#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_dust_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_dust_bc_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
ds1 = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_bc_dust_emission\wrfout_d01_2018-04-10_00%3A00%3A00")

RAINC23_1 = w.getvar(ds1, "RAINC", timeidx=0)
RAINNC23_1 = w.getvar(ds1, "RAINNC", timeidx=0)
RAINC766_1 = w.getvar(ds1, "RAINC", timeidx=240)
RAINNC766_1 = w.getvar(ds1, "RAINNC", timeidx=240)
RAINC1=RAINC766_1 - RAINC23_1
RAINNC1 =RAINNC766_1 - RAINNC23_1
RAIN1=RAINC1  + RAINNC1

#DIFF= (RAINNC1/RAIN1)*100
#index=ll_to_xy(ds1,latitude=30.3996734619,longitude=93.194915771484375)
#print(index)

DIFF = ((RAIN1-RAIN))#/RAIN)*100
#precp1 = DIFF[157:204, 182:263]
precp1 = DIFF[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = DIFF[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = DIFF[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = DIFF[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = DIFF[168:178, 174:216]
print(np.nanmean(precp1))

'''
##########################################
#ds = Dataset("G:\WRF_Chem_Output\April\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset("G:\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_dust_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"G:\WRF_Chem_Output\202\\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"G:\WRF_Chem_Output\202\\NO_BC_ABS\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"G:\WRF_Chem_Output\202\\NO_BC_ABS\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"G:\WRF_Chem_Output\202\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")

RAINC23 = w.getvar(ds, "RAINC", timeidx=0)
RAINNC23 = w.getvar(ds, "RAINNC", timeidx=0)
RAINC766 = w.getvar(ds, "RAINC", timeidx=240)
RAINNC766 = w.getvar(ds, "RAINNC", timeidx=240)
RAINC= RAINC766 - RAINC23
RAINNC =RAINNC766 - RAINNC23
RAIN=RAINC + RAINNC

#RAIN= (RAINC/RAIN)*100

#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_dust_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_dust_bc_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"G:\WRF_Chem_Output\202\\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"G:\WRF_Chem_Output\202\\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")

RAINC23_1 = w.getvar(ds1, "RAINC", timeidx=0)
RAINNC23_1 = w.getvar(ds1, "RAINNC", timeidx=0)
RAINC766_1 = w.getvar(ds1, "RAINC", timeidx=240)
RAINNC766_1 = w.getvar(ds1, "RAINNC", timeidx=240)
RAINC1=RAINC766_1 - RAINC23_1
RAINNC1 =RAINNC766_1 - RAINNC23_1
RAIN1=RAINC1  + RAINNC1

#DIFF= (RAINNC1/RAIN1)*100
#index=ll_to_xy(ds1,latitude=30.3996734619,longitude=93.194915771484375)
#print(index)

DIFF1 = ((RAIN1-RAIN)/RAIN)*100
###############################################

DIFF2 = DIFF1-DIFF
'''
'''
for i in range(0,231,1):
    for j in range(0,182,1):
        DIFF[i,j]=np.nan

for i in range(0,123,1):
    for j in range(182,273,1):
        DIFF[i,j]=np.nan

for i in range(215,231,1):
    for j in range(182,273,1):
        DIFF[i,j]=np.nan

for i in range(0,231,1):
    for j in range(273,299,1):
        DIFF[i,j]=np.nan
ds1 = open_dataset(r"F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\2018-Apr_regridded.nc")
for i in range(0,231,1):
    for j in range(0,299,1):
        Rain = ds1.variables["RAINFALL"][0,i,j]
        if np.isnan(Rain).any() == True:
            DIFF[i,j]=np.nan
            #print(i)
            #print(j)

'''
#print(np.nanmean(DIFF))
#ds2 = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
#RAINC231 = w.getvar(ds2, "RAINC", timeidx=0)
lats,lons =w.latlon_coords(RAINNC766_1)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
#cart_proj = get_cartopy(RAINC231)
ax = plt.axes(projection=ccrs.PlateCarree())
#ax = plt.axes(projection=cart_proj)
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
#ax.set_xlim(cartopy_xlim(RAINC231))
#ax.set_ylim(cartopy_ylim(RAINC231))

'''
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='black', alpha=0.25, linestyle='-')
gl.xlabels_top = False
gl.xlabels_bottom = True
gl.ylabels_left = True
gl.ylabels_right = False
gl.xlines = True
gl.ylines = True
gl.xlocator = mticker.FixedLocator([88, 90, 92, 94, 96, 98])
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

#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(DIFF),np.max(DIFF),21) ##COLORBAR LIMITS

#v = np.linspace(-250,250,16) ##RAINNC LIMITS
#v = np.linspace(-100,100,16) ##RAINC LIMITS
#v = np.linspace(-100,100,21) ##TOTAL LIMITS
v = np.linspace(-100,100,21) ##TOTAL LIMITS

plt.contourf(lons, lats, DIFF,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.contourf(lons, lats, DIFF,v, cmap=plt.get_cmap("jet"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
plt.colorbar(ax=ax, shrink=.93,ticks=v)
#plt.title("RAINCcontribution_change_aerfeedback - no_bc_absorbtion rainfall difference")
#plt.figure(figsize=(8,6))
plt.savefig(r'F:\WRF-CHEM ANALYSIS Chap 5\rainfall\effectwise rainfall\RAIN_201-202 NOR-NOBCABS.png',dpi=600,bbox_inches='tight')
plt.show()