import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter
#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset(r"G:\WRF_Chem_Output\202\NO_BC_ABS\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")

time = ALL_TIMES
#time = 60
#level = range(0,25,1)
C = getvar(wrf_file, "RAINC", timeidx=1)
#CLD = getvar(wrf_file, "SWDNBC", timeidx=time)
#CLD1 = getvar(wrf_file, "SWDNB", timeidx=time)
#CLD2 = getvar(wrf_file1, "SWDNBC", timeidx=time)
#CLD3 = getvar(wrf_file1, "SWDNB", timeidx=time)

CLD = getvar(wrf_file, "SWDNBC", timeidx=time)
CLD1 = getvar(wrf_file, "SWDNB", timeidx=time)
CLD2 = getvar(wrf_file1, "SWDNBC", timeidx=time)
CLD3 = getvar(wrf_file1, "SWDNB", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")
#lats, lons = w.latlon_coords(CLD)

#CLD_DIFF = (CLD-CLD1)-(CLD2-CLD3)
CLD_DIFF = (CLD2 - CLD)
lats,lons =w.latlon_coords(CLD_DIFF)
#print(lats)
#print(lons)


precp1 = CLD_DIFF[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = CLD_DIFF[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = CLD_DIFF[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = CLD_DIFF[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = CLD_DIFF[168:178, 174:216]
print(np.nanmean(precp1))

font = {'family': 'serif',
     'color':  'black',
     'weight': 'bold',
     'size': 12,
       }
cart_proj = get_cartopy(C)
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

#ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())

ax.set_xticks([75, 80, 85, 90, 95,100], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())

#ax.set_xlim(cartopy_xlim(RAINNC766_1))
#ax.set_ylim(cartopy_ylim(RAINNC766_1))
#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),20) ##COLORBAR LIMITS
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF), 30)
v = np.linspace(-30.0,30.0,21) ##RAINC LIMITS
plt.contourf(lons, lats, CLD_DIFF,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
plt.colorbar(ax=ax, shrink=.81,ticks=v)
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
plt.savefig(r"F:\WRF-CHEM ANALYSIS NEW\radiation\RADIATION_AEROSOL_NOR-NOBCABS.png",dpi=600,bbox_inches='tight')
plt.show()
CLD=None
CLD1=None
CLD2=None
CLD3=None