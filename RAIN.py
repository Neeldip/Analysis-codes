import wrf as w
from wrf import omp_set_num_threads, to_np,getvar,get_cartopy,cartopy_xlim,cartopy_ylim
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
from xarray import open_dataset
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter

#'''
#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrfout_d01_2018-04-10_00%3A00%3A00")

#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_No_NE_BC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_No_NE_2xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_Only_NE_BC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_Only_NE_2xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
ds =wrf =  Dataset(r"K:\WRF_Chem_Output\202\April\no_emiss_ne\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\wrf_lin\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds=wrf=Dataset("F:\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds =wrf =  Dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds =wrf =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds=wrf=Dataset(r'G:\WRF_Chem_Output\202\GF\wrfout_d01_2018-04-10_00%3A00%3A00')
#ds = wrf_file = [Dataset(r"G:\WRF_Chem_Output\202new\March\wrfout_d01_2018-03-11_00%3A00%3A00"),
#ds = wrf_file = [Dataset(r"G:\WRF_Chem_Output\202new\March\wrfout_d01_2018-03-11_00%3A00%3A00"),
RAINC23 = w.getvar(ds, "RAINC", timeidx=0)
RAINNC23 = w.getvar(ds, "RAINNC", timeidx=0)
RAINC766 = w.getvar(ds, "RAINC", timeidx=240)
RAINNC766 = w.getvar(ds, "RAINNC", timeidx=240)
lats = w.getvar(ds, "XLAT", timeidx=0)
#print(lats)
RAINC=RAINC766 - RAINC23
RAINNC =RAINNC766 - RAINNC23
#RAIN=(RAINNC/(RAINC + RAINNC))*100
RAIN=RAINC+RAINNC


precp1 = RAIN[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RAIN[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RAIN[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RAIN[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RAIN[168:178, 174:216]
print(np.nanmean(precp1))
'''
#################ARRAY MODIFICATION to select region################
for i in range(0,231,1):
    for j in range(0,169,1):
        RAIN[i,j]=np.nan

for i in range(0,123,1):
    for j in range(169,273,1):
        RAIN[i,j]=np.nan

for i in range(215,231,1):
    for j in range(169,273,1):
        RAIN[i,j]=np.nan

for i in range(0,231,1):
    for j in range(273,299,1):
        RAIN[i,j]=np.nan

ds1 = open_dataset(r"F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\2018-Apr_regridded.nc")
for i in range(0,231,1):
    for j in range(0,299,1):
        Rain = ds1.variables["RAINFALL"][0,i,j]
        if np.isnan(Rain).any() == True:
            RAIN[i,j]=np.nan
            print(i)
            print(j)
'''


######
lats,lons =w.latlon_coords(RAIN)

#cart_proj = w.get_cartopy(RAINC23)
#plt.figure(figsize=(6.692,6))
#states = NaturalEarthFeature(category="cultural", scale="50m",
#                             facecolor="none",
#                             name="admin_1_states_provinces_shp")
#plt.contour(lons, lats, RAIN, 10, colors="black")
#v = np.linspace(0,np.max(RAIN),21) ##COLORBAR LIMITS
v = np.linspace(0,237.0,21) ##COLORBAR LIMITS
#plt.contourf(lons, lats, RAIN, v,cmap=plt.get_cmap("jet"))
#plt.contourf(lons, lats, RAIN, v,cmap=plt.get_cmap("jet"),extend='max')
font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }

ds2 = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
#ds2 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")


RAINNC766_1 = getvar(ds2, "RAINNC", timeidx=0)
#RAINNC766_1 = wrf.variables["value"][:,:]
#plt.figure(figsize=(10,6))
cart_proj = get_cartopy(RAINNC766_1)
ax = plt.axes(projection=cart_proj)
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8,edgecolor="black")
#ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.COASTLINE,edgecolor="black")
#ax.add_feature(cartopy.feature.LAKES)
#ax.add_feature(cartopy.feature.RIVERS)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-',edgecolor="black")
ax.set_xlim(cartopy_xlim(RAINNC766_1))
ax.set_ylim(cartopy_ylim(RAINNC766_1))



#ax.set_xticks([75, 80, 85, 90, 95, 100], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())


ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())

##SET TICK LABELS PROPERTIES####
#ax=plt.gca()
#ax.set_xticklabels(ax.get_xticks(),font)
#ax.set_yticklabels(ax.get_yticks(),font)
#########
#plt.colorbar(shrink=.90)

#current_cmap=plt.cm.jet
#x=current_cmap.set_bad(color='white')

plt.contourf(lons, lats, RAIN, v,cmap=plt.get_cmap("jet"),transform=cartopy.crs.PlateCarree(),extend='max')
plt.colorbar(fraction=0.046, pad=0.04,ticks=v).set_label(label='',size=15,weight='bold')
#plt.xlabel("Longitude",fontdict=font)
#plt.ylabel("Latitude",fontdict=font)
#plt.title("Rainfall-bc_directON",weight='bold',fontdict=font)
plt.savefig(r'F:\WRF-CHEM ANALYSIS Chap 6\rainfall\RAIN_no_emiss_ne.png',dpi=600, bbox_inches='tight')
#plt.savefig(r'F:\WRF-CHEM ANALYSIS\rainfall\RAIN_NOR.png',dpi=600, bbox_inches='tight')
plt.show()








