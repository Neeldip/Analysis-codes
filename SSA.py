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
#wrf_file = Dataset(r"G:\WRF_Chem_Output\201\April\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf = Dataset(r"G:\WRF_Chem_Output\202\GF\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")


time = ALL_TIMES
#TAU300 = getvar(wrf_file, "TAUAER1", timeidx=time)
TAU400 = getvar(wrf, "WAER2", timeidx=time)
TAU600 = getvar(wrf, "WAER3", timeidx=time)
#TAU999 = getvar(wrf_file, "TAUAER4", timeidx=time)

Z = TAU400
X = TAU400
Y = TAU600
if time == ALL_TIMES:
        X = X.mean("Time")
        Y = Y.mean("Time")

lats, lons = w.latlon_coords(Z)

##AOD_550nm calculation
a = -(np.log(X/Y))/(np.log(400/600))
TAU550 = X*((550/400)**(-a))
x = TAU550.mean(dim='bottom_top')

plt.figure(figsize=(10,6))
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

#ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())


ax.set_xticks([75, 80, 85, 90, 95, 100], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())


v = np.linspace(np.min(x),np.max(x),21)
#v = np.linspace(0,1,21)
plt.contourf(lons, lats, x, v,cmap=plt.get_cmap("jet"),extend='neither',transform=cartopy.crs.PlateCarree()) ###use if conc does not cahnege much between plots
#plt.contourf(lons, lats, CONC, cmap=plt.get_cmap("nipy_spectral"), extend='max')

#plt.xlabel("Longitude",fontdict=font)
#plt.ylabel("Latitude",fontdict=font)
#plt.title("(b) WRF-Chem 2xBC" ,weight='bold',fontdict=font)
#plt.title("(b) WRF-Chem 2xBC")

##SET TICK LABELS PROPERTIES####
ax=plt.gca()
#ax.set_xticklabels(ax.get_xticks(),font)
#ax.set_yticklabels(ax.get_yticks(),font)

#plt.colorbar(shrink=.90)
cb=plt.colorbar(fraction=0.046, pad=0.04, ticks=v).set_label(label='',size=15,weight='bold')
#cb = plt.colorbar(fraction=0.046, pad=0.04).set_label(label='', size=15, weight='bold')

plt.savefig("F:\WRF-CHEM ANALYSIS NEW\SSA\SSA_april_10-19.png",dpi=600,bbox_inches='tight')
plt.show()

