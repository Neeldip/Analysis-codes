import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")

wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
time = 202
#TAU300 = getvar(wrf_file, "TAUAER1", timeidx=time)
TAU400 = getvar(wrf_file, "TAUAER2", timeidx=time)
TAU600 = getvar(wrf_file, "TAUAER3", timeidx=time)
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
x = TAU550.sum(dim='bottom_top')
############NEW########

#wrf_file = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
time = 202
#TAU300 = getvar(wrf_file, "TAUAER1", timeidx=time)
TAU400 = getvar(wrf_file, "TAUAER2", timeidx=time)
TAU600 = getvar(wrf_file, "TAUAER3", timeidx=time)
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
m = TAU550.sum(dim='bottom_top')


####DIFFERENCE######
DIFF= ((x-m)/m)*100

font = {'family': 'serif',
     'color':  'black',
     'weight': 'bold',
     'size': 12,
       }
cart_proj = get_cartopy(TAU400)
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
v = np.linspace(-50,50,21)
plt.contourf(lons, lats, DIFF, v,cmap=plt.get_cmap("jet"),extend='both',transform=cartopy.crs.PlateCarree()) ###use if conc does not cahnege much between plots
#plt.contourf(lons, lats, CONC, cmap=plt.get_cmap("nipy_spectral"), extend='max')
font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
plt.xlabel("Longitude",fontdict=font)
plt.ylabel("Latitude",fontdict=font)
#plt.title("AOD_550nm_feedback-no_bc_absorbtion" ,weight='bold',fontdict=font)

##SET TICK LABELS PROPERTIES####
ax=plt.gca()
ax.set_xticklabels(ax.get_xticks(),font)
ax.set_yticklabels(ax.get_yticks(),font)

#plt.colorbar(shrink=.90)
cb=plt.colorbar(fraction=0.046, pad=0.04, ticks=v).set_label(label='',size=15,weight='bold')
#cb = plt.colorbar(fraction=0.046, pad=0.04).set_label(label='', size=15, weight='bold')

plt.savefig(r"F:\WRF-CHEM ANALYSIS\radiation\AOD_NORM-NOFEEDBACK%_time=202.png",dpi=300,bbox_inches='tight')
plt.show()