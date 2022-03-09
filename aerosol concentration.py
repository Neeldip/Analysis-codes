import numpy as np
from numpy import transpose
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES,destagger)
import cartopy.crs as crs
#from cartopy.feature import NaturalEarthFeature
#from mpl_toolkits.basemap import Basemap
wrf_file = Dataset("G:\WRF_Chem_Output\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
time = ALL_TIMES
level = 0
BC1 = getvar(wrf_file, "water_a01", timeidx=time)
BC2 = getvar(wrf_file, "water_a02", timeidx=time)
BC3 = getvar(wrf_file, "water_a03", timeidx=time)
BC4 = getvar(wrf_file, "water_a04", timeidx=time)
BC = BC1+BC2+BC3+BC4
#print(BC)
if time == ALL_TIMES:
 BC = BC.mean("Time")
BC = BC/0.8163 ###BC unit is ug/kg dry air;convert tu ug/m3; 1kg dry air = 0.8163 m3
BC=BC[level,:,:]
#w = getvar(wrf_file,'wa', timeidx=time)
#w = w.mean("Time")
#w = w[level,:,:]
#print(w)
#BC = w*100 #ug/m2
lats,lons = latlon_coords(BC)
v = np.linspace(0, 100, 30, endpoint=True) ##COLORBAR & plot LIMITS
plt.contourf(lons, lats, BC,v,cmap=plt.get_cmap("nipy_spectral"))
font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
plt.xlabel("Longitude",fontdict=font)
plt.ylabel("Latitude",fontdict=font)
plt.title("BC(ug/m^3)",weight='bold',fontdict=font)
##SET TICK LABELS PROPERTIES####
ax=plt.gca()
ax.set_xticklabels(ax.get_xticks(),font)
ax.set_yticklabels(ax.get_yticks(),font)
#########
#plt.colorbar(shrink=.90)
cb=plt.colorbar(fraction=0.046, pad=0.04,ticks=v).set_label(label='',size=15,weight='bold')
#plt.clim(0,200) ###plot LIMITS
#plt.savefig('rain.png',dpi=500)
plt.show()