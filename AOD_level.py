import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *
#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")

wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")

time = ALL_TIMES
#time = 60
level = [0,5,10,15,20,25,30]
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

for l in level:
    lats, lons = w.latlon_coords(Z)
    x = X[l,:,:]    ##unit = ug/kg dry air
    y = Y[l, :, :]
    a = -(np.log(x / y)) / (np.log(400 / 600))
    TAU550 = x * ((550 / 400) ** (-a))
    v = np.linspace(0,.2,30)
    plt.contourf(lons, lats, TAU550, v,cmap=plt.get_cmap("nipy_spectral"),extend='max') ###use if conc does not cahnege much between plots
    #plt.contourf(lons, lats, CONC, cmap=plt.get_cmap("nipy_spectral"), extend='max')
    font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
    plt.xlabel("Longitude",fontdict=font)
    plt.ylabel("Latitude",fontdict=font)
    plt.title("AOD_550nm at level="+str(l)+"time="+str(time),weight='bold',fontdict=font)

##SET TICK LABELS PROPERTIES####
    ax=plt.gca()
    ax.set_xticklabels(ax.get_xticks(),font)
    ax.set_yticklabels(ax.get_yticks(),font)

#plt.colorbar(shrink=.90)
    cb=plt.colorbar(fraction=0.046, pad=0.04, ticks=v).set_label(label='',size=15,weight='bold')
    #cb = plt.colorbar(fraction=0.046, pad=0.04).set_label(label='', size=15, weight='bold')

    plt.savefig("AOD_550nm at level="+str(l)+"time="+str(time)+".png",dpi=300)
    plt.show()
