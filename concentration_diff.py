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

#trad = Dataset("G:\WRF_Chem_Output\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\wrfout_d01_2018-04-07_00%3A00%3A00")
#trad = Dataset(r"G:\WRF_Chem_Output\202\wrf_trad_fields_d01_2018-04-07_00%3A00%3A00")
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
trad1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

time = 203
#time = 60
#level = [0,5,10,15,20,25,30]
level = [9]
BIN1 = getvar(wrf_file, "bc_a01", timeidx=time)
BIN2 = getvar(wrf_file, "bc_a02", timeidx=time)
BIN3 = getvar(wrf_file, "bc_a03", timeidx=time)
BIN4 = getvar(wrf_file, "bc_a04", timeidx=time)
T = getvar(trad, "TEMPERATURE", timeidx=time) ##sensible temperature,units = K
P = getvar(trad, "PRESSURE", timeidx=time)  ##units = Pascal
R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 + BIN2 + BIN3 + BIN4
if time == ALL_TIMES:
    T = T.mean("Time")
    P = P.mean("Time")
    X = X.mean("Time")
BIN1=None
BIN2=None
BIN3=None
BIN4=None

BIN1 = getvar(wrf_file1, "bc_a01", timeidx=time)
BIN2 = getvar(wrf_file1, "bc_a02", timeidx=time)
BIN3 = getvar(wrf_file1, "bc_a03", timeidx=time)
BIN4 = getvar(wrf_file1, "bc_a04", timeidx=time)
T1 = getvar(trad1, "TEMPERATURE", timeidx=time) ##sensible temperature,units = K
P1 = getvar(trad1, "PRESSURE", timeidx=time)  ##units = Pascal
X1 = BIN1 + BIN2 + BIN3 + BIN4
if time == ALL_TIMES:
    T1 = T1.mean("Time")
    P1 = P1.mean("Time")
    X1 = X1.mean("Time")


for y in level:
    lats, lons = w.latlon_coords(X)
    x = X[y,:,:]    ##unit = ug/kg dry air.
    x1 = X1[y,:,:]
    TEM = T[y,:,:]
    TEM1= T1[y,:,:]
    PRE = P[y,:,:]
    PRE1= P1[y,:,:]
    K = (R*TEM)/PRE     ### K unit = m^3*kg^-1 dry air
    K1=(R*TEM1)/PRE1
    CONC =x/K          ##unit = ug*m^-3
    CONC1=x1/K1

    DIFF= CONC-CONC1

    font = {'family': 'serif',
            'color': 'black',
            'weight': 'bold',
            'size': 12,
            }
    cart_proj = get_cartopy(BIN1)
    ax = plt.axes(projection=cart_proj)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    # ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.add_feature(cartopy.feature.LAKES)
    # ax.add_feature(cartopy.feature.RIVERS)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-')


    v = np.linspace(-2,2,21)
    plt.contourf(lons, lats, CONC, v,cmap=plt.get_cmap("nipy_spectral"),extend='max',transform=cartopy.crs.PlateCarree()) ###use if conc does not cahnege much between plots
    #plt.contourf(lons, lats, CONC, cmap=plt.get_cmap("nipy_spectral"), extend='max')
    font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
    plt.xlabel("Longitude",fontdict=font)
    plt.ylabel("Latitude",fontdict=font)
    #plt.title("oin concentration mdm (ug m^-3) at level="+str(y),weight='bold',fontdict=font)

##SET TICK LABELS PROPERTIES####
    ax=plt.gca()
    ax.set_xticklabels(ax.get_xticks(),font)
    ax.set_yticklabels(ax.get_yticks(),font)

#plt.colorbar(shrink=.90)
    cb=plt.colorbar(fraction=0.046, pad=0.04, ticks=v).set_label(label='',size=15,weight='bold')
    #cb = plt.colorbar(fraction=0.046, pad=0.04).set_label(label='', size=15, weight='bold')

    #plt.savefig("bcconc_avg_level="+str(y)+".png",dpi=300)
    plt.show()
