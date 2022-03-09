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

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\February\wrfout_d01_2018-02-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\February\wrf_trad_fields_d01_2018-02-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\March\wrf_trad_fields_d01_2018-03-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\May\wrfout_d01_2018-05-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\May\wrf_trad_fields_d01_2018-05-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\June\wrfout_d01_2018-06-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\June\wrf_trad_fields_d01_2018-06-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\July\wrfout_d01_2018-07-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\July\wrf_trad_fields_d01_2018-07-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\August\wrfout_d01_2018-08-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\August\wrf_trad_fields_d01_2018-08-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\September\wrfout_d01_2018-09-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\September\wrf_trad_fields_d01_2018-09-10_00%3A00%3A00")

wrf_file = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

time = ALL_TIMES
#time = 60
#level = [0,5,10,15,20,25,30]
level = [0]
BIN1 = getvar(wrf_file, "oc_a01", timeidx=time)
BIN2 = getvar(wrf_file, "oc_a02", timeidx=time)
BIN3 = getvar(wrf_file, "oc_a03", timeidx=time)
BIN4 = getvar(wrf_file, "oc_a04", timeidx=time)
T = getvar(trad, "TEMPERATURE", timeidx=time) ##sensible temperature,units = K
P = getvar(trad, "PRESSURE", timeidx=time)  ##units = Pascal
R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 + BIN2 + BIN3 + BIN4
#BIN1=None
BIN2=None
BIN3=None
BIN4=None
if time == ALL_TIMES:
    T = T.mean("Time")
    P = P.mean("Time")
    X = X.mean("Time")

for y in level:
    lats, lons = w.latlon_coords(X)
    x = X[y,:,:]    ##unit = ug/kg dry air
    TEM = T[y,:,:]
    PRE = P[y,:,:]
    K = (R*TEM)/PRE     ### K unit = m^3*kg^-1 dry air

    CONC =x/K          ##unit = ug*m^-3


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

    #ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
    #ax.xaxis.set_major_formatter(LongitudeFormatter())
    #ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
    #ax.yaxis.set_major_formatter(LatitudeFormatter())

    ax.set_xticks([75, 80, 85, 90, 95, 100], crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
    ax.yaxis.set_major_formatter(LatitudeFormatter())

    v = np.linspace(np.min(CONC),10,21)
    plt.contourf(lons, lats, CONC, v,cmap=plt.get_cmap("jet"),extend='max',transform=cartopy.crs.PlateCarree()) ###use if conc does not cahnege much between plots
    #plt.contourf(lons, lats, CONC, cmap=plt.get_cmap("nipy_spectral"), extend='max')
    font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
    #plt.xlabel("Longitude",fontdict=font)
    #plt.ylabel("Latitude",fontdict=font)
    #plt.title("oin concentration mdm (ug m^-3) at level="+str(y),weight='bold',fontdict=font)

##SET TICK LABELS PROPERTIES####
    #ax=plt.gca()
    #ax.set_xticklabels(ax.get_xticks(),font)
    #ax.set_yticklabels(ax.get_yticks(),font)

#plt.colorbar(shrink=.90)
    cb=plt.colorbar(fraction=0.046, pad=0.04, ticks=v).set_label(label='',size=15,weight='bold')
    #cb = plt.colorbar(fraction=0.046, pad=0.04).set_label(label='', size=15, weight='bold')

    plt.savefig(r'F:\WRF-CHEM ANALYSIS Chap 5\concentration and fraction\oc_level_'+str(y)+ '.png',dpi=600,bbox_inches='tight')
    plt.show()

T=None
P=None
X=None
CONC=None