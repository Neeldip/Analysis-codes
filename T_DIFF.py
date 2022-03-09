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

wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset("G:\WRF_Chem_Output\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

time = ALL_TIMES
#time = 60
level = [5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,30]
#level = range(0,25,1)
C = getvar(wrf_file, "TEMPERATURE", timeidx=1)
CLD = getvar(wrf_file, "TEMPERATURE", timeidx=time)
CLD1 = getvar(wrf_file1, "TEMPERATURE", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")

for y in level:
    #lats, lons = w.latlon_coords(CLD)
    x = CLD[y,:,:]
    z = CLD1[y, :, :]
    CLD_DIFF = x - z
    #v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),20)
    #plt.contourf(lons, lats, CLD_DIFF,v, cmap=plt.get_cmap("nipy_spectral"),extend='max') ###use if conc does not cahnege much between plots
    #plt.contourf(lons, lats, CONC, cmap=plt.get_cmap("nipy_spectral"), extend='max')
    #font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
    #plt.xlabel("Longitude",fontdict=font)
    #plt.ylabel("Latitude",fontdict=font)
    #plt.title("CLOUD FRACTION(%) increase at level="+str(y),weight='bold',fontdict=font)

##SET TICK LABELS PROPERTIES####
    #ax=plt.gca()
    #ax.set_xticklabels(ax.get_xticks(),font)
    #ax.set_yticklabels(ax.get_yticks(),font)
    #bx = plt.axes(projection=get_cartopy(CLD))
#plt.colorbar(shrink=.90)
    #cb=plt.colorbar(fraction=0.046, pad=0.04, ticks=v).set_label(label='',size=15,weight='bold')
    #cb = plt.colorbar(fraction=0.046, pad=0.04).set_label(label='', size=15, weight='bold')
    #bx.gridlines(color="black", linestyle="dotted")
    #plt.savefig("CLDFRA AT LEVEL ="+str(y)+".png",dpi=300)
    #plt.show()

    lats,lons =w.latlon_coords(CLD_DIFF)
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
    ax.add_feature(cartopy.feature.BORDERS, linestyle='--')
    #ax.set_xlim(cartopy_xlim(RAINNC766_1))
    #ax.set_ylim(cartopy_ylim(RAINNC766_1))
    #ax.gridlines(color="black", linestyle="dotted")
    #v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),30) ##COLORBAR LIMITS
    v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF), 30)
    plt.contourf(lons, lats, CLD_DIFF,v, cmap=plt.get_cmap("nipy_spectral"),antialiased=False,
                 transform=cartopy.crs.PlateCarree())
    plt.colorbar(ax=ax, shrink=.81)
    plt.title("Temperature difference at level="+str(y),weight='bold',fontdict=font)
    #plt.figure(figsize=(8,6))
    plt.savefig("Temperature difference AT LEVEL ="+str(y)+".png",dpi=300)
    plt.show()