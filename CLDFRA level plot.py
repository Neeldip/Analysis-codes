import wrf as w
from wrf import latlon_coords,omp_set_num_threads,get_cartopy,getvar,ALL_TIMES,cartopy_ylim,cartopy_xlim
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter

#from xarray import Dataset
# https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
# trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
# wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"K:\WRF_Chem_Output\202\April\only_emiss_ne\wrfout_d01_2018-04-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
# wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
# wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
# wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\NO_BC_ABS\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
time = ALL_TIMES

#time = 60
level = range(0, 44, 1)  # [5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,30]
#level = [11,13,15,17]
#level = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
#level = [15,16,17,18,19,20,21,22,23,24,25]
C = getvar(wrf_file, "V10", timeidx=1)
CLD = getvar(wrf_file, "CLDFRA", timeidx=time)
#CLD1 = getvar(wrf_file1, "TEMPERATURE", timeidx=time)
#print(CLD1)

if time == ALL_TIMES:
    CLD = CLD.mean(axis=0)
    #CLD1 = CLD1.mean(axis=0)
    # CLD = CLD.std("Time")
    # CLD1 = CLD1.std("Time")

for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    #z = CLD1[y, :, :]
    CLD_DIFF = x# ((x - z))/10**7   #* 100  # 0000#*1000000000 ###*100 FOR BETTER READABILITY

    lats, lons = latlon_coords(C)
    font = {'family': 'serif',
            'color': 'black',
            'weight': 'bold',
            'size': 12,
            }

    #ds2 = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
    ds2=Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
    RAINC231 = getvar(ds2, "RAINC", timeidx=0)

    # RAINNC766_1 = getvar(wrf, "RAINNC", timeidx=0)
    # RAINNC766_1 = wrf.variables["value"][:,:]
    # plt.figure(figsize=(10,6))
    cart_proj = get_cartopy(RAINC231)
    #ax = plt.axes(projection=ccrs.PlateCarree())
    ax = plt.axes(projection=cart_proj)
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    #ax.add_feature(states, linewidth=.5, edgecolor="black")
    ax.coastlines('50m', linewidth=0.8)
    # ax.add_feature(cartopy.feature.LAND)
    ax.add_feature(cartopy.feature.OCEAN)
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.add_feature(cartopy.feature.LAKES)
    # ax.add_feature(cartopy.feature.RIVERS)

    ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
    ax.set_xlim(cartopy_xlim(RAINC231))
    ax.set_ylim(cartopy_ylim(RAINC231))

    #ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
    #ax.xaxis.set_major_formatter(LongitudeFormatter())
    #ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
    #ax.yaxis.set_major_formatter(LatitudeFormatter())

    ax.set_xticks([75, 80, 85, 90, 95,100], crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
    ax.yaxis.set_major_formatter(LatitudeFormatter())


    # ax.gridlines(color="black", linestyle="dotted")
    #v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),21) ##COLORBAR LIMITS
    v = np.linspace(0.0, 1.0, 21)
    plt.contourf(lons, lats, CLD_DIFF, v, cmap=plt.get_cmap("jet"), antialiased=False,
                 transform=cartopy.crs.PlateCarree(), extend='max')
    plt.colorbar(ax=ax, shrink=.86, ticks=v)#,label='x10^7')
    # plt.title("CLDFRA change at level="+str(y),weight='bold',fontdict=font)
    # plt.figure(figsize=(8,6))

    plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 6\WRF level analysis\cldfra\CLDFRA_NOR_LEVEL =" + str(y) + ".png",dpi=600, bbox_inches='tight')
    plt.show()

CLD = None
CLD1 = None