import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
from xarray import open_dataset
#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
from netCDF4 import Dataset
import cartopy.crs as ccrs

wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")

#wrf_file1=Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_lat = wrf_file.variables["XLAT"][0,:,0]
wrf_lon = wrf_file.variables["XLAT"][0,0,:]

#time = ALL_TIMES
time = 0
BIN1 = getvar(wrf_file, "QCLOUD", timeidx=time)
#BIN2 = getvar(wrf_file, "bc_a02", timeidx=time)
#BIN3 = getvar(wrf_file, "bc_a03", timeidx=time)
#BIN4 = getvar(wrf_file, "bc_a04", timeidx=time)
#PH = getvar(wrf_file,"PH",timeidx=time)
#PHB = getvar(wrf_file,"PHB",timeidx=time)
#T = getvar(trad, "TEMPERATURE", timeidx=time) ##sensible temperature,units = K
#P = getvar(trad, "PRESSURE", timeidx=time)  ##units = Pascal
#R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 #+ BIN2 + BIN3 + BIN4
#BIN1=None
#BIN2=None
#BIN3=None
#BIN4=None
#GPH=(PH+PHB)/9.81
#PH=None
#PHB=None

if time == ALL_TIMES:
#    T = T.mean("Time")
#    P = P.mean("Time")
    X = X.mean("Time")
#    GPH= GPH.mean('Time')

#X=X.mean(axis=0)
#T=T.mean(axis=0)
#P=P.mean(axis=0)
#GPH=GPH.mean(axis=0)

CONC= X#*P)/(R*T)

#ds= Dataset(r"F:\WRF-CHEM ANALYSIS\CLDFRA_column_integrated_sum_NOR.nc",mode='w',format='NETCDF4')

#lat = ds.createDimension('lat', 231)
#lon = ds.createDimension('lon', 299)
#lats = ds.createVariable('lat', 'f4', ('lat',))
#lons = ds.createVariable('lon', 'f4', ('lon',))
#value = ds.createVariable('column mass sum', 'f4', ('lat', 'lon',),fill_value=np.nan)
#value.units = 'column_integrated_sum kg/kg'
#lats.units = 'degrees north'
#lons.units = 'degrees east'

#lats[:] = wrf_lat
#lons[:] = wrf_lon

value=np.empty((231,299))
for i in range(0,231,1):
    for j in range(0,299,1):
        column=[]
        for k in range(0,44,1):
            print(i)
            print(j)
            #print(k)
            conc = (CONC[k,i,j])#*(GPH[k+1,i,j]-GPH[k,i,j])   ###ug/m3 *m = ug/m2
            column.append(conc)
            if len(column)==44:
                column_sum= np.sum(column)
                value[i,j]=column_sum
                #print(value)
#ds.close()

CONC=None
X=None


font = {'family': 'serif',
            'color': 'black',
            'weight': 'bold',
            'size': 12,
            }
cart_proj = get_cartopy(BIN1)
ax = plt.axes(projection=ccrs.PlateCarree())
#ax = plt.axes(projection=cart_proj)
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
plt.contourf(wrf_lon, wrf_lat, value, v,cmap=plt.get_cmap("nipy_spectral"),extend='max',transform=cartopy.crs.PlateCarree()) ###use if conc does not cahnege much between plots
#plt.contourf(lons, lats, CONC, cmap=plt.get_cmap("nipy_spectral"), extend='max')
#font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
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

