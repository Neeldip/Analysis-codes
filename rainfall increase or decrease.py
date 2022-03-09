from xarray import open_dataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import netCDF4 as nc
from netCDF4 import Dataset
from wrf import getvar,get_cartopy,cartopy_xlim,cartopy_ylim
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter


#ds = open_dataset(r"F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\2018-Apr_regridded.nc")
#wrf =  open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf =  open_dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf =  open_dataset(r"G:\WRF_Chem_Output\202\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf =  open_dataset(r"G:\WRF_Chem_Output\202\NO_BC_ABS\wrfout_d01_2018-04-10_00%3A00%3A00",engine='netcdf4')
wrf1 =  open_dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00",engine='netcdf4')
#wrf1 =  open_dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
lat= wrf.variables["XLAT"][0,:,0]
lon = wrf.variables["XLONG"][0,0,:]

#fn = r'F:\WRF-CHEM ANALYSIS\apr(4NOR-NOBCABS)_rain_percentage_time_inc_or_dec .nc'   #<----------------
#ds1 = nc.Dataset(fn, 'w', format='NETCDF4')

#lat = ds1.createDimension('lat', 231)
#lon = ds1.createDimension('lon', 299)
#lats = ds1.createVariable('lat', 'f4', ('lat',))
#lons = ds1.createVariable('lon', 'f4', ('lon',))

#inc = ds1.createVariable('inc', 'f4', ('lat', 'lon',),fill_value=np.nan)
#dec = ds1.createVariable('dec', 'f4', ('lat', 'lon',),fill_value=np.nan)

#lats.units = 'degrees north'
#lons.units = 'degrees east'

#lats[:] = wrf_lat
#lons[:] = wrf_lon
inc=np.empty((231,299))
dec=np.empty((231,299))

for i in range(0,231,1):
    for j in range(0,299,1):
        #Rain = ds.variables["RAINFALL"][0,i,j]
        #if np.isnan(Rain).any() == True:
        #    inc[i, j] = np.nan
        #    dec[i, j] = np.nan
        #else:

            RAIN= wrf.variables['RAINC'][0:241,i,j]+wrf.variables['RAINNC'][0:241,i,j]
            RAIN1 = wrf1.variables['RAINC'][0:241, i, j]+ wrf1.variables['RAINNC'][0:241, i, j]
            rain= RAIN[1:241]-RAIN[0:240]
            rain1 = RAIN1[1:241] - RAIN1[0:240]
            rain_diff=rain1-rain
            #rain_acc_daily = []
            #for time in range(0, 240, 24):
            #    rain_daily = rain[time:time + 24]
            #    rain_daily_sum = np.sum(rain_daily)
            #   rain_acc_daily.append(rain_daily_sum)
            #    if len(rain_acc_daily) == 10:
            #        print(i)
            #        print(j)
            #        RMSE[i, j] = rmse(Rain, rain_acc_daily)
            #        IOA[i, j] = d(Rain, rain_acc_daily)
            #       MAE[i, j] = mae(Rain, rain_acc_daily)
            ##increase
            z1 = rain_diff[np.where(rain_diff > 0)]
            z1 = np.size(z1)
            z3 = rain_diff[np.where(rain_diff < 0)]
            z3=np.size(z3)
            inc[i, j] = (z1/ 240) * 100
            dec[i,j] = (z3 / 240) * 100
            print(i)
            print(j)

#ds1.close()
precp1 = inc[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = inc[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = inc[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = inc[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = inc[168:178, 174:216]
print(np.nanmean(precp1))

precp1 = dec[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = dec[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = dec[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = dec[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = dec[168:178, 174:216]
print(np.nanmean(precp1))

font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
#ds2 = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")

#RAINC231 = getvar(ds2, "RAINC", timeidx=0)
#RAINC231 = getvar(wrf, "RAINC", timeidx=0)
#RAINNC766_1 = getvar(wrf, "RAINNC", timeidx=0)
#RAINNC766_1 = wrf.variables["value"][:,:]
#plt.figure(figsize=(10,6))
#cart_proj = get_cartopy(var=RAINC231)
ax = plt.axes(projection=ccrs.PlateCarree())
#ax = plt.axes(projection=cart_proj)
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

#ax.set_xlim(cartopy_xlim(RAINC231))
#ax.set_ylim(cartopy_ylim(RAINC231))
#ax.set_xlim(left=87.65714,right=98.07542)
#ax.set_ylim(bottom=21.82169,top=30.20221)

#ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())

ax.set_xticks([75, 80, 85, 90, 95, 100], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())

#v = np.linspace(np.min(precp),np.max(precp),21)
#v = [-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0]
#v = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]
#v = [-1.96,0,1.96]
#v= np.linspace(-30.0,30,11)
#v= np.linspace(0,50,21)
v= np.linspace(0,50,21)
#print(lon)
#print(lat)

plt.contourf(lon,lat,inc,v,cmap="jet",orientation="vertical",transform=cartopy.crs.PlateCarree(),extend='max',antialiased='False')
#plt.pcolormesh(lon,lat,precp,vmin=-1.96,vmax=1.96,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
plt.colorbar(ax=ax,ticks=v,shrink=0.96)
#plt.title("(d) decrease 0-0.21 mm/hr")
plt.savefig(r'F:\WRF-CHEM ANALYSIS Chap 5\rainfall\RAIN_inc%.png',dpi=1200,bbox_inches='tight')
#plt.show()

plt.contourf(lon,lat,dec,v,cmap="jet",orientation="vertical",transform=cartopy.crs.PlateCarree(),extend='max',antialiased='False')
#plt.pcolormesh(lon,lat,precp,vmin=-1.96,vmax=1.96,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
#plt.colorbar(ax=ax,ticks=v,shrink=0.96)
#plt.title("(d) decrease 0-0.21 mm/hr")
plt.savefig(r'F:\WRF-CHEM ANALYSIS Chap 5\rainfall\RAIN_dec%.png',dpi=1200,bbox_inches='tight')
plt.show()
