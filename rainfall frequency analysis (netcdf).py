from xarray import open_dataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import netCDF4 as nc
#ds = open_dataset(r"F:\WRF-CHEM ANALYSIS Chap 3\WRF-Chem rainfall evaluation\2018-Apr_regridded.nc")
#wrf =  open_dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf =  open_dataset(r"G:\WRF_Chem_Output\202\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf1 =  open_dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf1 =  open_dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")

wrf =  open_dataset(r"G:\WRF_Chem_Output\202\NO_BC_ABS\wrfout_d01_2018-04-10_00%3A00%3A00",engine='netcdf4')
wrf1 =  open_dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00",engine='netcdf4')

wrf_lat= wrf.variables["XLAT"][0,:,0]
wrf_lon = wrf.variables["XLONG"][0,0,:]

fn = r'F:\WRF-CHEM ANALYSIS Chap 5\rainfall\apr_rainfall_RAINC_frequency_analysis-NOR-NOBCABS.nc'   #<----------------
ds1 = nc.Dataset(fn, 'w', format='NETCDF4')

lat = ds1.createDimension('lat', 231)
lon = ds1.createDimension('lon', 299)
lats = ds1.createVariable('lat', 'f4', ('lat',))
lons = ds1.createVariable('lon', 'f4', ('lon',))

low_inc = ds1.createVariable('inc_0_5', 'f4', ('lat', 'lon',),fill_value=np.nan)
moderate_inc= ds1.createVariable('inc_5_10', 'f4', ('lat', 'lon',),fill_value=np.nan)
heavy_inc=ds1.createVariable('inc_>10', 'f4', ('lat', 'lon',),fill_value=np.nan)
low_dec = ds1.createVariable('dec_0_5', 'f4', ('lat', 'lon',),fill_value=np.nan)
moderate_dec= ds1.createVariable('dec_5_10', 'f4', ('lat', 'lon',),fill_value=np.nan)
heavy_dec=ds1.createVariable('dec_>10', 'f4', ('lat', 'lon',),fill_value=np.nan)
no_change=ds1.createVariable('no_change', 'f4', ('lat', 'lon',),fill_value=np.nan)

lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = wrf_lat
lons[:] = wrf_lon

for i in range(0,231,1):
    for j in range(0,299,1):
    #    Rain = ds.variables["RAINFALL"][10:20,i,j]
    #    if np.isnan(Rain).any() == True:
    #        low_inc[i, j] = np.nan
    #        moderate_inc[i, j] = np.nan
    #        heavy_inc[i, j] = np.nan
    #        low_dec[i, j] = np.nan
    #        moderate_dec[i, j] = np.nan
    #        heavy_dec[i, j] = np.nan
    #        no_change[i, j] = np.nan
    #    else:
            RAIN= wrf.variables['RAINC'][0:241,i,j]#+wrf.variables['RAINNC'][0:241,i,j]
            RAIN1 = wrf1.variables['RAINC'][0:241, i, j]# + wrf1.variables['RAINNC'][0:241, i, j]
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
            z2 = rain_diff[np.where(rain_diff > 0.21)]
            z2 = np.size(z2)
            z6 = rain_diff[np.where(rain_diff > 0.42)]
            z6 = np.size(z6)
            z3 = rain_diff[np.where(rain_diff < 0)]
            z3 = np.size(z3)
            z4 = rain_diff[np.where(rain_diff < -0.21)]
            z4 = np.size(z4)
            z7 = rain_diff[np.where(rain_diff < -0.42)]
            z7 = np.size(z7)
            z5 = rain_diff[np.where(rain_diff == 0)]
            z5 = np.size(z5)

            low_inc[i,j] = ((z1 - z2) / 240) * 100
            moderate_inc[i,j] = ((z2 - z6) / 240) * 100
            heavy_inc[i,j] = (z6 / 240) * 100
            low_dec[i,j] = ((z3 - z4) / 240) * 100
            moderate_dec[i,j] = ((z4 - z7) / 240) * 100
            heavy_dec[i,j] = (z7 / 240) * 100
            no_change[i,j] = (z5 / 240) * 100
            print(i)
            print(j)

print("low inc")
precp1 = low_inc[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = low_inc[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = low_inc[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = low_inc[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = low_inc[168:178, 174:216]
print(np.nanmean(precp1))

print("moderate inc")
precp1 = moderate_inc[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = moderate_inc[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = moderate_inc[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = moderate_inc[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = moderate_inc[168:178, 174:216]
print(np.nanmean(precp1))

print("high inc")
precp1 = heavy_inc[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = heavy_inc[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = heavy_inc[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = heavy_inc[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = heavy_inc[168:178, 174:216]
print(np.nanmean(precp1))

print("low dec")
precp1 = low_dec[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = low_dec[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = low_dec[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = low_dec[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = low_dec[168:178, 174:216]
print(np.nanmean(precp1))

print("moderate dec")
precp1 = moderate_dec[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = moderate_dec[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = moderate_dec[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = moderate_dec[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = moderate_dec[168:178, 174:216]
print(np.nanmean(precp1))

print("high dec")
precp1 = heavy_dec[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = heavy_dec[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = heavy_dec[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = heavy_dec[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = heavy_dec[168:178, 174:216]
print(np.nanmean(precp1))
ds1.close()
