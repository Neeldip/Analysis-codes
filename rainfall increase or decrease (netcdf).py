from xarray import open_dataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import netCDF4 as nc
ds = open_dataset(r"F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\2018-Apr_regridded.nc")
#wrf =  open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf =  open_dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf1 =  open_dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_lat= wrf.variables["XLAT"][0,:,0]
wrf_lon = wrf.variables["XLONG"][0,0,:]

fn = r'F:\WRF-CHEM ANALYSIS NEW\rainfall\apr(NOR-NOFEED)_rain_percentage_time_inc_or_dec .nc'   #<----------------
ds1 = nc.Dataset(fn, 'w', format='NETCDF4')

lat = ds1.createDimension('lat', 231)
lon = ds1.createDimension('lon', 299)
lats = ds1.createVariable('lat', 'f4', ('lat',))
lons = ds1.createVariable('lon', 'f4', ('lon',))

inc = ds1.createVariable('inc', 'f4', ('lat', 'lon',),fill_value=np.nan)
dec = ds1.createVariable('dec', 'f4', ('lat', 'lon',),fill_value=np.nan)

lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = wrf_lat
lons[:] = wrf_lon


for i in range(0,231,1):
    for j in range(0,299,1):
        Rain = ds.variables["RAINFALL"][0,i,j]
        if np.isnan(Rain).any() == True:
            inc[i, j] = np.nan
            dec[i, j] = np.nan
        else:
            RAIN= wrf.variables['RAINC'][0:241,i,j]+wrf.variables['RAINNC'][0:241,i,j]
            RAIN1 = wrf1.variables['RAINC'][0:241, i, j] + wrf1.variables['RAINNC'][0:241, i, j]
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

ds1.close()
