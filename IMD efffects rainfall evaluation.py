import pandas as pd
from xarray import open_dataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import netCDF4 as nc
from goodness_of_fit import d,rmse,me,r_pearson,mae
from netCDF4 import Dataset
###First regrid the imd dataset with IMD regridding with time axis.py
##cdo remapbil,template.nc 2018-Apr.nc 2018-Apr_regridded.nc

ds = open_dataset(r"F:\WRF-CHEM ANALYSIS Chap 2\WRF-Chem rainfall evaluation\2018-Apr_regridded_bicubic.nc") #<----------------
#wrf =  open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf =  open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf =  open_dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf =  Dataset(r"K:\WRF_Chem_Output\201\April\YSU\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf1 = open_dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf1 = open_dataset(r"G:\WRF_Chem_Output\202\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")

wrf_lat= wrf.variables["XLAT"][0,:,0]
wrf_lon = wrf.variables["XLONG"][0,0,:]

fn = r'F:\WRF-CHEM ANALYSIS Chap 5\rainfall\apr_rainfall_evaluation_NOFEED-NOCHEM.nc'   #<----------------
ds1 = nc.Dataset(fn, 'w', format='NETCDF4')

lat = ds1.createDimension('lat', 231)
lon = ds1.createDimension('lon', 299)
lats = ds1.createVariable('lat', 'f4', ('lat',))
lons = ds1.createVariable('lon', 'f4', ('lon',))
RMSE = ds1.createVariable('rmse', 'f4', ('lat', 'lon',),fill_value=np.nan)
IOA= ds1.createVariable('ioa', 'f4', ('lat', 'lon',),fill_value=np.nan)
ME= ds1.createVariable('mean error', 'f4', ('lat', 'lon',),fill_value=np.nan)
MAE= ds1.createVariable('mae', 'f4', ('lat', 'lon',),fill_value=np.nan)
CC= ds1.createVariable('cc', 'f4', ('lat', 'lon',),fill_value=np.nan)

lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = wrf_lat
lons[:] = wrf_lon

for i in range(0,231,1):
    for j in range(0,299,1):
        Rain = ds.variables["RAINFALL"][9:19,i,j]
        if np.isnan(Rain).any() == True:
           RMSE[i,j] = np.nan
           IOA[i, j] = np.nan
           ME[i, j] = np.nan
           MAE[i, j] = np.nan
           CC[i, j] = np.nan
        else:
            RAIN= wrf.variables['RAINC'][0:241,i,j]+wrf.variables['RAINNC'][0:241,i,j]
            RAIN1= wrf1.variables['RAINC'][0:241,i,j]+wrf1.variables['RAINNC'][0:241,i,j]
            diff= RAIN-RAIN1
            rain= (diff[1:241]-diff[0:240])
            rain_acc_daily=[]
            for time in range(0,240,24):
                rain_daily=rain[time:time+24]
                rain_daily_sum = np.sum(rain_daily)
                rain_acc_daily.append(rain_daily_sum)
                if len(rain_acc_daily)==10:
                    print(i)
                    print(j)
                    RMSE[i, j] = rmse(rain_acc_daily,Rain)
                    IOA[i, j] = d(rain_acc_daily,Rain)
                    ME[i, j] = me(rain_acc_daily,Rain)
                    MAE[i, j] = mae(rain_acc_daily,Rain)
                    CC[i, j] = r_pearson(rain_acc_daily,Rain)
ds1.close()
