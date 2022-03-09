import pandas as pd
from xarray import open_dataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import netCDF4 as nc
from goodness_of_fit import mae,d,rmse,me

###First regrid the imd dataset with IMD regridding with time axis.py
##cdo remapbil,template.nc 2018-Apr.nc 2018-Apr_regridded.nc

ds = open_dataset(r"F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\2018-Apr_regridded_bicubic.nc") #<----------------
#wrf =  open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf =  open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_lat= wrf.variables["XLAT"][0,:,0]
wrf_lon = wrf.variables["XLONG"][0,0,:]
lats= ds.variables["lat"][:]
lons=ds.variables['lon'][:]

sel_lat =25.47598
sel_lon =94.2562
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)
#print(i)
#print(j)

Rain = ds.variables["RAINFALL"][9:19,i,j]
#print(np.reshape(np.asarray(Rain), (10, 1)))
RAINC=wrf.variables['RAINC'][0:241,i,j]
RAINNC=wrf.variables['RAINNC'][0:241,i,j]
#RAIN= wrf.variables['RAINC'][0:241,i,j]+wrf.variables['RAINNC'][0:241,i,j]
RAIN=RAINC+RAINNC
rain= (RAIN[1:241]-RAIN[0:240])
rainc=(RAINC[1:241]-RAINC[0:240])
rainnc=(RAINNC[1:241]-RAINNC[0:240])
rain_acc_daily=[]
rainc_acc_daily=[]
rainnc_acc_daily=[]
for time in range(0,240,24):
        rain_daily=rain[time:time+24]
        rain_daily_sum = np.sum(rain_daily)
        rain_acc_daily.append(rain_daily_sum)
        rainc_daily=rainc[time:time+24]
        rainc_daily_sum = np.sum(rainc_daily)
        rainc_acc_daily.append(rainc_daily_sum)
        rainnc_daily=rainnc[time:time+24]
        rainnc_daily_sum = np.sum(rainnc_daily)
        rainnc_acc_daily.append(rainnc_daily_sum)
        if len(rain_acc_daily)==10:
            #print(np.reshape(np.asarray(rain_acc_daily), (10, 1)))
            print(np.reshape(np.asarray(rainc_acc_daily),(10,1)))
            print(np.reshape(np.asarray(rainnc_acc_daily), (10, 1)))
            RMSE = rmse(rain_acc_daily,Rain)
            IOA = d(rain_acc_daily,Rain)
            ME = me(rain_acc_daily,Rain)
            print(RMSE)
            print(IOA)
            print(ME)


