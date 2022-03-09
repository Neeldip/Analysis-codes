from xarray import open_dataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from goodness_of_fit import mae,d,rmse,me
import netCDF4 as nc
ds = open_dataset(r"F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\apr_rainfall_evaluation_nofeedback_new.nc")
wrf_nor =  open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_nofeed =  open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")

#This program finds whether overall rainfall decreased or increased dut to aerosols at locations where performance passes criteria
#rainfall statistics is read from 'ds'. only values fulfilling the crietria is only considered an mean of those values is calculated

rmse = ds.variables["rmse"][123:231, 182:299]
ioa = ds.variables["ioa"][123:231, 182:299]
perc=50  ##percentile value input
ioa_perc= np.nanpercentile(ioa,perc)    ##ioa percentile calculation
rmse_perc= np.nanpercentile(rmse,100-perc) ## since lower rmse values are better, 100-perc
#print(ioa_perc)
#print(rmse_perc)
array1=np.empty((231,299))  #empty array of 231*299
wrf_lat= wrf_nor.variables["XLAT"][0,:,0]
wrf_lon = wrf_nor.variables["XLONG"][0,0,:]

fn = r'F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\apr_ME_evaluation_50_percentile_nofeedback.nc'   #<----------------
ds1 = nc.Dataset(fn, 'w', format='NETCDF4')

lat = ds1.createDimension('lat', 231)
lon = ds1.createDimension('lon', 299)
lats = ds1.createVariable('lat', 'f4', ('lat',))
lons = ds1.createVariable('lon', 'f4', ('lon',))
array=ds1.createVariable('mean error', 'f4', ('lat', 'lon',),fill_value=np.nan)

lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = wrf_lat
lons[:] = wrf_lon


for i in range(0,231,1):
    for j in range(0,299,1):
        IOA = ds.variables["ioa"][i,j]
        RMSE= ds.variables["rmse"][i,j]
        IOA=np.asarray(IOA)
        RMSE=np.asarray(RMSE)
        IOA=np.asscalar(IOA)
        RMSE=np.asscalar(RMSE)
        if IOA==np.nan and RMSE==np.nan:
           array[i,j] = np.nan
           array1[i, j] = np.nan
        elif IOA<ioa_perc and RMSE>rmse_perc:   ##location selection is then done based on these percentile values
            array[i,j]=np.nan
            array1[i, j] = np.nan
        elif IOA<ioa_perc and RMSE<rmse_perc:
            array[i,j]=np.nan
            array1[i, j] = np.nan
        elif IOA>ioa_perc and RMSE>rmse_perc:
            array[i,j]=np.nan
            array1[i, j] = np.nan
        elif IOA==0 and RMSE==0:
            array[i,j]=np.nan
            array1[i, j] = np.nan
        elif IOA>=ioa_perc and RMSE <=rmse_perc:
            RAIN_nor= wrf_nor.variables['RAINC'][0:241,i,j]+wrf_nor.variables['RAINNC'][0:241,i,j]
            RAIN_nofeed = wrf_nofeed.variables['RAINC'][0:241, i, j] + wrf_nofeed.variables['RAINNC'][0:241, i, j]
            rain_nor= (RAIN_nor[1:241]-RAIN_nor[0:240])
            rain_nofeed = (RAIN_nofeed[1:241] - RAIN_nofeed[0:240])
            rain_nor_acc_daily=[]
            rain_nofeed_acc_daily = []
            for time in range(0,240,24):
                rain_nor_daily=rain_nor[time:time+24]
                rain_nor_daily_sum = np.sum(rain_nor_daily)
                rain_nor_acc_daily.append(rain_nor_daily_sum)
                rain_nofeed_daily=rain_nofeed[time:time+24]
                rain_nofeed_daily_sum = np.sum(rain_nofeed_daily)
                rain_nofeed_acc_daily.append(rain_nofeed_daily_sum)
                if len(rain_nor_acc_daily)==10:
                    array[i, j] = me(rain_nor_acc_daily, rain_nofeed_acc_daily)
                    array1[i, j] = me(rain_nor_acc_daily, rain_nofeed_acc_daily)
        else:
            array[i, j] = np.nan
            array1[i, j] = np.nan

for i in range(0, 231, 1):         #convert to nan at locations outside NE region
    for j in range(0, 182, 1):
        array[i, j] = np.nan
        array1[i, j] = np.nan
for i in range(0, 123, 1):
    for j in range(182, 273, 1):
        array[i, j] = np.nan
        array1[i, j] = np.nan
for i in range(215, 231, 1):
   for j in range(182, 273, 1):
        array[i, j] = np.nan
        array1[i, j] = np.nan
for i in range(0, 231, 1):
  for j in range(273, 299, 1):
        array[i, j] = np.nan
        array1[i, j] = np.nan

a=array1[np.where(array1>0)]
b=array1[np.where(array1==0)]
c=array1[np.where(array1<0)]
print(np.size(a))
#print(np.size(b))
print(np.size(c))
#print(np.count_nonzero(np.isnan(array)))
print(np.nanmean(array1)) #unit mm/day per grid point

ds1.close()

