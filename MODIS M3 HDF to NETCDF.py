import netCDF4 as nc
import numpy as np
from pyhdf.SD import SD, SDC
import pandas as pd
#from osgeo import gdal
#https://hdfeos.org/software/pyhdf.php
#https://hdfeos.org/examples/c_grid_lonlat.php

fn = r'F:\AOD_ANALYSIS\APR.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')
ds.set_fill_on()
list=[r'I:\MOD08_M3\2001\APR.hdf',
      r'I:\MOD08_M3\2002\APR.hdf',
      r'I:\MOD08_M3\2003\APR.hdf',
      r'I:\MOD08_M3\2004\APR.hdf',
      r'I:\MOD08_M3\2005\APR.hdf',
      r'I:\MOD08_M3\2006\APR.hdf',
      r'I:\MOD08_M3\2007\APR.hdf',
      r'I:\MOD08_M3\2008\APR.hdf',
      r'I:\MOD08_M3\2009\APR.hdf',
      r'I:\MOD08_M3\2010\APR.hdf',
      r'I:\MOD08_M3\2011\APR.hdf',
      r'I:\MOD08_M3\2012\APR.hdf',
      r'I:\MOD08_M3\2013\APR.hdf',
      r'I:\MOD08_M3\2014\APR.hdf',
      r'I:\MOD08_M3\2015\APR.hdf',
      r'I:\MOD08_M3\2016\APR.hdf',
      r'I:\MOD08_M3\2017\APR.hdf',
      r'I:\MOD08_M3\2018\APR.hdf',
      r'I:\MOD08_M3\2019\APR.hdf',
      r'I:\MOD08_M3\2020\APR.hdf',]
hdf = SD(r'I:\MOD08_M3\2001\APR.hdf', SDC.READ)

#Lat = 'YDim'
#Lon = 'XDim'

LAT = hdf.select('YDim')
LON = hdf.select('XDim')


latitude = LAT[:]
latitude = latitude[::-1] #reverse the list
longitude = LON[:]

lat = ds.createDimension('lat', 180)
lon = ds.createDimension('lon', 360)
time = ds.createDimension('time', 20)

lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
Time = ds.createVariable('time','f4', ('time',))

AOD = ds.createVariable('AOD', 'f4', ('time', 'lat','lon'),fill_value=np.nan)
AOD.units = 'AEROSOL OPTICAL DEPTH'
Cer = ds.createVariable('CER', 'f4', ('time', 'lat','lon'),fill_value=np.nan)
Cer.units = 'CLOUD EFFECTIVE LIQUID DROPLET RADIUS (um)'
Cldfra = ds.createVariable('CLDFRA', 'f4', ('time', 'lat','lon'),fill_value=np.nan)
Cldfra.units = 'CLOUD FRACTION (0 to 1)'
COD = ds.createVariable('COD', 'f4', ('time', 'lat','lon'),fill_value=np.nan)
COD.units = 'CLOUD OPTICAL DEPTH'
ClwP = ds.createVariable('CLWP', 'f4', ('time', 'lat','lon'),fill_value=np.nan)
ClwP.units = 'CLOUD LIQUID WATER PATH g/m2'


lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = latitude
lons[:] = longitude
Time[:] = ['2001','2002','2003','2004','2005','2006','2007','2008','2009','2010',
           '2011','2012','2013','2014','2015','2016','2017','2018','2019','2020']


# index starts at upper left corner as (1,1)
# ***********************
###LATITUDES AND LONGITUDES ARE OF UPPER LEFT CORNER AND LOWER RIGHT CORNER
# ***********************

TIME=np.arange(0,20,1)


for FILE_NAME,T in zip(list,TIME):
    hdf = SD(FILE_NAME, SDC.READ)

    DATAFIELD_NAME = 'AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean'
    AOD_550 = hdf.select(DATAFIELD_NAME)
    data = AOD_550[:, :]
    ds1 = np.asarray(data)
    ds1 = ds1.astype('float')
    ds1[ds1 <= -9999] = np.NAN  ##replace fill values by NaN
    ds1 = ds1 / 1000  ##aod

    DATAFIELD_NAME1='Cloud_Optical_Thickness_Liquid_Mean_Mean'
    COT = hdf.select(DATAFIELD_NAME1)
    data1 = COT[:, :]
    ds2 = np.asarray(data1)
    ds2 = ds2.astype('float')
    ds2[ds2 <= -9999] = np.NAN  ##replace fill values by NaN
    ds2 = ds2*0.009999999776482582  ##aod

    DATAFIELD_NAME2="Cloud_Water_Path_Liquid_Mean_Mean"
    CLWP = hdf.select(DATAFIELD_NAME2)
    data2 = CLWP[:, :]
    ds3 = np.asarray(data2)
    ds3 = ds3.astype('float')
    ds3[ds3 <= -9999] = np.NAN  ##replace fill values by NaN
    ds3 = ds3

    DATAFIELD_NAME3='Cloud_Fraction_Mean_Mean'
    CLDFRA = hdf.select(DATAFIELD_NAME3)
    data3 = CLDFRA[:, :]
    ds4 = np.asarray(data3)
    ds4 = ds4.astype('float')
    ds4[ds4 <= -9999] = np.NAN  ##replace fill values by NaN
    ds4 = ds4*9.999999747378752/100000

    DATAFIELD_NAME4 = 'Cloud_Effective_Radius_Liquid_Mean_Mean'
    CER = hdf.select(DATAFIELD_NAME4)
    data4 = CER[:, :]
    ds5 = np.asarray(data4)
    ds5 = ds5.astype('float')
    ds5[ds5 <= -9999] = np.NAN  ##replace fill values by NaN
    ds5 = ds5 * 0.009999999776482582  ##aod

    # ds1 = ds1*0.009999999776482582 #COD & droplet radius
    # ds1 = ds1*9.999999747378752/100000 #cloud fraction
    #print(ds)
    #for i,l in zip(range(179,-1,-1),range(0,180,1)):

    for i,l in zip(range(179,-1,-1),range(0,180,1)): ###starts from lower left corner
        for j in range(0,360,1):
            AOD[T,l,j]= ds1[i,j]
            COD[T, l, j] = ds2[i, j]
            ClwP[T, l, j] = ds3[i, j]
            Cldfra[T, l, j] = ds4[i, j]
            Cer[T, l, j] = ds5[i, j]

    print("###")

ds.close()