import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import h5py
from pyhdf.SD import SD, SDC
import pandas as pd
from osgeo import gdal
#https://hdfeos.org/software/pyhdf.php
#https://hdfeos.org/examples/c_grid_lonlat.php
#'''

#FILE_NAME = r'I:\MYD08_M3\2018\MAR.hdf'
list=[r'I:\MYD08_M3\2001\JAN.hdf',
      r'I:\MYD08_M3\2002\JAN.hdf',
      r'I:\MYD08_M3\2003\JAN.hdf',
      r'I:\MYD08_M3\2004\JAN.hdf',
      r'I:\MYD08_M3\2005\JAN.hdf',
      r'I:\MYD08_M3\2006\JAN.hdf',
      r'I:\MYD08_M3\2007\JAN.hdf',
      r'I:\MYD08_M3\2008\JAN.hdf',
      r'I:\MYD08_M3\2009\JAN.hdf',
      r'I:\MYD08_M3\2010\JAN.hdf',
      r'I:\MYD08_M3\2011\JAN.hdf',
      r'I:\MYD08_M3\2012\JAN.hdf',
      r'I:\MYD08_M3\2013\JAN.hdf',
      r'I:\MYD08_M3\2014\JAN.hdf',
      r'I:\MYD08_M3\2015\JAN.hdf',
      r'I:\MYD08_M3\2016\JAN.hdf',
      r'I:\MYD08_M3\2017\JAN.hdf',
      r'I:\MYD08_M3\2018\JAN.hdf',
      r'I:\MYD08_M3\2019\JAN.hdf',
      r'I:\MYD08_M3\2020\JAN.hdf',]

for FILE_NAME in list:

    hdf = SD(FILE_NAME, SDC.READ)

    # List available SDS datasets.
    #print(hdf.datasets())


    # Read dataset with PANDAS
    DATAFIELD_NAME='AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean'
    Lat='YDim'
    Lon='XDim'
    AOD = hdf.select(DATAFIELD_NAME)
    lat = hdf.select(Lat)
    lon = hdf.select(Lon)
    data = AOD[:,:]
    #print(data)
    ds = pd.DataFrame(data)
    data1 = ds.replace(-9999,np.NaN)
    data2 = data1/1000
    #ds.to_csv("AOD.csv")

    # Read geolocation dataset and plotting with pandas
    latitude = lat[:]
    longitude = lon[:]
    m = Basemap(projection='cyl', resolution='h', llcrnrlat=5, urcrnrlat = 40, llcrnrlon=60, urcrnrlon = 100)
    m.drawcoastlines(linewidth=0.5)
    #m.drawparallels(np.arange(-90., 120., 30.), labels=[1, 0, 0, 0])
    #m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1])
    x, y = m(longitude, latitude)
    #m.pcolormesh(x, y, data,shading='auto')
    v = np.linspace(0,1,20)
    plt.contourf(longitude,latitude,data2,v,color="jet",extend='max')
    plt.colorbar(ticks=v)
    #plt.show()

    ###extarction at a location cant be done with pandas
    #so numpy is used
    #Read geolocation dataset and plotting with numpy
    sel_lat = 26.11
    sel_lon = 91.74
    latitude=np.reshape(latitude,(180,1))
    #longitude=np.reshape(longitude,(360,1))
    #print(latitude)
    a = abs(latitude-sel_lat)+abs(longitude-sel_lon)
    #print(a)
    i,j = np.unravel_index(a.argmin(), a.shape)
    #print(i)
    #print(j)
    ds = np.asarray(data)
    ds = ds.astype('float')
    ds[ds<=-9999]=np.NAN ##replace fill values by NaN
    ds = ds/1000
    #print(ds)
    aod = ds[i,j]
    print(aod)
    #******************PLOTTING with NUMPY***************************
    #m = Basemap(projection='cyl', resolution='h', llcrnrlat=5, urcrnrlat = 40, llcrnrlon=60, urcrnrlon = 100)
    #m.drawcoastlines(linewidth=0.5)
    #plt.contourf(longitude,latitude,ds,v,color="jet",extend='max')
    #plt.colorbar(ticks=v)
    #plt.show()


    #'''

    #file = gdal.Open(r'I:\MYD08_M3\2018\APR.hdf')
    #for path, desc in file.GetSubDatasets():
    #    print(desc)
