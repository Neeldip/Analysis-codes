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

#FILE_NAME = r'I:\MYD08_M3\2018\APR.hdf'
list=[r'I:\MYD08_M3\2001\APR.hdf',
      r'I:\MYD08_M3\2002\APR.hdf',
      r'I:\MYD08_M3\2003\APR.hdf',
      r'I:\MYD08_M3\2004\APR.hdf',
      r'I:\MYD08_M3\2005\APR.hdf',
      r'I:\MYD08_M3\2006\APR.hdf',
      r'I:\MYD08_M3\2007\APR.hdf',
      r'I:\MYD08_M3\2008\APR.hdf',
      r'I:\MYD08_M3\2009\APR.hdf',
      r'I:\MYD08_M3\2010\APR.hdf',
      r'I:\MYD08_M3\2011\APR.hdf',
      r'I:\MYD08_M3\2012\APR.hdf',
      r'I:\MYD08_M3\2013\APR.hdf',
      r'I:\MYD08_M3\2014\APR.hdf',
      r'I:\MYD08_M3\2015\APR.hdf',
      r'I:\MYD08_M3\2016\APR.hdf',
      r'I:\MYD08_M3\2017\APR.hdf',
      r'I:\MYD08_M3\2018\APR.hdf',
      r'I:\MYD08_M3\2019\APR.hdf',
      r'I:\MYD08_M3\2020\APR.hdf',]

for FILE_NAME in list:

    hdf = SD(FILE_NAME, SDC.READ)

    # List available SDS datasets.
    #print(hdf.datasets())


    # Read dataset with PANDAS
    #DATAFIELD_NAME='Cloud_Optical_Thickness_Liquid_Mean_Mean'
    DATAFIELD_NAME='AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean'
    #DATAFIELD_NAME='Cloud_Fraction_Mean_Mean'
    #DATAFIELD_NAME='Cloud_Effective_Radius_Liquid_Mean_Mean'
    Lat='YDim'
    Lon='XDim'
    AOD = hdf.select(DATAFIELD_NAME)
    lat = hdf.select(Lat)
    lon = hdf.select(Lon)
    data = AOD[:,:]
    print(np.shape(data))
    latitude = lat[:]
    print(latitude)
    longitude = lon[:]
    print(longitude)
    #index starts at upper left corner as (1,1)
    #***********************
    ###LATITUDES AND LONGITUDES ARE OF UPPER LEFT CORNER AND LOWER RIGHT CORNER
    #***********************


    #REGION1
    #sel_lat = 26.5
    #sel_lon = 89.5
    #sel_lat = 25.5
    #sel_lon = 91.5

    #REGION2
    #sel_lat = 28.5
    #sel_lon = 92.5
    #sel_lat = 27.5
    #sel_lon = 96.5

    #REGION3
    #sel_lat = 26.5
    #sel_lon = 92.5
    #sel_lat = 25.5
    #sel_lon = 95.5

    #REGION4
    #sel_lat = 24.5
    #sel_lon = 91.5
    #sel_lat = 22.5
    #sel_lon = 94.5

    #IGP1
    #sel_lat = 26.5
    #sel_lon = 83.5
    #sel_lat = 24.5
    #sel_lon = 88.5

    #IGP2
    #sel_lat = 26.5
    #sel_lon = 83.5
    #sel_lat = 24.5
    #sel_lon = 86.5

    #BELOW IGP
    #sel_lat = 23.5
    #sel_lon = 83.5
    #sel_lat = 19.5
    #sel_lon = 86.5

    sel_lat = 33.5
    sel_lon = 86.5
    latitude=np.reshape(latitude,(180,1))
    #longitude=np.reshape(longitude,(360,1))
    #print(latitude)
    a = abs(latitude-sel_lat)+abs(longitude-sel_lon)
    i,j = np.unravel_index(a.argmin(), a.shape)
    print(i)
    print(j)

    sel_lat = 18.5
    sel_lon = 102.5
    #latitude=np.reshape(latitude,(180,1))
    #longitude=np.reshape(longitude,(360,1))
    #print(latitude)
    #latitude=latitude[::-1]
    #print(latitude)
    #print(longitude)
    c = abs(latitude-sel_lat)+abs(longitude-sel_lon)
    k,l = np.unravel_index(c.argmin(), c.shape)
    print(k)
    print(l)
    ds = np.asarray(data)
    ds = ds.astype('float')
    ds[ds<=-9999]=np.NAN ##replace fill values by NaN
    ds = ds/1000 ##aod
    #ds = ds*0.009999999776482582 #COD & droplet radius
    #ds = ds*9.999999747378752/100000 #cloud fraction
    #print(ds)
    aod = ds[i:k+1,j:l+1]
    #print(aod)
    #print(aod)
    mean = np.nanmean(aod)
    #print(mean)

