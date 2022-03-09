import netCDF4 as nc
import wrf as w
import numpy as np
import math as m
import scipy as sc
import xarray as xr
import matplotlib as mp
import pandas as pd
#np.set_printoptions(threshold=np.inf)
wrf = nc.Dataset("G:\WRF_Outputs\Premonsoon2018_with_nudging_ACM2\wrfout_d03_2018-03-31_01%3A00%3A00.nc")

lat = w.getvar(wrf,'XLAT')
lon = w.getvar(wrf,'XLONG')
x =[93,180,259,304,321,520,603,770,789,902,988,1497,1537,1968,2010,2094,2461,2942,3116,3151,3257,3352,3790,3876,4194,4351,4576,5043,5760,5928,6342,6751,7420,7456,7787,8094,8592,9083,9440,9810,10670,11014,11493,11611,11881,12100,12165,12813,13140,13254,13764,13930,14055,14403,14719,15152,15826,16400,16570,16638,16769,17278,17657,17806,18117,18530,18572,18703,18944,19093,19160,19456,19842]
#print(lat,lon)
#print(type(lat))
for timeidx in range(24,48,1):
    for y in x:
           U = w.getvar(wrf,'ua', timeidx)
           V = w.getvar(wrf,'va', timeidx)
           SPD = U[:,:,:]**2+V[:,:,:]**2
           SPD = np.sqrt(SPD)
           PHB = w.getvar(wrf,"PHB", timeidx)
           PH = w.getvar(wrf,"PH", timeidx)
           P1 = w.destagger(PHB, stagger_dim=0, meta=True)
           P2 = w.destagger(PH, stagger_dim=0, meta=True)
           P = (P1+P2)/9.81
           SPD = w.interplevel(field3d=U, vert=P, desiredlev=y)
           SPD = SPD.rename({'south_north':'latitude','west_east':"longitude"})
           print(SPD[156,131])
           print(DIR[156,131])
    #xr.plot.pcolormesh(SPD)
    #SPD.rename_dims(SPD, dims_dict=['west_east':'latitude','sou'])
    #SPD.to_netcdf('SPD', mode='w')
    #print(type(INTSPD))
#nc.Dataset("WRF.nc","w", format="NETCDF4")
################
