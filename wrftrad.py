############THIS DOES NOT WORK#############################

import netCDF4 as nc
import wrf as w
import numpy as np
import scipy as sc
import xarray as xr
import matplotlib as mp
import pandas as pd

# np.set_printoptions(threshold=np.inf)
wrf = nc.Dataset(r"G:\WRF_Outputs\Premonsoon2001_with_nudging_QNSE\wrf_trad_fields_d03_2001-03-31_00%3A00%3A00.nc")
wrf1 = nc.Dataset("G:\WRF_Outputs\Premonsoon2018_with_nudging_ACM2\wrfout_d03_2018-03-31_01%3A00%3A00.nc")
# lat = w.getvar(wrf,'XLAT')
# lon = w.getvar(wrf,'XLONG')

# print(lat,lon)
# print(type(lat))
for timeidx in range(0, 1, 1):
    SPD = w.getvar(wrf, 'SPEED', timeidx)
    SPD1 = SPD[:,:,:]   ###[level,y,x]
    #print(SPD[:,157,132])
    #DIR = w.getvar(wrf, 'DIR', timeidx)
    GHT = w.getvar(wrf, "GEOHEIGHT", timeidx)
    PHB = w.getvar(wrf1,"PHB", timeidx)
    PH = w.getvar(wrf1,"PH", timeidx)
    P1 = w.destagger(PHB, stagger_dim=0, meta=True)
    P2 = w.destagger(PH, stagger_dim=0, meta=True)
    P = (P1+P2)/9.81
    #SPD2 = w.interplevel(field3d=SPD1, vert=GHT, desiredlev=[80,93,180,259,304,321,520,603,770,789,902,988,1497,1537,1968,2010,2094,2461,
    #                                                         2942,3116,3151,3257,3352,3790,3876,4194,4351,4576,5043,5760,5928,6342,6751,7420,
    #                                                         7456,7787,8094,8592,9083,9440,9810,10670,11014,11493,11611,11881,12100,12165,12813,
    #                                                         13140,13254,13764,13930,14055,14403,14719,15152,15826,16400,16570,16638,16769,17278,
    #                                                         17657,17806,18117,18530,18572,18703,18944,19093,19160,19456,19842,20098]
    SPD2 = w.interplevel(field3d=SPD, vert=P, desiredlev=10100, squeeze=False, meta=False)
    SPD2 = SPD.rename({'south_north':'latitude','west_east':"longitude"})
    # SPD = SPD.assign_attrs({'variables': })
    #print(SPD.attrs)
    # SPD = xr.DataArray.expand_dims(SPD,dim=range(0,747,1))
    # sc.interp()
    # SPD = np.interp(SPD, 26.184, 91.693)
    print(SPD2)
    print(SPD2.sizes)
    # xr.plot.pcolormesh(SPD)
    # SPD.rename_dims(SPD, dims_dict=['west_east':'latitude','sou'])
    # SPD.to_netcdf('SPD', mode='w')
    # print(type(INTSPD))
# nc.Dataset("WRF.nc","w", format="NETCDF4")
################
   # VINTP = w.vinterp(wrf,field=SPD, vert_coord='ght_msl', interp_levels=[0.08,0.1], extrapolate=True, timeidx=[1])

    #print(VINTP)
###############
# VINTP1 = w.vinterp(wrf,field=SPD, vert_coord='pressure', interp_levels=[1000])
# print(VINTP1)
