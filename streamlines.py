import numpy as np
from numpy import *
from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES, enable_pyngl)
import netCDF4 as nc
#wrf_file = [Dataset(r"H:\wrfout_d01_2018-03-11_00%3A00%3A00"),
#            Dataset(r"H:\wrfout_d01_2018-03-12_00%3A00%3A00"),
#            Dataset(r"H:\wrfout_d01_2018-03-13_00%3A00%3A00"),
#            Dataset(r"H:\wrfout_d01_2018-03-14_00%3A00%3A00"),
#            Dataset(r"H:\wrfout_d01_2018-03-15_00%3A00%3A00")]

wrf_file = nc.Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
time = ALL_TIMES
#time = 1
speed_plot =True
level = [0,5,10,13,15,20,25,30]
U = getvar(wrf_file,"ua",timeidx=time)
V = getvar(wrf_file,"va",timeidx=time)
if time == ALL_TIMES:
     U = U.mean("Time")
     V = V.mean("Time")

x = np.linspace(71.679886,100.43091,299)
y = np.linspace(10.64794,31.228256,231)
#x = np.linspace(67.5047,102.5753,359)
#y = np.linspace(10.82174,34.78333,269)
#x = np.linspace(87.65714,98.07542,348)
#y = np.linspace(21.82169,30.20221,312)


for t in level:
    #u = U[:,:]
    #v = V[:,:]
    u = U[t,:,:]
    v = V[t,:,:]
    if speed_plot ==True:
     SPD = sqrt(u**2+v**2)
     SPD=asarray(SPD)
     fig = pyplot.figure(figsize=(9, 7))
     k = np.linspace(np.min(SPD), np.max(SPD), 21)
     strm=pyplot.streamplot(x,y,u,v,density=4,arrowsize=1,integration_direction='forward',color=SPD,cmap=pyplot.get_cmap("jet"))
     fig.colorbar(strm.lines,ticks=k)
    else:
        pyplot.streamplot(x, y, u, v, density=8, arrowsize=1, integration_direction='forward',cmap=pyplot.get_cmap("jet"))

    #pyplot.savefig("D:\PHD\My PhD Reports\Manuscript3\IMAGES\STREAMLINES\STREAMLINES AT level=" + str(t) + ".png",dpi=600, bbox_inches='tight')
    pyplot.show()

