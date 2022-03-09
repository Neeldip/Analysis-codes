import numpy as np
from numpy import *

from matplotlib import pyplot,colors, axes,gridspec,figure
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES, enable_pyngl)

#wrf_file = Dataset(r"G:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"G:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
##wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\July\wrfout_d01_2018-07-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\July\wrf_trad_fields_d01_2018-07-10_00%3A00%3A00")
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")

trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")


time = ALL_TIMES
level = [0,5,10,13,15,20,25,30]
VAR = wrf_file.variables["QVAPOR"][:,level,:,:]
#print(VAR)
ua = getvar(wrf_file,"ua",timeidx=time)
print("ua complete")
va = getvar(wrf_file,"va",timeidx=time)
print("va complete")
#VAR = getvar(wrf_file, "QVAPOR", timeidx=time)

#T = getvar(trad, "TEMPERATURE", timeidx=time)  ##sensible temperature,units = K
T = trad.variables["TEMPERATURE"][:,level,:,:]
#P = getvar(trad, "PRESSURE", timeidx=time)     ##units = Pascal
P = trad.variables["PRESSURE"][:,level,:,:]
X = VAR

R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1
x = np.linspace(71.679886,100.43091,299)
y = np.linspace(10.64794,31.228256,231)
#x = np.linspace(87.65714,98.07542,348)
#y = np.linspace(21.82169,30.20221,312)
print("reading complete")


UA=ua[:,level,:,:]
VA=va[:,level,:,:]
UA = UA.mean("Time")
VA = VA.mean("Time")
G = X.mean(axis=0)
H = T.mean(axis=0)
I = P.mean(axis=0)


U = UA
V = VA
#s = G.squeeze(axis=0)
#A = H.squeeze(axis=0)
#B = I.squeeze(axis=0)
K = (R * H) / I
L = G/K
u = U#*L
v = V#*L
#u = u.squeeze(axis=0)
#v = v.squeeze(axis=0)
SPD = sqrt(u**2+v**2)
SPD = asarray(SPD)
#SPD = SPD.squeeze(axis=0)
print(np.shape(u))
print(np.shape(SPD))
fig = pyplot.figure(figsize=(9, 7))
k = np.linspace(np.min(SPD), np.max(SPD), 21, endpoint=True)
#k = np.linspace(np.min(SPD), 50, 30, endpoint=True)
strm=pyplot.streamplot(x, y, u, v, density=7, arrowsize=0.8, integration_direction='forward', color=SPD, cmap="nipy_spectral",linewidth=1.5)
#pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
fig.colorbar(strm.lines,ticks=k)
pyplot.savefig("D:\PHD\My PhD Reports\Manuscript3\IMAGES\STREAMLINES\STREAMLINES AT level="+str(level)+"time="+str(time)+".png",dpi=600, bbox_inches='tight')
pyplot.show()
g=None
H=None
I=None
X=None
T=None
P=None
U=None
V=None
UA=None
VA=None
s=None
A=None
B=None
K=None
L=None
SPD=None