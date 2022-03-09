import numpy as np
from numpy import *

from matplotlib import pyplot,colors, axes,gridspec,figure
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES, enable_pyngl)


wrf_file = Dataset(r"H:\WRF_Chem_Output\201\May\wrfout_d01_2018-05-10_00%3A00%3A00")
trad = Dataset(r"H:\WRF_Chem_Output\201\may\wrf_trad_fields_d01_2018-05-10_00%3A00%3A00")
#LAT= wrf_file.variables["XLAT"][0,:,0]
#LON= wrf_file.variables["XLONG"][0,0,:]

time = ALL_TIMES
#time = 0
level = [10]
#level = [0,5,10,15,20,25,30]
#0=20m
#5=160m
#10=520m
#15=1152m
#20=2337m
#25=4851m
#30=8111m
#32=10497m
###########################################BC#######################################
ua = getvar(wrf_file,"ua",timeidx=time)
va = getvar(wrf_file,"va",timeidx=time)
BIN1 = getvar(wrf_file, "bc_a01", timeidx=time)
BIN2 = getvar(wrf_file, "bc_a02", timeidx=time)
BIN3 = getvar(wrf_file, "bc_a03", timeidx=time)
BIN4 = getvar(wrf_file, "bc_a04", timeidx=time)

T = getvar(trad, "TEMPERATURE", timeidx=time)  ##sensible temperature,units = K
P = getvar(trad, "PRESSURE", timeidx=time)     ##units = Pascal
X = BIN1 + BIN2 + BIN3 + BIN4

BIN1=None
BIN2=None
BIN3=None
BIN4=None

if time == ALL_TIMES:
     ua = ua.mean("Time")
     va = va.mean("Time")
     X = X.mean("Time")
     T = T.mean("Time")
     P = P.mean("Time")
print("##")

R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1
#x=LON
#y=LAT
x = np.linspace(71.679886,100.43091,299)
y = np.linspace(10.64794,31.228256,231)
#x = np.linspace(np.min(lon),np.max(lon),359)
#y = np.linspace(np.min(lat),np.max(lat),269)
height =[520]
#height =[20,160,520,1150,2340,4850,8100]

for t,h in zip(level,height):
    U = ua[t,:,:]
    V = va[t,:,:]
    s = X[t,:,:]
    a = T[t,:,:]
    b = P[t,:,:]
    K = (R * a) / b
    f = s/K
    u = U*f
    v = V*f
    SPD = sqrt(u**2+v**2)
    SPD = asarray(SPD)
    fig = pyplot.figure(figsize=(9, 7))
    k = np.linspace(np.min(SPD), np.max(SPD), 15, endpoint=True)
    #k = np.linspace(0, 30, 20, endpoint=True)
    strm=pyplot.streamplot(x, y, U, V, density=6, arrowsize=1, integration_direction='forward', color=SPD, cmap="jet",linewidth=1.2)
    #pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
    fig.colorbar(strm.lines,ticks=k)
    #fig.colorbar(ticks=k)
    pyplot.ylabel("Latitude")
    pyplot.xlabel("Longitude")
    pyplot.title("October - BC aerosol flux at "+str(h)+" m")
    #pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\2Dfluxstreamline_October_BC1 at "+str(h)+" m.png",dpi=300,bbox_inches="tight")
    pyplot.show()


s=None
a=None
b=None
K=None
f=None
u=None
v=None
X=None
U=None
V=None

'''
###########################################OC#######################################

BIN1 = getvar(wrf_file, "oc_a01", timeidx=time)
BIN2 = getvar(wrf_file, "oc_a02", timeidx=time)
BIN3 = getvar(wrf_file, "oc_a03", timeidx=time)
BIN4 = getvar(wrf_file, "oc_a04", timeidx=time)

X = BIN1 + BIN2 + BIN3 + BIN4

BIN1=None
BIN2=None
BIN3=None
BIN4=None

if time == ALL_TIMES:
     X = X.mean("Time")

print("##")
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1

for t,h in zip(level,height):
    U = ua[t,:,:]
    V = va[t,:,:]
    s = X[t,:,:]
    a = T[t,:,:]
    b = P[t,:,:]
    K = (R * a) / b
    f = s/K
    u = U*f
    v = V*f
    SPD = sqrt(u**2+v**2)
    SPD = asarray(SPD)
    fig = pyplot.figure(figsize=(9, 7))
    k = np.linspace(np.min(SPD), np.max(SPD), 15, endpoint=True)
    #k = np.linspace(0, 30, 20, endpoint=True)
    strm=pyplot.streamplot(x, y, U, V, density=6, arrowsize=1, integration_direction='forward', color=SPD, cmap="jet",linewidth=1.2)
    #pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
    fig.colorbar(strm.lines,ticks=k)
    #fig.colorbar(ticks=k)
    pyplot.ylabel("Latitude")
    pyplot.xlabel("Longitude")
    pyplot.title("October - OC aerosol flux at "+str(h)+" m")
    pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\2Dfluxstreamline_October_OC at "+str(h)+" m.png",dpi=300,bbox_inches="tight")
    #pyplot.show()


s=None
a=None
b=None
K=None
f=None
u=None
v=None
X=None
U=None
V=None

###########################################OIN#######################################

BIN1 = getvar(wrf_file, "oin_a01", timeidx=time)
BIN2 = getvar(wrf_file, "oin_a02", timeidx=time)
BIN3 = getvar(wrf_file, "oin_a03", timeidx=time)
BIN4 = getvar(wrf_file, "oin_a04", timeidx=time)

X = BIN1 + BIN2 + BIN3 + BIN4

BIN1=None
BIN2=None
BIN3=None
BIN4=None

if time == ALL_TIMES:
     X = X.mean("Time")

print("##")
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1

for t,h in zip(level,height):
    U = ua[t,:,:]
    V = va[t,:,:]
    s = X[t,:,:]
    a = T[t,:,:]
    b = P[t,:,:]
    K = (R * a) / b
    f = s/K
    u = U*f
    v = V*f
    SPD = sqrt(u**2+v**2)
    SPD = asarray(SPD)
    fig = pyplot.figure(figsize=(9, 7))
    k = np.linspace(np.min(SPD), np.max(SPD), 15, endpoint=True)
    #k = np.linspace(0, 30, 20, endpoint=True)
    strm=pyplot.streamplot(x, y, U, V, density=6, arrowsize=1, integration_direction='forward', color=SPD, cmap="jet",linewidth=1.2)
    #pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
    fig.colorbar(strm.lines,ticks=k)
    #fig.colorbar(ticks=k)
    pyplot.ylabel("Latitude")
    pyplot.xlabel("Longitude")
    pyplot.title("October - Dust aerosol flux at "+str(h)+" m")
    pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\2Dfluxstreamline_October_Dust at "+str(h)+" m.png",dpi=300,bbox_inches="tight")
    #pyplot.show()

s=None
a=None
b=None
K=None
f=None
u=None
v=None
X=None
U=None
V=None

###########################################so4#######################################

BIN1 = getvar(wrf_file, "so4_a01", timeidx=time)
BIN2 = getvar(wrf_file, "so4_a02", timeidx=time)
BIN3 = getvar(wrf_file, "so4_a03", timeidx=time)
BIN4 = getvar(wrf_file, "so4_a04", timeidx=time)

X = BIN1 + BIN2 + BIN3 + BIN4

BIN1=None
BIN2=None
BIN3=None
BIN4=None

if time == ALL_TIMES:
     X = X.mean("Time")

print("##")
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1

for t,h in zip(level,height):
    U = ua[t,:,:]
    V = va[t,:,:]
    s = X[t,:,:]
    a = T[t,:,:]
    b = P[t,:,:]
    K = (R * a) / b
    f = s/K
    u = U*f
    v = V*f
    SPD = sqrt(u**2+v**2)
    SPD = asarray(SPD)
    fig = pyplot.figure(figsize=(9, 7))
    k = np.linspace(np.min(SPD), np.max(SPD), 15, endpoint=True)
    #k = np.linspace(0, 30, 20, endpoint=True)
    strm=pyplot.streamplot(x, y, U, V, density=6, arrowsize=1, integration_direction='forward', color=SPD, cmap="jet",linewidth=1.2)
    #pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
    fig.colorbar(strm.lines,ticks=k)
    #fig.colorbar(ticks=k)
    pyplot.ylabel("Latitude")
    pyplot.xlabel("Longitude")
    pyplot.title("October - Sulfate aerosol flux at "+str(h)+" m")
    pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\2Dfluxstreamline_October_Sulfate at "+str(h)+" m.png",dpi=300,bbox_inches="tight")
    #pyplot.show()

s=None
a=None
b=None
K=None
f=None
u=None
v=None
X=None
U=None
V=None

###########################################no3#######################################

BIN1 = getvar(wrf_file, "no3_a01", timeidx=time)
BIN2 = getvar(wrf_file, "no3_a02", timeidx=time)
BIN3 = getvar(wrf_file, "no3_a03", timeidx=time)
BIN4 = getvar(wrf_file, "no3_a04", timeidx=time)

X = BIN1 + BIN2 + BIN3 + BIN4

BIN1=None
BIN2=None
BIN3=None
BIN4=None

if time == ALL_TIMES:
     X = X.mean("Time")

print("##")
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1

for t,h in zip(level,height):
    U = ua[t,:,:]
    V = va[t,:,:]
    s = X[t,:,:]
    a = T[t,:,:]
    b = P[t,:,:]
    K = (R * a) / b
    f = s/K
    u = U*f
    v = V*f
    SPD = sqrt(u**2+v**2)
    SPD = asarray(SPD)
    fig = pyplot.figure(figsize=(9, 7))
    k = np.linspace(np.min(SPD), np.max(SPD), 15, endpoint=True)
    #k = np.linspace(0, 30, 20, endpoint=True)
    strm=pyplot.streamplot(x, y, U, V, density=6, arrowsize=1, integration_direction='forward', color=SPD, cmap="jet",linewidth=1.2)
    #pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
    fig.colorbar(strm.lines,ticks=k)
    #fig.colorbar(ticks=k)
    pyplot.ylabel("Latitude")
    pyplot.xlabel("Longitude")
    pyplot.title("October - Nitrate aerosol flux at "+str(h)+" m")
    pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\2Dfluxstreamline_October_Nitrate at "+str(h)+" m.png",dpi=300,bbox_inches="tight")
    #pyplot.show()


s=None
a=None
b=None
K=None
f=None
u=None
v=None
X=None
U=None
V=None

###########################################nh4#######################################

BIN1 = getvar(wrf_file, "nh4_a01", timeidx=time)
BIN2 = getvar(wrf_file, "nh4_a02", timeidx=time)
BIN3 = getvar(wrf_file, "nh4_a03", timeidx=time)
BIN4 = getvar(wrf_file, "nh4_a04", timeidx=time)

X = BIN1 + BIN2 + BIN3 + BIN4

BIN1=None
BIN2=None
BIN3=None
BIN4=None

if time == ALL_TIMES:
     X = X.mean("Time")

print("##")
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1

for t,h in zip(level,height):
    U = ua[t,:,:]
    V = va[t,:,:]
    s = X[t,:,:]
    a = T[t,:,:]
    b = P[t,:,:]
    K = (R * a) / b
    f = s/K
    u = U*f
    v = V*f
    SPD = sqrt(u**2+v**2)
    SPD = asarray(SPD)
    fig = pyplot.figure(figsize=(9, 7))
    k = np.linspace(np.min(SPD), np.max(SPD), 15, endpoint=True)
    #k = np.linspace(0, 30, 20, endpoint=True)
    strm=pyplot.streamplot(x, y, U, V, density=6, arrowsize=1, integration_direction='forward', color=SPD, cmap="jet",linewidth=1.2)
    #pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
    fig.colorbar(strm.lines,ticks=k)
    #fig.colorbar(ticks=k)
    pyplot.ylabel("Latitude")
    pyplot.xlabel("Longitude")
    pyplot.title("October - Ammonium aerosol flux at "+str(h)+" m")
    pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\2Dfluxstreamline_October_Ammonium at "+str(h)+" m.png",dpi=300,bbox_inches="tight")
    #pyplot.show()

s=None
a=None
b=None
K=None
f=None
u=None
v=None
X=None
U=None
V=None
T=None
P=None


###########################################pm2.5#######################################

BIN1 = getvar(wrf_file, "PM2_5_DRY", timeidx=time)
#BIN2 = getvar(wrf_file, "nh4_a02", timeidx=time)
#BIN3 = getvar(wrf_file, "nh4_a03", timeidx=time)
#BIN4 = getvar(wrf_file, "nh4_a04", timeidx=time)

X = BIN1 #+ BIN2 + BIN3 + BIN4

BIN1=None
#BIN2=None
#BIN3=None
#BIN4=None

if time == ALL_TIMES:
     X = X.mean("Time")

print("##")
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1

for t,h in zip(level,height):
    U = ua[t,:,:]
    V = va[t,:,:]
    s = X[t,:,:]
    f = s
    u = U*f
    v = V*f
    SPD = sqrt(u**2+v**2)
    SPD = asarray(SPD)
    fig = pyplot.figure(figsize=(9, 7))
    k = np.linspace(np.min(SPD), np.max(SPD), 15, endpoint=True)
    #k = np.linspace(0, 30, 20, endpoint=True)
    strm=pyplot.streamplot(x, y, U, V, density=6, arrowsize=1, integration_direction='forward', color=SPD, cmap="jet",linewidth=1.2)
    #pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
    fig.colorbar(strm.lines,ticks=k)
    #fig.colorbar(ticks=k)
    pyplot.ylabel("Latitude")
    pyplot.xlabel("Longitude")
    pyplot.title("October - PM2.5 aerosol flux at "+str(h)+" m")
    pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\2Dfluxstreamline_October_PM2.5 at "+str(h)+" m.png",dpi=300,bbox_inches="tight")
    #pyplot.show()

s=None
a=None
b=None
K=None
f=None
u=None
v=None
X=None
U=None
V=None
ua=None
va=None
'''
