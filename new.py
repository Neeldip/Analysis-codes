import numpy as np
from numpy import *
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES,destagger)
import matplotlib.patches as mpatch
#wrf_file = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
#wrf_file1 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrfout_d03.nc")
#wrf_file2 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYJ\wrfout_d03.nc")
#wrf_file3 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrfout_d03.nc")
#wrf_file4 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_YSU\wrfout_d03.nc")
#wrf_file5 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_QNSE\wrfout_d03.nc")
wrf_file = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrf_trad_d03.nc")
wrf_file1 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrf_trad_d03.nc")
wrf_file2 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYJ\wrf_trad_d03.nc")
wrf_file3 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrf_trad_d03.nc")
wrf_file4 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_YSU\wrf_trad_d03.nc")
wrf_file5 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_QNSE\wrf_trad_d03.nc")
lats = wrf_file4.variables['XLAT'][0,:,:] ##231*299
lons = wrf_file4.variables['XLONG'][0,:,:] ##231*299

##SELECT TIME AND LOCATION####
time = ALL_TIMES
if time == ALL_TIMES:
    stagger_dim = 0  ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
else:
    stagger_dim = 0
sel_lat =26.1
sel_lon = 91.74
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)
#####WRF_FILE####################

T = getvar(wrf_file, "QVAPOR", timeidx=time)
PH = getvar(wrf_file, "PH", timeidx=0)
PHB = getvar(wrf_file, "PHB", timeidx=0)
PH = destagger(PH, stagger_dim=stagger_dim,meta=False)  ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)
#print(PH[:,10,10])
if time == ALL_TIMES:
     T = T.mean("Time")

GPH = (PH + PHB) / 9.81
t = T[:, i, j]*1000
y = GPH[:, i, j]

#####WRF_FILE1####################
T = getvar(wrf_file1, "QVAPOR", timeidx=time)
PH = getvar(wrf_file1, "PH", timeidx=0)
PHB = getvar(wrf_file1, "PHB", timeidx=0)
PH = destagger(PH, stagger_dim=stagger_dim,meta=False)  ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)
if time == ALL_TIMES:
     T = T.mean("Time")

GPH = (PH + PHB) / 9.81
t1 = T[:, i, j]*1000
y1 = GPH[:, i, j]

#####WRF_FILE2####################
T = getvar(wrf_file2, "QVAPOR", timeidx=time)
PH = getvar(wrf_file2, "PH", timeidx=0)
PHB = getvar(wrf_file2, "PHB", timeidx=0)
PH = destagger(PH, stagger_dim=stagger_dim,meta=False)  ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)

if time == ALL_TIMES:
     T = T.mean("Time")

GPH = (PH + PHB) / 9.81
t2 = T[:, i, j]*1000
y2 = GPH[:, i, j]

#####WRF_FILE3####################
T = getvar(wrf_file3, "QVAPOR", timeidx=time)
PH = getvar(wrf_file3, "PH", timeidx=0)
PHB = getvar(wrf_file3, "PHB", timeidx=0)
PH = destagger(PH, stagger_dim=stagger_dim,meta=False)  ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)
if time == ALL_TIMES:
     T = T.mean("Time")

GPH = (PH + PHB) / 9.81
t3 = T[:, i, j]*1000
y3 = GPH[:, i, j]

#####WRF_FILE4####################
T = getvar(wrf_file4, "QVAPOR", timeidx=time)
PH = getvar(wrf_file4, "PH", timeidx=0)
PHB = getvar(wrf_file4, "PHB", timeidx=0)
PH = destagger(PH, stagger_dim=stagger_dim,meta=False)  ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)
if time == ALL_TIMES:
     T = T.mean("Time")

GPH = (PH + PHB) / 9.81
t4 = T[:, i, j]*1000
y4 = GPH[:, i, j]
#####WRF_FILE5####################
T = getvar(wrf_file5, "QVAPOR", timeidx=time)
PH = getvar(wrf_file5, "PH", timeidx=0)
PHB = getvar(wrf_file5, "PHB", timeidx=0)
PH = destagger(PH, stagger_dim=stagger_dim,meta=False)  ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)
if time == ALL_TIMES:
     T = T.mean("Time")

GPH = (PH + PHB) / 9.81
t5= T[:, i, j]*1000
y5 = GPH[:, i, j]

plt.figure(figsize=(5.5,10))
plt.plot(t,y,color='dodgerblue')
plt.plot(t1,y1,color='coral')
plt.plot(t2,y2,color='gold')
plt.plot(t3,y3,color='lime')
plt.plot(t4,y4,color='aqua')
plt.plot(t5,y5,color='grey')
plt.ylim(0,5000)
    #plt.legend()
    #plt.xlim(260,300)
#plt.savefig('guwahati.png',dpi=500)
plt.xticks(None)
#plt.axvline(x=0,color="black",linestyle="--")
#magenta = mpatch.Patch(color="magenta",label = "Radiosonde")
blue = mpatch.Patch(color="dodgerblue", label="ACM2")
coral = mpatch.Patch(color="coral", label="HONG")
gold = mpatch.Patch(color="gold", label="MYJ")
chocolate = mpatch.Patch(color="lime", label="MYNN3")
aqua = mpatch.Patch(color="aqua", label="QNSE")
grey = mpatch.Patch(color="grey", label="YSU")
plt.legend(handles=[blue,coral,gold,chocolate,aqua,grey])
plt.ylabel("Height (m)")
plt.xlabel("QVAPOR (g/kg)")
#plt.title("(a) April")
#plt.xlim(275,325)
plt.savefig("QV_APRILavg_profile.png",dpi=300,bbox_inches='tight')
plt.show()


'''

###Plot2 Lucknow
sel_lat =26.76
sel_lon = 80.88
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
x=X[:,i,j]
y=GPH[:,i,j]
plt.plot(x,y)
plt.savefig('lucknow.png',dpi=500)
plt.show()



###Plot3 Kolkata
sel_lat =22.6457
sel_lon = 88.4467
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
x=X[:,i,j]
y=GPH[:,i,j]
plt.plot(x,y)
plt.savefig('kolkata.png',dpi=500)
plt.show()




###Plot4 Bhubaneswar
sel_lat =20.2443
sel_lon =85.8177
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
x=X[:,i,j]
y=GPH[:,i,j]
plt.plot(x,y)
plt.savefig('bhubaneswar.png',dpi=500)
plt.show()



###Plot5 shillong
sel_lat =25.5788
sel_lon =91.8933
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
x=X[:,i,j]
y=GPH[:,i,j]
plt.plot(x,y)
plt.savefig('shillong.png',dpi=500)
plt.show()



###Plot6 agartala
sel_lat =23.8293
sel_lon =91.277
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
x=X[:,i,j]
y=GPH[:,i,j]
plt.plot(x,y)
plt.savefig('agartala.png',dpi=500)
plt.show()



###Plot7 dibrugarh
sel_lat =27.4728
sel_lon =94.9120
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
x=X[:,i,j]
y=GPH[:,i,j]
plt.plot(x,y)
plt.savefig('dibrugarh.png',dpi=500)
plt.show()



###Plot8 imphal
sel_lat =24.8170
sel_lon =93.9368
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
x=X[:,i,j]
y=GPH[:,i,j]
plt.plot(x,y)
plt.savefig('imphal.png',dpi=500)
plt.show()

'''