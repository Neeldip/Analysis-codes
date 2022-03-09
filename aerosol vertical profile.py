import numpy as np
from numpy import *
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES,destagger)

##ACM2
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

##MYNN3
#wrf_file = Dataset("G:\WRF_Chem_Output\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")


#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#'''
#lats = wrf_file.variables['XLAT'][0,:,:] ##231*299
#lons = wrf_file.variables['XLONG'][0,:,:] ##231*299
lats = getvar(wrf_file,"XLAT",timeidx=0) ##231*299
lons = getvar(wrf_file,"XLONG",timeidx=0) ##231*299

##SELECT TIME AND LOCATION####
time = ALL_TIMES

if time == ALL_TIMES:
 stagger_dim = 1 ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
else:
    stagger_dim = 0


####
#print(i) = 169 for sel_lat =26.1 sel_lon = 91.7
#print(j) =208
###

BIN1 = getvar(wrf_file, "bc_a01", timeidx=time)
BIN2 = getvar(wrf_file, "bc_a02", timeidx=time)
BIN3 = getvar(wrf_file, "bc_a03", timeidx=time)
BIN4 = getvar(wrf_file, "bc_a04", timeidx=time)

PH = getvar(wrf_file,"PH",timeidx=time)
PHB = getvar(wrf_file,"PHB",timeidx=time)
PH = destagger(PH, stagger_dim=stagger_dim, meta=False)   ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)

T = getvar(trad, "TEMPERATURE", timeidx=time)  ##sensible temperature,units = K
P = getvar(trad, "PRESSURE", timeidx=time)     ##units = Pascal

#print(np.shape(PH)) = 44*231*299
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1
X = BIN1+BIN2+BIN3+BIN4
BIN4=None
BIN2=None
BIN1=None
BIN3=None

if time == ALL_TIMES:
 X = X.mean("Time")
 PH = PH.mean(0)
 PHB = PHB.mean(0)
 T = T.mean("Time")
 P = P.mean("Time")

K = (R*T)/P
T=None
P=None
CONC = X/K             ##unit = ug*m^-3
X=None
GPH = (PH+PHB)/9.81
PH=None
PHB=None
###Plot1 Guwahati
sel_lat =26.67129
sel_lon = 89.4322
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
x=X[:,i,j]
y=GPH[:,i,j]
print(x)
print(y)
plt.plot(x,y)
#plt.savefig('guwahati.png',dpi=500)
plt.ylim(0,5000)
plt.show()
X=None
GPH=None


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