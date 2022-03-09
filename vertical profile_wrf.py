import numpy as np
from numpy import *
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES,destagger)

##ACM2
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
##MYNN3
#wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset("G:\WRF_Chem_Output\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffOFF\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"H:\GF\wrfout_d01_2018-04-10_00%3A00%3A00")
lats = wrf_file.variables['XLAT'][0,:,:] ##231*299
lons = wrf_file.variables['XLONG'][0,:,:] ##231*299


#wrf_file1 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrf_trad_d03.nc")
#wrf_file2 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYJ\wrf_trad_d03.nc")
#wrf_file3 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrf_trad_d03.nc")
#wrf_file4 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_YSU\wrf_trad_d03.nc")
#wrf_file5 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_QNSE\wrf_trad_d03.nc")

##SELECT TIME AND LOCATION####
sel_lat =26.67129#24.5#26.5#28.0
sel_lon =89.4322#93.0# 92.2#96.0
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)
print(i)
print(j)

i=176
j=183
time =60
if time == ALL_TIMES:
 stagger_dim = 1 ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
else:
    stagger_dim = 0


T = getvar(wrf_file, "QCLOUD", timeidx=time)
T1 = getvar(wrf_file1, "QCLOUD", timeidx=time)
PH = getvar(wrf_file,"PH",timeidx=time)
PHB = getvar(wrf_file,"PHB",timeidx=time)
PH = destagger(PH, stagger_dim=stagger_dim, meta=False)   ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)
PH1 = getvar(wrf_file1,"PH",timeidx=time)
PHB1 = getvar(wrf_file1,"PHB",timeidx=time)
PH1 = destagger(PH1, stagger_dim=stagger_dim, meta=False)   ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
PHB1 = destagger(PHB1, stagger_dim=stagger_dim, meta=False)

if time == ALL_TIMES:
 T = T.mean("Time")
 PH = PH.mean(0)
 PHB = PHB.mean(0)
 T1 = T1.mean("Time")
 PH1 = PH1.mean(0)
 PHB1 = PHB1.mean(0)

T = T
T1 = T1
GPH = (PH+PHB)/9.81
GPH1 = (PH1+PHB1)/9.81

###Plot1

t=T[:,i,j]
y=GPH[:,i,j]
t1=T1[:,i,j]
y1=GPH1[:,i,j]
plt.plot(t,y)
plt.plot(t1,y1)
plt.ylim(0,16000)
#plt.xlim(300,320)
plt.title("blue=NORM, yellow=nofeedback")
#plt.savefig(r"F:\WRF-CHEM ANALYSIS\rainfall\LOCATION ANALYSIS\26.67129 89.4322\THM_time=61.png",dpi=600,bbox_inches='tight')
plt.show()
T=None
T1=None
PH=None
PH1=None
PHB=None
PHB1=None
GPH=None
GPH1=None
###
'''
sel_lat =23.5
sel_lon = 88.75
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
t=T[:,i,j]
y=GPH[:,i,j]
t1=T1[:,i,j]
y1=GPH1[:,i,j]
plt.plot(t,y)
plt.plot(t1,y1)
plt.ylim(0,5000)
plt.xlim(290,330)
#plt.savefig('Potential_Temperature_profile_lat='+str(sel_lat)+'_lon='+str(sel_lon)+'.png',dpi=500)
plt.show()



###
sel_lat =24.85
sel_lon = 93.25
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
t=T[:,i,j]
y=GPH[:,i,j]
t1=T1[:,i,j]
y1=GPH1[:,i,j]
plt.plot(t,y)
plt.plot(t1,y1)
plt.ylim(0,5000)
plt.xlim(290,330)
#plt.savefig('Potential_Temperature_profile_lat='+str(sel_lat)+'_lon='+str(sel_lon)+'.png',dpi=500)
plt.show()




###Plot4 Bhubaneswar
sel_lat =26.2
sel_lon =89.13
a = abs(lats-sel_lat)+abs(lons-sel_lon)
#print(shape(a)) = 231*299
i,j = np.unravel_index(a.argmin(), a.shape)
t=T[:,i,j]
y=GPH[:,i,j]
t1=T1[:,i,j]
y1=GPH1[:,i,j]
plt.plot(t,y)
plt.plot(t1,y1)
plt.ylim(0,5000)
plt.xlim(290,330)
#plt.savefig('Potential_Temperature_profile_lat='+str(sel_lat)+'_lon='+str(sel_lon)+'.png',dpi=500)
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