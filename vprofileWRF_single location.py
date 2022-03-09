import numpy as np
from numpy import *
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES,destagger)
import matplotlib.patches as mpatch
##ACM2
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

##MYNN3
#wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset("G:\WRF_Chem_Output\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffOFF\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"H:\GF\wrfout_d01_2018-04-10_00%3A00%3A00")
##MYNN3
#wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset("G:\WRF_Chem_Output\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffOFF\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrf_trad_d03.nc")
#wrf_file1 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrf_trad_d03.nc")
#wrf_file2 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYJ\wrf_trad_d03.nc")
#wrf_file3 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrf_trad_d03.nc")
#wrf_file4 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_YSU\wrf_trad_d03.nc")
#wrf_file5 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_QNSE\wrf_trad_d03.nc")
#ds = pd.read_excel(r"I:\observation_data\APRIL\DATA\Radiosonde_stations_NE\APRIL\WYOMING\New.xlsx",sheet_name="Sheet2")
#lats = wrf_file.variables['XLAT'][0,:,:] ##231*299
#lons = wrf_file.variables['XLONG'][0,:,:] ##231*299


wrf_file = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
wrf_file1 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrfout_d03.nc")
wrf_file2 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYJ\wrfout_d03.nc")
wrf_file3 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrfout_d03.nc")
wrf_file4 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_QNSE\wrfout_d03.nc")
wrf_file5 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_YSU\wrfout_d03.nc")
lats = wrf_file.variables['XLAT'][0,:,:] ##231*299
lons = wrf_file.variables['XLONG'][0,:,:] ##231*299
print(lons)
##SELECT TIME AND LOCATION####
time = 0


if time == ALL_TIMES:
 stagger_dim = 1 ###use stagger_dim=0 when running for single time; use 1 for ALL_TIMES
else:
    stagger_dim = 0

sel_lat =26.1
sel_lon = 91.74
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)


T =  wrf_file.variables['T'][:,:,i,j]  ##if time = ALL_TIMES then [:,:,156,136]
T1 = wrf_file1.variables['T'][:,:,i,j]
T2 = wrf_file2.variables['T'][:,:,i,j]
T3 = wrf_file3.variables['T'][:,:,i,j]
T4 = wrf_file4.variables['T'][:,:,i,j]
T5 = wrf_file5.variables['T'][:,:,i,j]

PH = wrf_file.variables['PH'][:,:,i,j]
PH1 = wrf_file1.variables['PH'][:,:,i,j]
PH2 = wrf_file2.variables['PH'][:,:,i,j]
PH3 = wrf_file3.variables['PH'][:,:,i,j]
PH4 = wrf_file4.variables['PH'][:,:,i,j]
PH5 = wrf_file5.variables['PH'][:,:,i,j]

PHB = wrf_file.variables['PHB'][:,:,i,j]
PHB1 = wrf_file1.variables['PHB'][:,:,i,j]
PHB2 = wrf_file2.variables['PHB'][:,:,i,j]
PHB3 = wrf_file3.variables['PHB'][:,:,i,j]
PHB4 = wrf_file4.variables['PHB'][:,:,i,j]
PHB5 = wrf_file5.variables['PHB'][:,:,i,j]

PH = destagger(PH, stagger_dim=stagger_dim, meta=False)
PH1 = destagger(PH1, stagger_dim=stagger_dim, meta=False)
PH2 = destagger(PH2, stagger_dim=stagger_dim, meta=False)
PH3 = destagger(PH3, stagger_dim=stagger_dim, meta=False)
PH4 = destagger(PH4, stagger_dim=stagger_dim, meta=False)
PH5 = destagger(PH5, stagger_dim=stagger_dim, meta=False)
PHB = destagger(PHB, stagger_dim=stagger_dim, meta=False)
PHB1 = destagger(PHB1, stagger_dim=stagger_dim, meta=False)
PHB2 = destagger(PHB2, stagger_dim=stagger_dim, meta=False)
PHB3 = destagger(PHB3, stagger_dim=stagger_dim, meta=False)
PHB4 = destagger(PHB4, stagger_dim=stagger_dim, meta=False)
PHB5 = destagger(PHB5, stagger_dim=stagger_dim, meta=False)

if time == ALL_TIMES:
 T = T.mean(axis=0)
 T1 = T1.mean(axis=0)
 T2 = T2.mean(axis=0)
 T3 = T3.mean(axis=0)
 T4 = T4.mean(axis=0)
 T5 = T5.mean(axis=0)
 PH = PH.mean(axis=0)
 PH1= PH1.mean(axis=0)
 PH2 = PH2.mean(axis=0)
 PH3 = PH3.mean(axis=0)
 PH4 = PH4.mean(axis=0)
 PH5 = PH5.mean(axis=0)
 PHB = PHB.mean(axis=0)
 PHB1 = PHB1.mean(axis=0)
 PHB2 = PHB2.mean(axis=0)
 PHB3 = PHB3.mean(axis=0)
 PHB4 = PHB4.mean(axis=0)
 PHB5 = PHB5.mean(axis=0)

T = T+300
T1 = T1+300
T2 = T2+300
T3 = T3+300
T4 = T4+300
T5 = T5+300
GPH = (PH+PHB)/9.81
GPH1 = (PH1+PHB1)/9.81
GPH2 = (PH2+PHB2)/9.81
GPH3 = (PH3+PHB3)/9.81
GPH4 = (PH4+PHB4)/9.81
GPH5 = (PH5+PHB5)/9.81



t=T
y=GPH
t1=T1
y1=GPH1
t2=T2
y2=GPH2
t3=T3
y3=GPH3
t4=T4
y4=GPH4
t5=T5
y5=GPH5
plt.figure(figsize=(7,10))
plt.plot(t,y,color='dodgerblue')
plt.plot(t1,y1,color='coral')
plt.plot(t2,y2,color='gold')
plt.plot(t3,y3,color='lime')
plt.plot(t4,y4,color='aqua')
plt.plot(t5,y5,color='grey')
#plt.plot(ds["rh"],ds["h"],color="black",linestyle="--")
blue = mpatch.Patch(color="dodgerblue", label="ACM2")
coral = mpatch.Patch(color="coral", label="HONG")
gold = mpatch.Patch(color="gold", label="MYJ")
chocolate = mpatch.Patch(color="lime", label="MYNN3")
aqua = mpatch.Patch(color="aqua", label="QNSE")
grey = mpatch.Patch(color="grey", label="YSU")
#black = mpatch.Patch(color="black", label="Radiosonde")
#plt.legend(handles=[blue,coral,gold,chocolate,aqua,grey,black])
plt.legend(handles=[blue,coral,gold,chocolate,aqua,grey])
plt.ylim(0,5000)
plt.xlim(300,320)
plt.xlabel("Theta (K)")
plt.ylabel("Height (m)")
#plt.savefig('T_profile_ghy.png',dpi=300,bbox_inches='tight')
plt.show()


'''
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