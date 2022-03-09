import numpy as np
from numpy import *
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES,destagger)
import pandas as pd
import matplotlib.patches as mpatch
##ACM2

#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

##MYNN3
#wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffON\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset("G:\WRF_Chem_Output\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\bc_directeffOFF\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrf_trad_d03.nc")
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
wrf_file1= Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrf_trad_d03.nc")
#wrf_file2 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYJ\wrf_trad_d03.nc")
#wrf_file3 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrf_trad_d03.nc")
#wrf_file4 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_YSU\wrf_trad_d03.nc")
#wrf_file5 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_QNSE\wrf_trad_d03.nc")
#ds = pd.read_excel(r"I:\observation_data\APRIL\DATA\Radiosonde_stations_NE\APRIL\WYOMING\New.xlsx",sheet_name="Sheet2")
lats = wrf_file.variables['XLAT'][0,:,:] ##231*299
lons = wrf_file.variables['XLONG'][0,:,:] ##231*299

##SELECT TIME AND LOCATION####
time =0 ##if time != ALL_TIMES, set any other value. eg, time = 0
sel_lat =26.67129
sel_lon = 89.4322
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)

print(i)
print(j)
#T1 = getvar(wrf_file1, "RH", timeidx=time)
#T2 = getvar(wrf_file2, "RH", timeidx=time)
#T3 = getvar(wrf_file3, "RH", timeidx=time)
#T4 = getvar(wrf_file4, "RH", timeidx=time)
#T5 = getvar(wrf_file5, "RH", timeidx=time)
T =  wrf_file.variables['TEMPERATURE'][time,:,i,j]  ##if time = ALL_TIMES then [:,:,156,136]
T1 = wrf_file1.variables['TEMPERATURE'][time,:,i,j]
#T2 = wrf_file2.variables['RH'][24,:,i,j]
#T3 = wrf_file3.variables['RH'][24,:,i,j]
#T4 = wrf_file4.variables['RH'][24,:,i,j]
#T5 = wrf_file5.variables['RH'][24,:,i,j]

#T=T[:,:,156,136]
#T1=T1[:,:,156,136]
#T2=T2[:,:,156,136]
#T3=T3[:,:,156,136]
#T4=T4[:,:,156,136]
#T5=T5[:,:,156,136]


#PH = getvar(wrf_file,"GEOHEIGHT",timeidx=0)
#PH1 = getvar(wrf_file1,"GEOHEIGHT",timeidx=0)
#PH2 = getvar(wrf_file2,"GEOHEIGHT",timeidx=0)
#PH3 = getvar(wrf_file3,"GEOHEIGHT",timeidx=0)
#PH4 = getvar(wrf_file4,"GEOHEIGHT",timeidx=0)
#PH5 = getvar(wrf_file5,"GEOHEIGHT",timeidx=0)
PH = wrf_file.variables['GEOHEIGHT'][time,:,i,j]
PH1 = wrf_file1.variables['GEOHEIGHT'][time,:,i,j]
#PH2 = wrf_file2.variables['GEOHEIGHT'][24,:,i,j]
#PH3 = wrf_file3.variables['GEOHEIGHT'][24,:,i,j]
#PH4 = wrf_file4.variables['GEOHEIGHT'][24,:,i,j]
#PH5 = wrf_file5.variables['GEOHEIGHT'][24,:,i,j]
if time == ALL_TIMES:
 #T = T.mean("Time")
 #T1 = T1.mean("Time")
 #T2 = T2.mean("Time")
 #T3 = T3.mean("Time")
 #T4 = T4.mean("Time")
 #T5 = T5.mean("Time")
 T = T.mean(axis=0)
 T1 = T1.mean(axis=0)
 #T2 = T2.mean(axis=0)
 #T3 = T3.mean(axis=0)
 #T4 = T4.mean(axis=0)
 #T5 = T5.mean(axis=0)
 #PH = PH.mean(axis=0)
 #PH1 = PH1.mean(axis=0)
 #PH2 = PH2.mean(axis=0)
 #PH3 = PH3.mean(axis=0)
 #PH4 = PH4.mean(axis=0)
 #PH5 = PH5.mean(axis=0)
T = T
T1 = T1
#T2 = T2
#T3 = T3
#T4 = T4
#T5 = T5
GPH = PH
GPH1 = PH1
#GPH2 = PH2
#GPH3 = PH3
#GPH4 = PH4
#GPH5 = PH5

###Plot1

#t=T[:,i,j]
#y=GPH[:,i,j]
#t1=T1[:,i,j]
#y1=GPH1[:,i,j]
#t2=T2[:,i,j]
#y2=GPH2[:,i,j]
#t3=T3[:,i,j]
#y3=GPH3[:,i,j]
#t4=T4[:,i,j]
#y4=GPH4[:,i,j]
#t5=T5[:,i,j]
#y5=GPH5[:,i,j]
t=T
y=GPH
t1=T1
y1=GPH1
#t2=T2
#y2=GPH2
#t3=T3
#y3=GPH3
#t4=T4
#y4=GPH4
#t5=T5
#y5=GPH5
#print(y)
#print(y1)
#print(y2)
#print(y3)
#print(y4)
#print(y5)
plt.figure(figsize=(5.5,10))
plt.plot(t,y,color='dodgerblue')
plt.plot(t1,y1,color='coral')
#plt.plot(t2,y2,color='gold')
#plt.plot(t3,y3,color='lime')
#plt.plot(t4,y4,color='aqua')
#plt.plot(t5,y5,color='grey')
#plt.plot(ds["rh"],ds["h"],color="black",linestyle="--")
#blue = mpatch.Patch(color="dodgerblue", label="ACM2")
#coral = mpatch.Patch(color="coral", label="HONG")
#gold = mpatch.Patch(color="gold", label="MYJ")
#chocolate = mpatch.Patch(color="lime", label="MYNN3")
#aqua = mpatch.Patch(color="aqua", label="QNSE")
#grey = mpatch.Patch(color="grey", label="YSU")
#black = mpatch.Patch(color="black", label="Radiosonde")
#plt.legend(handles=[blue,coral,gold,chocolate,aqua,grey,black])
#plt.legend(handles=[blue,coral,gold,chocolate,aqua,grey])
plt.ylim(0,16000)
plt.xlim(270,305)
#plt.xlabel("RH (%)")
plt.ylabel("Height (m)")
plt.title("blue=NORM, yellow=nofeedback")
#plt.savefig(r"F:\WRF-CHEM ANALYSIS\rainfall\LOCATION ANALYSIS\26.67129 89.4322\TEMP_time=60.png",dpi=600,bbox_inches='tight')
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