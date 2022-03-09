import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
from xarray import open_dataset
#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
from netCDF4 import Dataset
import cartopy.crs as ccrs

wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

#wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
trad1=Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf_lat = wrf_file.variables["XLAT"][0,:,0]
wrf_lon = wrf_file.variables["XLONG"][0,0,:]

#print(wrf_lat)
#print(wrf_lon)
time = ALL_TIMES
#time = 0
BIN1 = getvar(wrf_file, "QCLOUD", timeidx=time)
#BIN2 = getvar(wrf_file, "bc_a02", timeidx=time)
#BIN3 = getvar(wrf_file, "bc_a03", timeidx=time)
#BIN4 = getvar(wrf_file, "bc_a04", timeidx=time)
PH = getvar(wrf_file,"PH",timeidx=time)
PHB = getvar(wrf_file,"PHB",timeidx=time)
T = getvar(trad, "TEMPERATURE", timeidx=time) ##sensible temperature,units = K
P = getvar(trad, "PRESSURE", timeidx=time)  ##units = Pascal
R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 #+ BIN2 + BIN3 + BIN4
BIN1=None
#BIN2=None
#BIN3=None
#BIN4=None
GPH=(PH+PHB)/9.81
PH=None
PHB=None

if time == ALL_TIMES:
    T = T.mean(axis=0)
    P = P.mean(axis=0)
    X = X.mean(axis=0)
    GPH= GPH.mean(axis=0)

#X=X.mean(axis=0)
#T=T.mean(axis=0)
#P=P.mean(axis=0)
#GPH=GPH.mean(axis=0)

CONC= (X*P)/(R*T)
#print(CONC)
X=None
P=None
T=None
#CLWP=CONC.sum(axis=0)
#print(CLWP)

ds= Dataset(r"F:\WRF-CHEM ANALYSIS\CLWP_DIFF_NOR-NOFEED(1-2km).nc",mode='w',format='NETCDF4')

lat = ds.createDimension('lat', 231)
lon = ds.createDimension('lon', 299)
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('CLWP', 'f4', ('lat', 'lon',))
value.units = 'kg/m2'
#value1 = ds.createVariable('column mass sum', 'f4', ('lat', 'lon',),fill_value=np.nan)
#value1.units = 'column_integrated_sum kg/m2'
lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = wrf_lat
lons[:] = wrf_lon

value1=np.empty((231,299))
for i in range(0,231,1):
    for j in range(0,299,1):
        column=[]
        #for k in range(0, 23, 1):
        for k in range(13,18,1):
            #print(k)
            conc = (CONC[k,i,j])*(GPH[k+1,i,j]-GPH[k,i,j])   ###ug/m3 *m = ug/m2
            column.append(conc)
            if len(column)==5:
                #print(k)
                print(i)
                print(j)
                column_sum= np.sum(column)
                value1[i,j]=column_sum
                #print(value)

GPH=None
CONC=None

###2nd file

BIN1 = getvar(wrf_file1, "QCLOUD", timeidx=time)
#BIN2 = getvar(wrf_file, "bc_a02", timeidx=time)
#BIN3 = getvar(wrf_file, "bc_a03", timeidx=time)
#BIN4 = getvar(wrf_file, "bc_a04", timeidx=time)
PH = getvar(wrf_file1,"PH",timeidx=time)
PHB = getvar(wrf_file1,"PHB",timeidx=time)
T = getvar(trad1, "TEMPERATURE", timeidx=time) ##sensible temperature,units = K
P = getvar(trad1, "PRESSURE", timeidx=time)  ##units = Pascal
R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 #+ BIN2 + BIN3 + BIN4
BIN1=None
#BIN2=None
#BIN3=None
#BIN4=None
GPH1=(PH+PHB)/9.81
PH=None
PHB=None

if time == ALL_TIMES:
    T = T.mean(axis=0)
    P = P.mean(axis=0)
    X = X.mean(axis=0)
    GPH= GPH1.mean(axis=0)

CONC= (X*P)/(R*T)
#print(CONC)
X=None
P=None
T=None

value2=np.empty((231,299))
for i in range(0,231,1):
    for j in range(0,299,1):
        column=[]
        #for k in range(0, 23, 1):
        for k in range(13,18,1):
            #print(k)
            conc = (CONC[k,i,j])*(GPH[k+1,i,j]-GPH[k,i,j])   ###Kg/m3 *m = Kg/m2
            column.append(conc)
            if len(column)==5:
                #print(k)
                print(i)
                print(j)
                column_sum= np.sum(column)
                value2[i,j]=column_sum
                #print(value)

for i in range(0,231,1):
    for j in range(0,299,1):
        value[i,j]=value1[i,j]-value2[i,j]

precp1 = value[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = value[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = value[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = value[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = value[168:178, 174:216]
print(np.nanmean(precp1))

ds.close()
GPH=None
CONC=None
