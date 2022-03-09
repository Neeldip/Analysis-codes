# 1.FIRST DO REGRIDDING OF MERRA2 DATA USING MERRA2 REGRIDDING.PY
# 2.CDO REMAPBIL,TEMPLATE.NC INFILE OUTFILE
import pandas
from xarray import open_dataset, open_mfdataset
import numpy as np

from goodness_of_fit import rmse, d, me, r_pearson

#wrf_file = open_dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = open_dataset(r"G:\WRF_Chem_Output\202\NOR\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

#wrf_file = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf_file =  open_dataset(r"K:\WRF_Chem_Output\202\April\only_emiss_IGP\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = open_dataset(r"K:\WRF_Chem_Output\202\April\only_emiss_IGP\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
trad = open_dataset(r"G:\WRF_Chem_Output\202\MYNN3_WRF\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
wrf_lats = wrf_file.variables['XLAT'][0, :, 0]
wrf_lons = wrf_file.variables['XLONG'][0, 0, :]

#OBS= pandas.DataFrame(pandas.read_excel(r"I:\observation_data\WRF-Chem\interpolated\bangalore.xlsx",sheet_name="bangalore"))
#OBS= pandas.DataFrame(pandas.read_excel(r"I:\observation_data\WRF-Chem\interpolated\delhi.xlsx",sheet_name="delhi"))
#OBS= pandas.DataFrame(pandas.read_excel(r"I:\observation_data\WRF-Chem\interpolated\guwahati.xlsx",sheet_name="guwahati"))
#OBS= pandas.DataFrame(pandas.read_excel(r"I:\observation_data\WRF-Chem\interpolated\kolkata.xlsx",sheet_name="kolkata"))
#OBS= pandas.DataFrame(pandas.read_excel(r"I:\observation_data\WRF-Chem\interpolated\mumbai.xlsx",sheet_name="mumbai"))
OBS= pandas.DataFrame(pandas.read_excel(r"I:\observation_data\WRF-Chem\interpolated\patna.xlsx",sheet_name="patna"))

OBS= OBS.iloc[::2]  #select alternate rows with same timing as model output
#print(OBS)
Tobs= np.array(OBS.iloc[:,3]) ##select columns
RHobs= np.array(OBS.iloc[:,5])
WSobs= np.array(OBS.iloc[:,7])
WDobs=np.array(OBS.iloc[:,6])
size=np.size(Tobs)

lats = [13.20,     28.56,   26.10,    22.65,    19.10,   25.59]
lons = [77.70,     77.11,   91.58,    88.45,    72.86,   85.08]
#      bangalore,  delhi,   guwahati  kolkata,  mumbai,  patna,

sel_lat=25.59
sel_lon=86.08

a = abs(wrf_lats - sel_lat) + abs(wrf_lons - sel_lon)
i, j = np.unravel_index(a.argmin(), a.shape)


T = trad.variables["T2"][0:size, i, j]  ##sensible temperature,units = K
RH = trad.variables["RH"][0:size,0, i, j]
SPEED = trad.variables["SPEED"][0:size,0, i, j]
DIRECTION = trad.variables["DIR"][0:size,0, i, j]

##############Temperature################
RMSE_T = rmse(T, Tobs)
IOA_T = d(T, Tobs)
ME_T = me(T, Tobs)
CORR_T = r_pearson(T, Tobs)
print("Temperature")
print(RMSE_T)
print(IOA_T)
print(ME_T)
print(CORR_T)

##############RH################
#print(RH)
#print(RHobs)
RMSE_RH = rmse(RH, RHobs)
IOA_RH = d(RH, RHobs)
ME_RH = me(RH, RHobs)
CORR_RH = r_pearson(RH, RHobs)
print("RH")
print(RMSE_RH)
print(IOA_RH)
print(ME_RH)
print(CORR_RH)

###WND SPEED###########
##convert the 0 and null values to NaN first in observed values manually
#WSobs= obs['sknt']
#WSmod= model['ws']
WSobs=np.asarray(WSobs)
nanindex=np.argwhere(~np.isnan(WSobs)) #Find index of non-nan values
WSobs=WSobs[~np.isnan(WSobs)]*0.5144 ##select only the non-nan values and convert to m/s
#WSobs=pd.DataFrame(WSobs)
#print(np.shape(WSobs))
WSnonnan=[]
for i in nanindex:
    ws=SPEED[i]
    WSnonnan.append(ws)
WSmodnonnan= np.asarray(WSnonnan)
WSmodnonnan=np.squeeze(WSmodnonnan)
#print(WSmodnonnan)
#print(WSobs)
print("WS")
WSrmse=rmse(WSmodnonnan,WSobs)
print(WSrmse)
WSioa=d(WSmodnonnan,WSobs)
print(WSioa)
WSme=me(WSmodnonnan,WSobs)
print(WSme)
WScorr=r_pearson(WSmodnonnan,WSobs)
print(WScorr)

###WND DIE###########
##convert the 0 and null values to NaN first in observed values manually
#WDobs= obs['drct']
#WDmod= model['dir']
WDobs=np.asarray(WDobs)
nanindex=np.argwhere(~np.isnan(WDobs)) #Find index of non-nan values
WDobs=WDobs[~np.isnan(WDobs)] ##select only the non-nan values
#WDobs=pd.DataFrame(WDobs)
#print(np.shape(WDobs))
WDnonnan=[]
for i in nanindex:
    ws=DIRECTION[i]
    WDnonnan.append(ws)
WDmodnonnan= np.asarray(WDnonnan)
WDmodnonnan=np.squeeze(WDmodnonnan)

print("WD")
WDrmse=rmse(WDmodnonnan,WDobs)
print(WDrmse)
WDioa=d(WDmodnonnan,WDobs)
print(WDioa)
WDme=me(WDmodnonnan,WDobs)
print(WDme)
WDcorr=r_pearson(WDmodnonnan,WDobs)
print(WDcorr)

