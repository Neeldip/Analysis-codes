import numpy as np
from goodness_of_fit import rmse,me,d
import pandas as pd
obs = pd.read_excel(r'I:\observation_data\WRF-Chem\interpolated\patna.xlsx', sheet_name="patna")  # READ Obserevd EXCEL FILE
model = pd.read_excel(r'I:\observation_data\WRF-Chem\Model\wrfchem_patna.xlsx', sheet_name="patna")  # READ modek data FILE

Tobs= obs["tK"]
RHobs= obs['relh']

Tmod=model['t']
RHmod=model['rh']

##temperature
Tme=me(Tmod,Tobs)
print(Tme)
Trmse= rmse(Tmod,Tobs)
print(Trmse)
Tioa=d(Tmod,Tobs)
print(Tioa)

#RH
RHme=me(RHmod,RHobs)
print(RHme)
Rhrmse=rmse(RHmod,RHobs)
print(Rhrmse)
RHioa=d(RHmod,RHobs)
print(RHioa)

###WND SPEED###########
##convert the 0 and null values to NaN first in observed values manually
WSobs= obs['sknt']
WSmod= model['ws']
WSobs=np.asarray(WSobs)
nanindex=np.argwhere(~np.isnan(WSobs)) #Find index of non-nan values
WSobs=WSobs[~np.isnan(WSobs)]*0.5144 ##select only the non-nan values and convert to m/s
#WSobs=pd.DataFrame(WSobs)
#print(np.shape(WSobs))
WSnonnan=[]
for i in nanindex:
    ws=WSmod[i]
    WSnonnan.append(ws)
WSmodnonnan= np.asarray(WSnonnan)
WSmodnonnan=np.squeeze(WSmodnonnan)
WSme=me(WSmodnonnan,WSobs)
print(WSme)
WSrmse=rmse(WSmodnonnan,WSobs)
print(WSrmse)
WSioa=d(WSmodnonnan,WSobs)
print(WSioa)

###WND DIE###########
##convert the 0 and null values to NaN first in observed values manually
WDobs= obs['drct']
WDmod= model['dir']
WDobs=np.asarray(WDobs)
nanindex=np.argwhere(~np.isnan(WDobs)) #Find index of non-nan values
WDobs=WDobs[~np.isnan(WDobs)] ##select only the non-nan values
#WDobs=pd.DataFrame(WDobs)
#print(np.shape(WDobs))
WDnonnan=[]
for i in nanindex:
    ws=WDmod[i]
    WDnonnan.append(ws)
WDmodnonnan= np.asarray(WDnonnan)
WDmodnonnan=np.squeeze(WDmodnonnan)
WDme=me(WDmodnonnan,WDobs)
print(WDme)
WDrmse=rmse(WDmodnonnan,WDobs)
print(WDrmse)
WDioa=d(WDmodnonnan,WDobs)
print(WDioa)
