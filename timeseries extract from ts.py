import numpy as np
import pandas as pd

ts = pd.read_excel(r'I:\observation_data\WRF-Chem\d01.xlsx', sheet_name="kolkata")  # READ EXCEL FILE


Temp=ts['t']
RH = ts['rh']
spd= ts["speed"]
dir=ts["dir"]

T2=[]
Rh=[]
Spd=[]
Dir=[]
for i in range(0,21600,45):
    T= Temp[i]
    T2.append(T)
    rh=RH[i]
    Rh.append(rh)
    S=spd[i]
    Spd.append(S)
    DIR=dir[i]
    Dir.append(DIR)
#T2=np.reshape(T2,(480,))
#Rh=np.reshape(Rh,(480,))
#Spd=np.reshape(Spd,(480,))
#Dir=np.reshape(Dir,(480,))
T2=pd.DataFrame([T2,Rh,Spd,Dir])

#print(ts)
with pd.ExcelWriter(r'I:\observation_data\WRF-Chem\wrfchem_kolkata.xlsx',
                    ) as writer:
    T2.to_excel(writer, sheet_name="kolkata")  # WRITE TO EXCEL
