from xarray import open_dataset
import wrf as w
import numpy as np
import pandas as pd
import xarray as x
from netCDF4 import Dataset
#file =open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
file =open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")

sel_lat =25.47598
sel_lon =94.2562
lats= file.variables["XLAT"][0,:] #0-->only 1st step of time dimension
lons= file.variables["XLONG"][0,:]

a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)

print(i)
print(j)
RAINC= file.variables["RAINC"][:,i,j]
RAINNC= file.variables["RAINNC"][:,i,j]

RAINC=pd.DataFrame(RAINC)
RAINNC=pd.DataFrame(RAINNC)
RAIN=RAINC.append(RAINNC)

RAIN.to_csv('F:\WRF-CHEM ANALYSIS\RAINk.csv')

