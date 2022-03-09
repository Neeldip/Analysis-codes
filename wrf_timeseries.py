from xarray import open_mfdataset
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as x
from wrf import getvar,ll_to_xy,ALL_TIMES
file = Dataset("G:\WRF_Outputs\Monsoon\ACM2\wrfout_d03.nc")
df = pd.read_excel(r"I:\asos\GHY\ghyjuly2018.xlsx",sheet_name="Sheet5")
timestep= df["timestep in python"].to_list()
#timestep= df["hongtimesteps"].to_list() #for hong monsoon only
print(timestep)

lats= file.variables["XLAT"][0,:] #0-->only 1st step of time dimension
lons= file.variables["XLONG"][0,:]

##With ncks,ERA-Interim nearest gridpoint to 26.1,91.583 is at 26.125,91.625
#so selected here 26.125,91.625

sel_lat =26.1
sel_lon = 91.58
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)

#c=ll_to_xy(file,26.1,91.583,0)
#print(c)

u10 = getvar(file,"ua",timeidx=ALL_TIMES)
v10 = getvar(file,"va",timeidx=ALL_TIMES)
u10 = u10[timestep,0,i,j]
v10 = v10[timestep,0,i,j]
#u10 = file.variables["U10"][144:888:3,:,:] #for hong monsoon
#v10 = file.variables["V10"][144:888:3,:,:]


#print(u10)#-->latitude
#print(v10)#-->longitude
#u = u10[0,i,j]
#v = v10[0,i,j]
#u = np.asarray(u10)
#v = np.asarray(v10)
u=pd.DataFrame(u10)
v=pd.DataFrame(v10)
print(u)
print(v)
pd.DataFrame(u).to_csv("U_lev0_jul_ACM2.csv")
pd.DataFrame(v).to_csv("V_lev0_jul_ACM2.csv")
#u.to_csv("U10_apr_acm2.csv")
#v.to_csv("V10_apr_acm2.csv")

