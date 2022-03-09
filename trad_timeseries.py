from xarray import open_mfdataset
import numpy as np
import pandas as pd
import xarray as x
file = open_mfdataset("G:\WRF_Outputs\Monsoon\Hong\wrf_trad_d03.nc")
df = pd.read_excel(r"I:\asos\GHY\ghyjuly2018.xlsx",sheet_name="Sheet5")
#timesteps= df["timestep in python"].to_list()
timesteps= df["hongtimesteps"].to_list()
print(timesteps)

lats= file.variables["XLAT"][0,:] #0-->only 1st step of time dimension
lons= file.variables["XLONG"][0,:]

##With ncks,ERA-Interim nearest gridpoint to 26.1,91.583 is at 26.125,91.625
#so selected here 26.125,91.625

sel_lat =26.1
sel_lon = 91.58
a = abs(lats-sel_lat)+abs(lons-sel_lon)

i,j = np.unravel_index(a.argmin(), a.shape)
u10 = file.variables["U10"][timesteps,i,j] #starting from 1_*_2018_00 for every 3rd step
#print(u10)
v10 = file.variables["V10"][timesteps,i,j]
#print(v10)
#u10 = file.variables["U10"][144:888:3,:,:] #for hong monsoon
#v10 = file.variables["V10"][144:888:3,:,:]

#print(i)#-->latitude
#print(j)#-->longitude
#u = u10[:,i,j]
#v = v10[:,i,j]
u = np.asarray(u10)
v = np.asarray(v10)
#u=pd.DataFrame(u10)
#v=pd.DataFrame(v10)
#print(u)
#print(v)
pd.DataFrame(u).to_csv("U10_jul_hong.csv")
pd.DataFrame(v).to_csv("V10_jul_hong.csv")
#u.to_csv("U10_apr_acm2.csv")
#v.to_csv("V10_apr_acm2.csv")