from xarray import open_mfdataset
import numpy as np
import pandas as pd
file = open_mfdataset(r"I:\observation_data\APRIL\DATA\Wind_data\4.ERA5-2018APR_U_V_3hr.nc")
#u10 = file.variables["u10"][:,:,:]
#v10 = file.variables["v10"][:,:,:]
u10 = file.variables["u"][:,88,:,:]
v10 = file.variables["v"][:,88,:,:]
print(u10)
#lats= file.variables["latitude"][:]
#lons= file.variables["longitude"][:]
lats= file.variables["latitude"][:]
lats=lats[::-1]
lons= file.variables["longitude"][:]
print(lats)
print(lons)
sel_lat =26.1
sel_lon = 91.583
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)

#print(i)-->latitude
#print(j)-->longitude
u = u10[:,i,j]
v = v10[:,i,j]
#u = np.asarray(u)
#v = np.asarray(v)
u=pd.DataFrame(u)
v=pd.DataFrame(v)
print(u)
print(v)
u.to_csv("ERA5_Model_U10_apr.csv")
v.to_csv("ERA5_Model_V10_apr.csv")

