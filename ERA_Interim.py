from xarray import open_mfdataset
import numpy as np
import pandas as pd
file = open_mfdataset("I:\ERA5wind\ERA5-apr-u10_v10.nc4")
#u10 = file.variables["u10"][:,:,:]
#v10 = file.variables["v10"][:,:,:]
u10 = file.variables["var165"][:,:,:]
v10 = file.variables["var166"][:,:,:]
print(u10)
#lats= file.variables["latitude"][:]
#lons= file.variables["longitude"][:]
lats= file.variables["lat"][:]
lons= file.variables["lon"][:]
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
u.to_csv("ERA5U10_apr.csv")
v.to_csv("ERA5V10_apr.csv")