import netCDF4 as nc
import numpy as np

fn = 'F:\IMD\z_modified_statistics_mar.nc'
ds = nc.Dataset(fn, 'w', format='NETCDF4')
import pandas as pd
ds.set_fill_on()

lat = ds.createDimension('lat', 129)
lon = ds.createDimension('lon', 135)

lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('value', 'f4', ('lat', 'lon',),fill_value=np.nan)
value.units = 'Unknown'
lats.units = 'degrees north'
lons.units = 'degrees east'
lats[:] = np.arange(6.5, 38.75, 0.25)
lons[:] = np.arange(66.5, 100.25, 0.25)

df = pd.read_csv(r"F:\IMD\modified mann kendall\2001-2020\z_modified_statistics_mar.txt",header=None)
df = np.asarray(df)
#print(np.nanmax(df))
#df = pd.DataFrame(df)
#print(df)
df = np.reshape(df,(129,135))

for i in range(0,129,1):
    for j in range(0,135,1):
        if np.isnan(df[i,j]).any() == True:
           value[i,j] = np.nan

        else:
           value[i,j] = df[i,j]
ds.close()