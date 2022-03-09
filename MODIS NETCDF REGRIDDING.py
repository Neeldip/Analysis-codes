
###RUN MODIS HDF TO NETCDF.PY FIRST
###use this program to cREATE THE TEMPLATE FOR REGRIDDING
import netCDF4 as nc
import numpy as np
fn = r'F:\AOD_ANALYSIS\template.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')
ds.set_fill_on()

latitude = np.arange(6.5,38.625,0.25)
longitude = np.arange(66.5,100.25,0.25)

lat = ds.createDimension('lat', 129)
lon = ds.createDimension('lon', 135)
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))

random = ds.createVariable('random', 'f4', ('lat','lon'),fill_value=np.nan)

lats.units = 'degrees north'
lons.units = 'degrees east'
lats[:] = latitude
lons[:] = longitude
ds.close()

#RUN cdo remapcon2,template.nc input output

