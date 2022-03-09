from xarray import open_dataset
ds1= open_dataset(r"H:\WRF_Chem_Output\201\April\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_lat= ds1.variables["XLAT"][0,:,0]
wrf_lon = ds1.variables["XLONG"][0,0,:]
print(wrf_lat)
###cREATE THE TEMPLATE FOR REGRIDDING WITH CDO
import netCDF4 as nc
import numpy as np
fn = r'I:\Merra 2 M2T1NXAER\template.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')
ds.set_fill_on()

#latitude = np.arange(6.5,38.625,0.25)
#longitude = np.arange(66.5,100.25,0.25)

lat = ds.createDimension('lat', 231)
lon = ds.createDimension('lon', 299)
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))

random = ds.createVariable('random', 'f4', ('lat','lon'),fill_value=np.nan)

lats.units = 'degrees north'
lons.units = 'degrees east'
lats[:] = wrf_lat
lons[:] = wrf_lon
ds.close()

#RUN cdo remapcon2,template.nc input output

