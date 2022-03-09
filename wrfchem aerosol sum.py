import numpy as np
from  xarray import open_dataset
import netCDF4 as nc
from netCDF4 import Dataset
from wrf import getvar,ALL_TIMES
wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
LAT = wrf_file.variables["XLAT"][0,:,0]
LON = wrf_file.variables["XLONG"][0,0,:]

fn = 'F:\WRF-CHEM ANALYSIS\WRF_CHEM_BC_April.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')
ds.set_fill_on()

lat = ds.createDimension('lat', 231)
lon = ds.createDimension('lon', 299)
lev = ds.createDimension('lev', 44)
time = ds.createDimension('time', 242)

lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
level = ds.createVariable('lev','f4', ('lev',))
times = ds.createVariable('time','f4', ('time',))
value = ds.createVariable('BC', 'f4', ('time','lev','lat', 'lon',),fill_value=np.nan)
value.units = 'BC conc ug/kg dry air'

times[:]= np.arange(0,242,1)
level[:] = np.arange(0,44,1)
lons[:] = LON
lats[:] = LAT
BIN1 = getvar(wrf_file,"bc_a01",timeidx=ALL_TIMES)
BIN2 = getvar(wrf_file,"bc_a02",timeidx=ALL_TIMES)
BIN3 = getvar(wrf_file,"bc_a03",timeidx=ALL_TIMES)
BIN4 = getvar(wrf_file,"bc_a04",timeidx=ALL_TIMES)
BC = BIN1+BIN2+BIN3+BIN4
BIN1=None
BIN2=None
BIN3=None
BIN4=None
for t in range(0, 242, 1):
    for k in range(0,44,1):
        for i in range(0,231,1):
            for j in range(0,299,1):
                    value[t,k,i,j]=BC[t,k,i,j]
                    print(j)

ds.close()
BC=None