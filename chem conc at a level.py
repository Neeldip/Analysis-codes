
#1.FIRST DO REGRIDDING OF MERRA2 DATA USING MERRA2 REGRIDDING.PY
#2.CDO REMAPBIL,TEMPLATE.NC INFILE OUTFILE

from xarray import open_dataset,open_mfdataset
import numpy as np
import netCDF4 as nc
ds = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_lat= ds.variables["XLAT"][0,:,0]
wrf_lon = ds.variables["XLONG"][0,0,:]
ds2 = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#ds1 =open_mfdataset([r"I:\Merra 2 M2T1NXAER\regridded\10.nc4",
#                     r"I:\Merra 2 M2T1NXAER\regridded\11.nc4",
#                     r"I:\Merra 2 M2T1NXAER\regridded\12.nc4",
#                     r"I:\Merra 2 M2T1NXAER\regridded\13.nc4",
#                    r"I:\Merra 2 M2T1NXAER\regridded\14.nc4",
#                    r"I:\Merra 2 M2T1NXAER\regridded\15.nc4",
#                     r"I:\Merra 2 M2T1NXAER\regridded\16.nc4",
#                     r"I:\Merra 2 M2T1NXAER\regridded\17.nc4",
#                     r"I:\Merra 2 M2T1NXAER\regridded\18.nc4",
#                     r"I:\Merra 2 M2T1NXAER\regridded\19.nc4",])

#T = ds2.variables["TEMPERATURE"][:, 0, :, :]
#PRESS = ds2.variables["PRESSURE"][:, 0, :, :]
#T = np.mean(T,axis=0)
#PRESS = np.mean(PRESS,axis=0)
bc1 = ds.variables["bc_a01"][:,0,:,:]
bc1= np.mean(bc1,axis=0)
bc2 = ds.variables["bc_a02"][:,0,:,:]
bc2=np.mean(bc2,axis=0)
bc3 = ds.variables["bc_a03"][:,0,:,:]
bc3=np.mean(bc3,axis=0)
bc4 = ds.variables["bc_a04"][:,0,:,:]
bc4=np.mean(bc4,axis=0)
bc = bc1+bc2+bc3+bc4

num1=(bc1/bc)*100
num2=(bc2/bc)*100
num3=(bc3/bc)*100
num4=(bc4/bc)*100


fn = r'F:\WRF-CHEM ANALYSIS\num_conc.nc'   #<----------------
ds3 = nc.Dataset(fn, 'w', format='NETCDF4')
ds3.set_fill_on()

lat = ds3.createDimension('lat', 231)
lon = ds3.createDimension('lon', 299)
lats = ds3.createVariable('lat', 'f4', ('lat',))
lons = ds3.createVariable('lon', 'f4', ('lon',))
NUM1 = ds3.createVariable('num_a01', 'f4', ('lat','lon'),fill_value=np.nan)
NUM2 = ds3.createVariable('num_a02', 'f4', ('lat','lon'),fill_value=np.nan)
NUM3 = ds3.createVariable('num_a03', 'f4', ('lat','lon'),fill_value=np.nan)
NUM4 = ds3.createVariable('num_a04', 'f4', ('lat','lon'),fill_value=np.nan)


lats.units = 'degrees north'
lons.units = 'degrees east'
NUM1.units = 'all variables units in #'
lats[:] = wrf_lat
lons[:] = wrf_lon

for i in range(0,231,1):
    for j in range(0,299,1):
        NUM1[i,j] = num1[i, j]
        NUM2[i, j] = num2[i, j]
        NUM3[i, j] = num3[i, j]
        NUM4[i,j]=  num4[i,j]
ds3.close()
