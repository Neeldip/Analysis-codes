
#1.FIRST DO REGRIDDING OF MERRA2 DATA USING MERRA2 REGRIDDING.PY
#2.CDO REMAPBIL,TEMPLATE.NC INFILE OUTFILE

from xarray import open_dataset,open_mfdataset
import numpy as np
import netCDF4 as nc
ds = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_lat= ds.variables["XLAT"][0,:,0]
wrf_lon = ds.variables["XLONG"][0,0,:]
#ds2 = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
ds1 =open_mfdataset([r"I:\Merra 2 M2T1NXAER\regridded\dec\10.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\11.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\12.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\13.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\14.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\15.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\16.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\17.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\18.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\dec\19.nc4",])
'''
T = ds2.variables["TEMPERATURE"][:, :, :, :]
PRESS = ds2.variables["PRESSURE"][:, :, :, :]
T = np.mean(T,axis=0)
PRESS = np.mean(PRESS,axis=0)
bc1 = ds.variables["bc_a01"][:,:,:,:]
bc2 = ds.variables["bc_a02"][:,:,:,:]
bc3 = ds.variables["bc_a03"][:,:,:,:]
bc4 = ds.variables["bc_a04"][:,:,:,:]
bc = bc1+bc2+bc3+bc4
bc1=None
bc2=None
bc3=None
bc4=None
BC = np.mean(bc,axis=0)
bc=None
BC = (BC*PRESS)/(287*T)

oin1 = ds.variables["oin_a01"][:,:,:,:]
oin2 = ds.variables["oin_a02"][:,:,:,:]
oin3 = ds.variables["oin_a03"][:,:,:,:]
oin4 = ds.variables["oin_a04"][:,:,:,:]
oin = oin1+oin2+oin3+oin4
oin1=None
oin2=None
oin3=None
oin4=None
OIN = np.mean(oin,axis=0)
oin=None
OIN = (OIN*PRESS)/(287*T)
OIN = OIN.sum(axis=0)

so41 = ds.variables["so4_a01"][:,:,:,:]
so42 = ds.variables["so4_a02"][:,:,:,:]
so43 = ds.variables["so4_a03"][:,:,:,:]
so44 = ds.variables["so4_a04"][:,:,:,:]
so4 = so41+so42+so43+so44
so41=None
so42=None
so43=None
so44=None
SO4 = np.mean(so4,axis=0)
so4=None
SO4 = (SO4*PRESS)/(287*T)

oc1 = ds.variables["oc_a01"][:,:,:,:]
oc2 = ds.variables["oc_a02"][:,:,:,:]
oc3 = ds.variables["oc_a03"][:,:,:,:]
oc4 = ds.variables["oc_a04"][:,:,:,:]
oc = oc1+oc2+oc3+oc4
oc1=None
oc2=None
oc3=None
oc4=None
OC = np.mean(oc,axis=0)
oc=None
OC = (OC*PRESS)/(287*T)
'''
bcsmass = ds1.variables["BCCMASS"][:,:,:]
bcsmass = np.nanmean(bcsmass,axis=0)*10**9 ##CONVERT FROM KG/M2 TO UG/M2
bcdiff = bcsmass

dusmass = ds1.variables["DUCMASS"][:,:,:]
dusmass = np.nanmean(dusmass,axis=0)*10**9
dudiff = dusmass

ocsmass = ds1.variables["OCCMASS"][:,:,:]
ocsmass = np.nanmean(ocsmass,axis=0)*10**9
ocdiff = ocsmass

so4smass = ds1.variables["SO4CMASS"][:,:,:]
so4smass = np.nanmean(so4smass,axis=0)*10**9
so4diff = so4smass

fn = r'F:\WRF-CHEM ANALYSIS\MERRA2_columnintegrated_aerosol_mass_dec.nc'   #<----------------
ds3 = nc.Dataset(fn, 'w', format='NETCDF4')
ds3.set_fill_on()

lat = ds3.createDimension('lat', 231)
lon = ds3.createDimension('lon', 299)
lats = ds3.createVariable('lat', 'f4', ('lat',))
lons = ds3.createVariable('lon', 'f4', ('lon',))
#BCDIFF = ds3.createVariable('bcdiff', 'f4', ('lat','lon'),fill_value=np.nan)
#BCCHEM = ds3.createVariable('bcchem', 'f4', ('lat','lon'),fill_value=np.nan)
BCMERRA = ds3.createVariable('bcmerra', 'f4', ('lat','lon'),fill_value=np.nan)
#OCDIFF = ds3.createVariable('ocdiff', 'f4', ('lat','lon'),fill_value=np.nan)
#OCCHEM = ds3.createVariable('occhem', 'f4', ('lat','lon'),fill_value=np.nan)
OCMERRA = ds3.createVariable('ocmerra', 'f4', ('lat','lon'),fill_value=np.nan)
#OINDIFF = ds3.createVariable('dustdiff', 'f4', ('lat','lon'),fill_value=np.nan)
#OINCHEM = ds3.createVariable('dustchem', 'f4', ('lat','lon'),fill_value=np.nan)
OINMERRA = ds3.createVariable('dustmerra', 'f4', ('lat','lon'),fill_value=np.nan)
#SO4DIFF = ds3.createVariable('so4diff', 'f4', ('lat','lon'),fill_value=np.nan)
#SO4CHEM = ds3.createVariable('so4chem', 'f4', ('lat','lon'),fill_value=np.nan)
SO4MERRA = ds3.createVariable('so4merra', 'f4', ('lat','lon'),fill_value=np.nan)

lats.units = 'degrees north'
lons.units = 'degrees east'
BCMERRA.units = 'column sum of bc all variables units in ug/m2'
lats[:] = wrf_lat
lons[:] = wrf_lon

for i in range(0,231,1):
    for j in range(0,299,1):
        #BCDIFF[i,j] = bcdiff[i, j]
        #BCCHEM[i, j] = BC[i, j]
        BCMERRA[i, j] = bcsmass[i, j]
        #OCDIFF[i,j]=ocdiff[i,j]
        #OCCHEM[i, j] = OC[i, j]
        OCMERRA[i, j] = ocsmass[i, j]
        #OINDIFF[i,j]=dudiff[i,j]
        #OINCHEM[i, j] = OIN[i, j]
        OINMERRA[i, j] = dusmass[i, j]
        #SO4DIFF[i,j]=so4diff[i,j]
        #SO4CHEM[i, j] = SO4[i, j]
        SO4MERRA[i, j] = so4smass[i, j]
ds3.close()
