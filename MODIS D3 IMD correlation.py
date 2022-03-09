import netCDF4 as nc
import numpy as np
from scipy.stats import spearmanr
from xarray import open_mfdataset


fn = 'F:\AOD_ANALYSIS\REGRIDDING MYD_L3\daily\IMD_MODIS_D3_apr_correlation.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')
import pandas as pd
ds.set_fill_on()

lat = ds.createDimension('lat', 129)
lon = ds.createDimension('lon', 135)

lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('Correlation', 'f4', ('lat', 'lon',),fill_value=np.nan)
value.units = 'Spearman correlation between IMD rainfall and MODIS MYD L3'
pvalue = ds.createVariable('pvalue', 'f4', ('lat', 'lon',),fill_value=np.nan)
pvalue.units = 'Spearman correlation p-value'
lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = np.arange(6.5, 38.75, 0.25)
lons[:] = np.arange(66.5, 100.25, 0.25)

for i in range(0,129,1):
    for j in range(0,135,1):
        imd = open_mfdataset(r"I:\observation_data\APRIL\DATA\Rainfall\IMD_gridded_datasets\RAW data\NETCDF\New folder\Working\apr\2018-Apr.nc")
        rain = imd.variables["RAINFALL"][:, i, j]
        if np.isnan(rain).any():
            value[i, j] = np.nan
            continue
        else:
           AOD = open_mfdataset(r"F:\AOD_ANALYSIS\REGRIDDING MYD_L3\daily\APR_2018_D3_regridded.nc")
           AOD = AOD.variables["AOD"][:,i,j]
           if np.count_nonzero(np.isnan(AOD)) > 10:  # count # of nan values
               value[i, j] = np.nan
           else:
               spearman_correl,p = spearmanr(AOD,rain,nan_policy='omit')
               #print(spearman_correl)
               value[i,j]= spearman_correl
               pvalue[i,j]=p
ds.close()


