import netCDF4 as nc
import numpy as np
from scipy.stats import spearmanr
from xarray import open_mfdataset


fn = 'F:\AOD_ANALYSIS\AOD_CER_apr_correlation.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')
import pandas as pd
ds.set_fill_on()

lat = ds.createDimension('lat', 129)
lon = ds.createDimension('lon', 135)

lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('Correlation AOD and RE', 'f4', ('lat', 'lon',),fill_value=np.nan)
value.units = 'Spearman correlation between AOD and cloud droplet radius'
pvalue = ds.createVariable('pvalue AOD RE', 'f4', ('lat', 'lon',),fill_value=np.nan)
pvalue.units = 'Spearman correlation p-value AOD and RE'

value1 = ds.createVariable('Correlation AOD and COD', 'f4', ('lat', 'lon',),fill_value=np.nan)
value1.units = 'Spearman correlation between AOD and COD'
pvalue1 = ds.createVariable('pvalue AOD and COD', 'f4', ('lat', 'lon',),fill_value=np.nan)
pvalue1.units = 'Spearman correlation p-value AOD and COD'
lats.units = 'degrees north'
lons.units = 'degrees east'

lats[:] = np.arange(6.5, 38.75, 0.25)
lons[:] = np.arange(66.5, 100.25, 0.25)

for i in range(0,129,1):
    for j in range(0,135,1):
        print(i)
        print(j)
        DS = open_mfdataset(r"F:\AOD_ANALYSIS\APR_regridded.nc")
        AOD = DS.variables["AOD"][:,i,j]
        CER = DS.variables["CER"][:, i, j]
        COD = DS.variables["COD"][:, i, j]
        if np.count_nonzero(np.isnan(AOD)) > 10 or np.count_nonzero(np.isnan(CER)) > 10 or np.count_nonzero(np.isnan(COD)) > 10  :  # count # of nan values
           value[i, j] = np.nan
           pvalue[i, j] = np.nan
        else:
           spearman_correl,p = spearmanr(AOD,CER,nan_policy='omit',axis=0)
           spearman_correl1, p1 = spearmanr(AOD, COD, nan_policy='omit', axis=0)
           #print(spearman_correl)
           value[i,j]= spearman_correl
           pvalue[i,j]=p
           value1[i,j]= spearman_correl1
           pvalue1[i,j]=p1
ds.close()


