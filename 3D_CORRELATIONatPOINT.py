import numpy as np
from numpy import transpose

from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES)
from scipy.stats import pearsonr,kendalltau
#import sys
#np.set_printoptions(threshold=sys.maxsize)
#np.set_printoptions(threshold=np.inf)


wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"G:\WRF_Chem_Output\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

#trad1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")

#####SELECT A 3D VARIABLE
VAR1 = getvar(wrf_file, "QCLOUD", timeidx=ALL_TIMES)
VAR2 = getvar(trad, "TEMPERATURE", timeidx=ALL_TIMES)


###SELECT LOCATION AND LEVEL OF VARIABLE####

lats = wrf_file.variables['XLAT'][0,:,:] ##231*299
lons = wrf_file.variables['XLONG'][0,:,:] ##231*299

sel_lat =24.83
sel_lon = 93.25
LEVEL= 0
a = abs(lats-sel_lat)+abs(lons-sel_lon)
i,j = np.unravel_index(a.argmin(), a.shape)
VAR1 = VAR1[:,LEVEL,1,1]*1000
VAR2 = VAR2[:,LEVEL,1,1]
corr_pearson, _ = pearsonr(VAR1, VAR2)
print(corr_pearson)
corr_kendall, _ = kendalltau(VAR1, VAR2)
print(corr_kendall)


