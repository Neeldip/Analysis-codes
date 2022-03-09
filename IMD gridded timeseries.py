import glob
list = glob.glob(r"I:\observation_data\APRIL\DATA\Rainfall\IMD_gridded_datasets\RAW data\NETCDF\*.nc")
from xarray import open_mfdataset
import numpy as np
import pandas as pd
#print(list)
import warnings
warnings.filterwarnings("ignore")
import sys
np.set_printoptions(threshold=sys.maxsize)
sys.stdout = open(r'F:\IMD\REGION4.txt', 'w')


#*****************concatenate************************
#ncrcat * IMD.nc


#For region 1
#lats=[26.32,26.40,26.47,26.18,26.31,26.44,26.17,25.28,25.56,26.02]
#lons=[89.45,90.27,90.54,90.62,91.00,91.44,91.75,91.56,91.87,89.97]

#For region 2
#lats=[27.22,27.48,27.48,27.28,27.92,28.06,27.47,27.58,27.66]
#lons=[94.09,94.92,95.36,95.65,96.14,95.32,94.55,91.86,95.86]

#For region 3
#lats=[26.34,26.62,25.91,26.76,25.65,26.98]
#lons=[92.68,92.78,93.72,94.21,94.09,94.63]

#For region 4
lats=[24.80,23.83,23.74,24.81,22.91]
lons=[93.93,91.28,92.73,92.75,92.74]

for sel_lat,sel_lon in zip(lats,lons):
    #print(sel_lat,sel_lon)
#FOR timeseries at a point
    for file in list:
        ds = open_mfdataset(file)
        latitude = ds.variables["LATITUDE"][:]
        longitude = ds.variables["LONGITUDE"][:]
        a = abs(latitude-sel_lat)+abs(longitude-sel_lon)
        i,j = np.unravel_index(a.argmin(), a.shape)
        precp = np.asarray(ds.variables["RAINFALL"][:, i, j])
        #print(precp.shape)
        precp = np.reshape(precp,(7305,1))
        #print(precp)

    print("\n \n \n \n \n \n \n \n")
np.savetxt(sys.stdout.buffer, precp,fmt="%.3f")
sys.stdout.close()
