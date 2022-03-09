
#use nco
import glob
list = glob.glob(r"I:\GPM\New folder\Jan\*.nc4")
from xarray import open_mfdataset
import numpy as np
import pandas as pd
#print(list)
import warnings
warnings.filterwarnings("ignore")
import sys
#sys.stdout = open(r'F:\GPM_ANALYSIS\REGION1_dec.txt', 'w')

# REGION1                # REGION2             # REGION3             # REGION4
# sel_lat = 25.05        # sel_lat = 27.05     # sel_lat = 25.05     # sel_lat = 22.05
# sel_lon = 89.05        # sel_lon = 92.05     # sel_lon = 92.05     # sel_lon = 91.05
# sel_lat = 26.95        # sel_lat = 28.95     # sel_lat = 26.95     # sel_lat = 24.95
# sel_lon = 91.95        # sel_lon = 96.95     # sel_lon = 95.95     # sel_lon = 94.95


#FOR timeseries in a region
for file in list:
    ds = open_mfdataset(file)
    latitude = ds.variables["lat"][:]
    longitude = ds.variables["lon"][:]
    sel_lat = 26.05
    sel_lon = 89.05
    a = abs(latitude-sel_lat)+abs(longitude-sel_lon)
    i,j = np.unravel_index(a.argmin(), a.shape)
    print(i)
    print(j)
    sel_lat = 26.95
    sel_lon = 91.95
    c = abs(latitude-sel_lat)+abs(longitude-sel_lon)
    k,l = np.unravel_index(c.argmin(), c.shape)
    print(k)
    print(l)
    precp = ds.variables["precipitationCal"][0, :, :]
    precp = precp.transpose("lat", "lon")
    precp = precp[i:k+1,j:l+1]
    #print(precp)
    precp= precp.values
    #print(precp)
    precp[precp <= -9999] = np.NAN
    precp= np.nansum(precp)
    print(precp)
