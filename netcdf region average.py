from xarray import open_dataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")


ds = open_dataset(r"F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\apr_rainfall_percentage_time_inc_or_dec.nc") #<----------------
rmse = np.nanmean(ds.variables["inc"][123:231, 182:299])
ioa = np.nanmean(ds.variables["dec"][123:231, 182:299])
#mae = np.nanmean(ds.variables["inc_>10"][123:231, 182:299])

print(rmse)
print(ioa)
#print(mae)