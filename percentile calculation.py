from xarray import open_dataset
import numpy as np
import warnings
warnings.filterwarnings("ignore")


ds = open_dataset(r"F:\WRF-CHEM ANALYSIS\WRF-Chem rainfall evaluation\apr_rainfall_evaluation_nor.nc") #<----------------
rmse = ds.variables["rmse"][123:231, 182:299]
ioa = ds.variables["ioa"][123:231, 182:299]
perc=90
ioa_perc= np.nanpercentile(ioa,perc)
rmse_perc= np.nanpercentile(rmse,100-perc)
print(ioa_perc)
print(rmse_perc)
