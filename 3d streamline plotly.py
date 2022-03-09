import numpy as np
from numpy import *
from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES, enable_pyngl)

wrf_file = Dataset("G:\WRF_Chem_Output\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
# wrf_file = Dataset("G:\WRF_Chem_Output\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
time = 1
U = getvar(wrf_file,"ua",timeidx=time)
V = getvar(wrf_file,"va",timeidx=time)
W = getvar(wrf_file,"wa",timeidx=time)
