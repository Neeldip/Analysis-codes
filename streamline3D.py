#https://docs.enthought.com/mayavi/mayavi/auto/mlab_helper_functions.html#flow
from mayavi import mlab
from numpy import random
from netCDF4 import Dataset
from wrf import getvar,ALL_TIMES
@mlab.show
def image():
    mlab.imshow(random.random((1, 1)))
import numpy as np
from mayavi.mlab import *
wrf_file = Dataset("G:\WRF_Chem_Output\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
U = getvar(wrf_file,"ua",timeidx=1)
V = getvar(wrf_file,"va",timeidx=1)
W = getvar(wrf_file,"wa",timeidx=1)
#U = U.mean("Time")
#V = V.mean("Time")
#W = w.mean("Time")
#x, y, z = np.mgrid[-4:4:40j, -4:4:40j, 0:4:20j]
#r = np.sqrt(x ** 2 + y ** 2 + z ** 2 + 0.1)
#u = y * np.sin(r) / r
#v = -x * np.sin(r) / r
#w = np.ones_like(z)*0.05
obj = flow(U, V, W)
show()