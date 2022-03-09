import wrf as w
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from netCDF4 import Dataset
file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#file=Dataset(r"G:\WRF_Outputs\Monsoon\ACM2\wrfout_d03.nc")
HGT = w.getvar(file,'HGT',timeidx=0)
HGT=HGT/1000

lats= file.variables["XLAT"][0,:,0] #0-->only 1st step of time dimension
lons= file.variables["XLONG"][0,0,:]
#lats=np.flip(lats)

#print(lats)
X,Y=np.meshgrid(lats,lons,indexing='ij')
print(X)
print(Y)
fig = plt.figure()
ax = plt.axes(projection='3d')
#ax.set_zlim3d(0, 3000)
#ax.axis('equal')
ax.plot_surface(Y,X,HGT,cmap='jet')
plt.savefig("terrain.jpg",dpi=1200,bbox_inches='tight')
plt.show()
