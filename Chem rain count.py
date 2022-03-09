import pandas as pd
import wrf as w
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim,ll_to_xy
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from matplotlib import cbook
from netCDF4 import Dataset
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from mpl_toolkits.basemap import *

#ds = Dataset("G:\WRF_Chem_Output\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset("G:\WRF_Chem_Output\April\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset("G:\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
RAINC23 = w.getvar(ds, "RAINC", timeidx=0)
RAINNC23 = w.getvar(ds, "RAINNC", timeidx=0)
RAINC766 = w.getvar(ds, "RAINC", timeidx=239)
RAINNC766 = w.getvar(ds, "RAINNC", timeidx=239)
RAINC= RAINC766 - RAINC23
RAINNC =RAINNC766 - RAINNC23
RAIN = RAINC #+ RAINNC

#RAIN= (RAINC/RAIN)*100

#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
#ds1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
RAINC23_1 = w.getvar(ds1, "RAINC", timeidx=0)
RAINNC23_1 = w.getvar(ds1, "RAINNC", timeidx=0)
RAINC766_1 = w.getvar(ds1, "RAINC", timeidx=239)
RAINNC766_1 = w.getvar(ds1, "RAINNC", timeidx=239)
RAINC1=RAINC766_1 - RAINC23_1
RAINNC1 =RAINNC766_1 - RAINNC23_1
RAIN1=RAINC1 # + RAINNC1

#RAIN1= (RAINC1/RAIN1)*100
#index=ll_to_xy(ds1,latitude=30.3996734619,longitude=93.194915771484375)
#print(index)

RAIN = ((RAIN1-RAIN)/RAIN)*100

for i in range(0,231,1):
    for j in range(0,169,1):
        RAIN[i,j]=np.nan

for i in range(0,123,1):
    for j in range(169,273,1):
        RAIN[i,j]=np.nan

for i in range(215,231,1):
    for j in range(169,273,1):
        RAIN[i,j]=np.nan

for i in range(0,231,1):
    for j in range(273,299,1):
        RAIN[i,j]=np.nan

lats,lons =w.latlon_coords(RAINC23)

font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
cart_proj = get_cartopy(RAINNC766_1)
ax = plt.axes(projection=cart_proj)
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8)
#ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.add_feature(cartopy.feature.COASTLINE)
#ax.add_feature(cartopy.feature.LAKES)
#ax.add_feature(cartopy.feature.RIVERS)
ax.add_feature(cartopy.feature.BORDERS, linestyle='-')
ax.set_xlim(cartopy_xlim(RAINNC766_1))
ax.set_ylim(cartopy_ylim(RAINNC766_1))
#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(DIFF),np.max(DIFF),21) ##COLORBAR LIMITS
#v = np.linspace(-250,250,16) ##RAINNC LIMITS
#v = np.linspace(-100,100,16) ##RAINC LIMITS
#v = np.linspace(-300,300,16) ##TOTAL LIMITS
v = np.linspace(-100,100,16)
plt.contourf(lons, lats, RAIN,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
plt.colorbar(ax=ax, shrink=.81,ticks=v)
#plt.title("RAINCcontribution_change_aerfeedback - no_bc_absorbtion rainfall difference")
#plt.figure(figsize=(8,6))
plt.savefig(r'F:\WRF-CHEM ANALYSIS\rainfall\RAINC_2XBC-NORM%.png',dpi=600,bbox_inches='tight')
plt.show()