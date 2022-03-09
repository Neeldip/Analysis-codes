import numpy as np
from numpy import transpose

from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES,destagger)
#import sys
#np.set_printoptions(threshold=sys.maxsize)
#np.set_printoptions(threshold=np.inf)


#wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrfout_d03.nc")
# Define the cross section start and end points
cross_start = CoordPair(lat=22.0, lon=91.74)
cross_end = CoordPair(lat=28.0, lon=91.74)
time= 0

#Get the WRF variables
ht = getvar(wrf_file, "z", timeidx=0)
print(ht)
ter = getvar(wrf_file, "ter", timeidx=0)
print(ter)
BIN1 = getvar(wrf_file, "W", timeidx=time)
RAINC = getvar(wrf_file,"RAINC",timeidx=0)
BIN1 = destagger(BIN1, stagger_dim=0, meta=False)
X = BIN1


if time == ALL_TIMES:
    X = X.mean("Time")

CONC = X*100
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#z_cross.plot()
print(z_cross)
# To remove the slight gap between the dbz contours and terrain due to the
# contouring of gridded data, a new vertical grid spacing, and model grid
# staggering, fill in the lower grid cells with the first non-missing value
# for each column.

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
dbz_cross_filled = np.ma.copy(to_np(z_cross))
#print(dbz_cross_filled)
#np.savetxt("BCflux.csv",dbz_cross_filled,delimiter=",")
# For each cross section column, find the first index with non-missing
# values and copy these to the missing elements below.
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    # Let's find the lowest index that isn't filled. The nonzero function
    # finds all unmasked values greater than 0. Since 0 is a valid value
    # for dBZ, let's change that threshold to be -200 dBZ instead.
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]

# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(RAINC)


fig = pyplot.figure(figsize=(20, 8))
ax_cross = pyplot.axes()
v = np.linspace(np.min(dbz_cross_filled), np.max(dbz_cross_filled), 45, endpoint=True)
#v = np.linspace(np.min(dbz_cross_filled), 0.02013, 45, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["dim_0"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="nipy_spectral")
ax_cross.tick_params(labelsize=12)
# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v)
cb_dbz.ax.tick_params(labelsize=12)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(z_cross.coords["dim_1"])
x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]

# Set the desired number of x ticks below
num_ticks = 5
thin = int((len(x_ticks) / num_ticks) + .5)
ax_cross.set_xticks(x_ticks[::thin])
ax_cross.set_xticklabels(x_labels[::thin], rotation=0, fontsize=12)

# Set the x-axis and  y-axis labels
ax_cross.set_xlabel("Latitude, Longitude", fontsize=14)
ax_cross.set_ylabel("Height (m)", fontsize=15)
ax_cross.set_ylim(0, 10000, {"fontsize":14})

# Add a title
ax_cross.set_title("THM (K) ", {"fontsize": 16})
#pyplot.savefig("THM_hong_"+str(cross_start)+str(cross_end)+".png",dpi=500)
pyplot.show()


