import numpy as np
from matplotlib import pyplot
from netCDF4 import Dataset

from wrf import (getvar, to_np, latlon_coords, vertcross,
                 interpline, CoordPair, ALL_TIMES)

wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")

# Define the cross section start and end points
cross_start = CoordPair(lat=24.0, lon=88)
cross_end = CoordPair(lat=28.16, lon=88)

# Get the WRF variables
ht = getvar(wrf_file, "z", timeidx=-1)
ter = getvar(wrf_file, "ter", timeidx=-1)
QCLOUD= getvar(wrf_file, "ua", timeidx=ALL_TIMES)
QCLOUD= QCLOUD.mean("Time")
print(np.ndim(QCLOUD))
# Compute the vertical cross-section interpolation.  Also, include the
# lat/lon points along the cross-section in the metadata by setting latlon
# to True.
z_cross = vertcross(QCLOUD, ht, wrfin=wrf_file,
                    start_point=cross_start,
                    end_point=cross_end,
                    latlon=True, meta=True)
z_cross.plot()
print(np.ndim(z_cross))
# To remove the slight gap between the dbz contours and terrain due to the
# contouring of gridded data, a new vertical grid spacing, and model grid
# staggering, fill in the lower grid cells with the first non-missing value
# for each column.

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
dbz_cross_filled = np.ma.copy(to_np(z_cross))
print(np.ndim(dbz_cross_filled))
# For each cross section column, find the first index with non-missing
# values and copy these to the missing elements below.
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:,i]
    # Let's find the lowest index that isn't filled. The nonzero function
    # finds all unmasked values greater than 0. Since 0 is a valid value
    # for dBZ, let's change that threshold to be -200 dBZ instead.
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]

# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,
                      end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(QCLOUD)

# Get the cartopy projection object
#cart_proj = get_cartopy(QCLOUD)

# Create the figure
fig = pyplot.figure(figsize=(20,6))
ax_cross = pyplot.axes()
v = np.linspace(-28, 40, 32, endpoint=True)
xs = np.arange(0,z_cross.shape[-1],1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs,ys,to_np(dbz_cross_filled),v,cmap="nipy_spectral")
# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v)
cb_dbz.ax.tick_params(labelsize=8)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),
                                facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(z_cross.coords["xy_loc"])
x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]

# Set the desired number of x ticks below
num_ticks = 5
thin = int((len(x_ticks) / num_ticks) + .5)
ax_cross.set_xticks(x_ticks[::thin])
ax_cross.set_xticklabels(x_labels[::thin], rotation=45, fontsize=8)

# Set the x-axis and  y-axis labels
ax_cross.set_xlabel("Latitude, Longitude", fontsize=12)
ax_cross.set_ylabel("Height (m)", fontsize=12)
ax_cross.set_ylim(0,16000)

# Add a title
ax_cross.set_title("U:along EAST_WEST (m/s) ", {"fontsize" : 14})
pyplot.savefig('U.png',dpi=500)
pyplot.show()