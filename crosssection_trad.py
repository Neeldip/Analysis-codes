import numpy as np
from numpy import transpose

from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES)
import matplotlib.pyplot as plt
#import sys
#np.set_printoptions(threshold=sys.maxsize)
#np.set_printoptions(threshold=np.inf)
from matplotlib.ticker import FormatStrFormatter

#wrf_file = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
'''
wrf_file = [Dataset(r"G:\WRF_Chem_Output\202new\March\wrfout_d01_2018-03-11_00%3A00%3A00"),
       Dataset(r"G:\WRF_Chem_Output\202new\March\wrfout_d01_2018-03-12_00%3A00%3A00"),
       Dataset(r"G:\WRF_Chem_Output\202new\March\wrfout_d01_2018-03-13_00%3A00%3A00"),
       Dataset(r"G:\WRF_Chem_Output\202new\March\wrfout_d01_2018-03-14_00%3A00%3A00"),
       Dataset(r"G:\WRF_Chem_Output\202new\March\wrfout_d01_2018-03-15_00%3A00%3A00"),
       ]
trad = [Dataset(r"G:\WRF_Chem_Output\202new\March\wrf_trad_fields_d01_2018-03-11_00%3A00%3A00"),
        Dataset(r"G:\WRF_Chem_Output\202new\March\wrf_trad_fields_d01_2018-03-12_00%3A00%3A00"),
        Dataset(r"G:\WRF_Chem_Output\202new\March\wrf_trad_fields_d01_2018-03-13_00%3A00%3A00"),
        Dataset(r"G:\WRF_Chem_Output\202new\March\wrf_trad_fields_d01_2018-03-14_00%3A00%3A00"),
        Dataset(r"G:\WRF_Chem_Output\202new\March\wrf_trad_fields_d01_2018-03-15_00%3A00%3A00")]
'''
'''
wrf_file = [Dataset(r"G:\WRF_Chem_Output\202new\April\wrfout_d01_2018-04-04_00%3A00%3A00"),
       Dataset(r"G:\WRF_Chem_Output\202new\April\wrfout_d01_2018-04-05_00%3A00%3A00"),
       Dataset(r"G:\WRF_Chem_Output\202new\April\wrfout_d01_2018-04-06_00%3A00%3A00"),
       Dataset(r"G:\WRF_Chem_Output\202new\April\wrfout_d01_2018-04-07_00%3A00%3A00"),
       Dataset(r"G:\WRF_Chem_Output\202new\April\wrfout_d01_2018-04-08_00%3A00%3A00"),
       ]
trad = [Dataset(r"G:\WRF_Chem_Output\202new\April\wrf_trad_fields_d01_2018-04-04_00%3A00%3A00"),
        Dataset(r"G:\WRF_Chem_Output\202new\April\wrf_trad_fields_d01_2018-04-05_00%3A00%3A00"),
        Dataset(r"G:\WRF_Chem_Output\202new\April\wrf_trad_fields_d01_2018-04-06_00%3A00%3A00"),
        Dataset(r"G:\WRF_Chem_Output\202new\April\wrf_trad_fields_d01_2018-04-07_00%3A00%3A00"),
        Dataset(r"G:\WRF_Chem_Output\202new\April\wrf_trad_fields_d01_2018-04-08_00%3A00%3A00")]
'''
#trad1 = Dataset("G:\WRF_Chem_Output\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#trad1 = Dataset(r"G:\WRF_Chem_Output\April\MYNN3_BC_no_absorbtion\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf_file = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
#trad = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrf_trad_d03.nc")
trad = Dataset("G:\WRF_Outputs\Monsoon\ACM2\wrf_trad_d03.nc")

#trad = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrf_trad_d03.nc")
#trad = Dataset("G:\WRF_Outputs\Monsoon\MYNN3\wrf_trad_d03.nc")

# Define the cross section start and end points
#cross_start = CoordPair(lat=24, lon=91.73)
#cross_end = CoordPair(lat=27.5, lon=91.73)

cross_start = CoordPair(lat=25.98, lon=89.84)
cross_end = CoordPair(lat=27.86, lon=95.94)
time= ALL_TIMES
#Get the WRF variables
ht = getvar(wrf_file, "z", timeidx=0)
ter = getvar(wrf_file, "ter", timeidx=0)
BIN1 = getvar(trad, "RH", timeidx=time)
#BIN2 = getvar(trad1, "TEMPERATURE", timeidx=time)

X = BIN1
#Y = BIN2

if time == ALL_TIMES:
    X = X.mean("Time")
    #Y = Y.mean("Time")

CONC = X #- Y
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)


dbz_cross_filled = np.ma.copy(to_np(z_cross))

for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]

# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)


fig = pyplot.figure(figsize=(13, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), np.max(dbz_cross_filled), 20, endpoint=True)
#v = np.linspace(0, 20, 20, endpoint=True)
v = np.linspace(0, 100, 20, endpoint=True)
#v = np.linspace(np.min(dbz_cross_filled), 0.02013, 45, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
#QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="nipy_spectral",extend='max')
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="nipy_spectral")
ax_cross.tick_params(labelsize=12)

# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v,extend='max')
cb_dbz.set_label("Relative humidity (%)",size=14)
cb_dbz.ax.tick_params(labelsize=14)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
#print(z_cross.coords["xy_loc"])
coord_pairs = to_np(z_cross.coords["xy_loc"])
x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]
#print(x_labels)
#x_labels= np.around(x_labels,decimals=2)
#print(x_labels)



# Set the desired number of x ticks below
#num_ticks = 5
#thin = int((len(x_ticks) / num_ticks) + .5)
#ax_cross.set_xticks(x_ticks[::thin])
#ax_cross.set_xticklabels(x_labels[::thin], rotation=0, fontsize=12)
#ax_cross.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax_cross.set_xticks(x_ticks[::thin])
#ax_cross.set_xticklabels(x_labels[::thin], rotation=0, fontsize=12)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

# Set the x-axis and  y-axis labels
#ax_cross.set_xlabel("Latitude, Longitude", fontsize=14)
ax_cross.set_ylabel("Height (m)", fontsize=15)
ax_cross.set_ylim(0, 5000, {"fontsize":15})

# Add a title
ax_cross.set_title("(d) July-ACM2", {"fontsize": 16})
#pyplot.savefig('TEMPERATURE(aer_feedback - bc_no_absorbtion) (K)'+str(cross_start)+str(cross_end)+'.png',dpi=500)
plt.savefig('RH-july-EW-acm2.png',dpi=300,bbox_inches='tight')
plt.show()


