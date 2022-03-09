import numpy as np
from numpy import transpose

from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES)
#import sys
#np.set_printoptions(threshold=sys.maxsize)
#np.set_printoptions(threshold=np.inf)

#ACM2
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

##MYNN3
#wrf_file = Dataset("G:\WRF_Chem_Output\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset("G:\WRF_Chem_Output\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
'''
wrf_file = [Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-10_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-11_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-12_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-13_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-14_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-15_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-16_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-17_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-18_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-19_00%3A00%3A00"),
            Dataset(r"H:WRF_Chem_Output\201\Sep\wrfout_d01_2018-09-20_00%3A00%3A00")]
'''
'''
trad = [Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-10_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-11_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-12_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-13_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-14_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-15_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-16_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-17_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-18_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-19_00%3A00%3A00"),
        Dataset(r"H:\WRF_Chem_Output\201\Sep\wrf_trad_fields_d01_2018-09-20_00%3A00%3A00"),]
'''
#'''
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

#Select time
time= ALL_TIMES
#time = ALL_TIMES
fig = True
t = 241   ##no of time steps in wrf file

##Select wind component to use
U = getvar(wrf_file,"ua",timeidx=time)
#U = getvar(wrf_file,"va",timeidx=time)
V = getvar(wrf_file,"va",timeidx=time)

# Get the WRF variables
ht = getvar(wrf_file, "z", timeidx=0)
ter = getvar(wrf_file, "ter", timeidx=0)
BIN1 = getvar(wrf_file, "bc_a01", timeidx=time)
BIN2 = getvar(wrf_file, "bc_a02", timeidx=time)
BIN3 = getvar(wrf_file, "bc_a03", timeidx=time)
BIN4 = getvar(wrf_file, "bc_a04", timeidx=time)
T = getvar(trad, "TEMPERATURE", timeidx=time)  ##sensible temperature,units = K
P = getvar(trad, "PRESSURE", timeidx=time)     ##units = Pascal

##Sum all bins
X = BIN1+BIN2+BIN3+BIN4
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1

##For average plot
if time == ALL_TIMES:
 X = X.mean("Time")
 U = U.mean("Time")
 V = V.mean("Time")
 T = T.mean("Time")
 P = P.mean("Time")

#Calculation
SPD= np.sqrt(U**2+V**2)
K = (R*T)/P                     ##Specific volume of dry air
y = (SPD*X)/K                     ##For single time or mean,use this, unit = ug/m^2*s
#y = (X*U*.000001*t*3600)/K      ##For total flux use this, unit = g/m^2

###################################################PLOT STARTS#################################################

# Define the cross section start and end points for plot 1
cross_start = CoordPair(lat=21, lon=90)
cross_end = CoordPair(lat=28, lon=90)


#################CROSS SECTION ###################
z_cross = vertcross(y, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#z_cross.plot()
#print(np.ndim(z_cross))
# To remove the slight gap between the dbz contours and terrain due to the
# contouring of gridded data, a new vertical grid spacing, and model grid
# staggering, fill in the lower grid cells with the first non-missing value
# for each column.

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
dbz_cross_filled = np.ma.copy(to_np(z_cross))
print(np.ndim(dbz_cross_filled))
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
print(np.mean(dbz_cross_filled))
print(np.average(dbz_cross_filled))
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(y)


# Create the figure
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
v = np.linspace(np.min(dbz_cross_filled), np.max(dbz_cross_filled), 30, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="bwr")
ax_cross.tick_params(labelsize=12)
# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v)
cb_dbz.ax.tick_params(labelsize=12)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(z_cross.coords["xy_loc"])
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
ax_cross.set_title("flux (ug m^-2 s-1) ", {"fontsize": 16})
#ax_cross.set_title("QCLOUD (g/kg) ", {"fontsize" : 14})

pyplot.savefig("flux_time="+str(time)+"coord="+str(cross_start)+"_"+str(cross_end)+".png",dpi=300,bbox_inches='tight')
pyplot.show()
#'''

# Define the cross section start and end points for plot 2
cross_start = CoordPair(lat=24.5, lon=85)
cross_end = CoordPair(lat=24.5, lon=96)

#################CROSS SECTION ###################
z_cross = vertcross(y, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#z_cross.plot()
print(np.ndim(z_cross))
# To remove the slight gap between the dbz contours and terrain due to the
# contouring of gridded data, a new vertical grid spacing, and model grid
# staggering, fill in the lower grid cells with the first non-missing value
# for each column.

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
dbz_cross_filled = np.ma.copy(to_np(z_cross))
print(np.ndim(dbz_cross_filled))
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
lats, lons = latlon_coords(y)


# Create the figure
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
v = np.linspace(np.min(dbz_cross_filled), np.max(dbz_cross_filled), 30, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="bwr")
ax_cross.tick_params(labelsize=12)
# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v)
cb_dbz.ax.tick_params(labelsize=12)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(z_cross.coords["xy_loc"])
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
ax_cross.set_ylim(0, 16000, {"fontsize":14})

# Add a title
ax_cross.set_title("flux (ug m^-2 s-1) ", {"fontsize": 16})
#ax_cross.set_title("QCLOUD (g/kg) ", {"fontsize" : 14})
pyplot.savefig("flux_time="+str(time)+"coord="+str(cross_start)+"_"+str(cross_end)+".png",dpi=300,bbox_inches='tight')
pyplot.show()


# Define the cross section start and end points for plot 3
cross_start = CoordPair(lat=26.1, lon=85)
cross_end = CoordPair(lat=26.1, lon=96)


#################CROSS SECTION ###################
z_cross = vertcross(y, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#z_cross.plot()
print(np.ndim(z_cross))
# To remove the slight gap between the dbz contours and terrain due to the
# contouring of gridded data, a new vertical grid spacing, and model grid
# staggering, fill in the lower grid cells with the first non-missing value
# for each column.

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
dbz_cross_filled = np.ma.copy(to_np(z_cross))
print(np.ndim(dbz_cross_filled))
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
lats, lons = latlon_coords(y)


# Create the figure
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
v = np.linspace(np.min(dbz_cross_filled), np.max(dbz_cross_filled), 30, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="bwr")
ax_cross.tick_params(labelsize=12)
# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v)
cb_dbz.ax.tick_params(labelsize=12)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(z_cross.coords["xy_loc"])
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
ax_cross.set_ylim(0, 16000, {"fontsize":14})

# Add a title
ax_cross.set_title("flux (ug m^-2 s-1) ", {"fontsize": 16})
#ax_cross.set_title("QCLOUD (g/kg) ", {"fontsize" : 14})
pyplot.savefig("flux_time="+str(time)+"coord="+str(cross_start)+"_"+str(cross_end)+".png",dpi=300)
pyplot.show()

# Define the cross section start and end points for plot 4
cross_start = CoordPair(lat=27, lon=85)
cross_end = CoordPair(lat=27, lon=96)


#################CROSS SECTION ###################
z_cross = vertcross(y, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#z_cross.plot()
print(np.ndim(z_cross))
# To remove the slight gap between the dbz contours and terrain due to the
# contouring of gridded data, a new vertical grid spacing, and model grid
# staggering, fill in the lower grid cells with the first non-missing value
# for each column.

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
dbz_cross_filled = np.ma.copy(to_np(z_cross))
print(np.ndim(dbz_cross_filled))
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
lats, lons = latlon_coords(y)


# Create the figure
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
v = np.linspace(np.min(dbz_cross_filled), np.max(dbz_cross_filled), 30, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="bwr")
ax_cross.tick_params(labelsize=12)
# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v)
cb_dbz.ax.tick_params(labelsize=12)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(z_cross.coords["xy_loc"])
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
ax_cross.set_ylim(0, 16000, {"fontsize":14})

# Add a title
ax_cross.set_title("flux (ug m^-2 s-1) ", {"fontsize": 16})
#ax_cross.set_title("QCLOUD (g/kg) ", {"fontsize" : 14})
pyplot.savefig("flux_time="+str(time)+"coord="+str(cross_start)+"_"+str(cross_end)+".png",dpi=300,bbox_inches='tight')
pyplot.show()


# Define the cross section start and end points for plot 5
cross_start = CoordPair(lat=28.16, lon=85)
cross_end = CoordPair(lat=28.16, lon=96)


#################CROSS SECTION ###################
z_cross = vertcross(y, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#z_cross.plot()
print(np.ndim(z_cross))
# To remove the slight gap between the dbz contours and terrain due to the
# contouring of gridded data, a new vertical grid spacing, and model grid
# staggering, fill in the lower grid cells with the first non-missing value
# for each column.

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
dbz_cross_filled = np.ma.copy(to_np(z_cross))
print(np.ndim(dbz_cross_filled))
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
lats, lons = latlon_coords(y)


# Create the figure
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
v = np.linspace(np.min(dbz_cross_filled), np.max(dbz_cross_filled), 30, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="bwr")
ax_cross.tick_params(labelsize=12)
# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v)
cb_dbz.ax.tick_params(labelsize=12)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(z_cross.coords["xy_loc"])
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
ax_cross.set_ylim(0, 16000, {"fontsize":14})

# Add a title
ax_cross.set_title("flux (ug m^-2 s-1) ", {"fontsize": 16})
#ax_cross.set_title("QCLOUD (g/kg) ", {"fontsize" : 14})
pyplot.savefig("flux_time="+str(time)+"coord="+str(cross_start)+"_"+str(cross_end)+".png",dpi=300,bbox_inches='tight')
pyplot.show()

# Define the cross section start and end points for plot 6
cross_start = CoordPair(lat=21, lon=96)
cross_end = CoordPair(lat=28.16, lon=96)


#################CROSS SECTION ###################
z_cross = vertcross(y, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#z_cross.plot()
print(np.ndim(z_cross))
# To remove the slight gap between the dbz contours and terrain due to the
# contouring of gridded data, a new vertical grid spacing, and model grid
# staggering, fill in the lower grid cells with the first non-missing value
# for each column.

# Make a copy of the z cross data. Let's use regular numpy arrays for this.
dbz_cross_filled = np.ma.copy(to_np(z_cross))
print(np.ndim(dbz_cross_filled))
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
lats, lons = latlon_coords(y)


# Create the figure
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
v = np.linspace(np.min(dbz_cross_filled), np.max(dbz_cross_filled), 30, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="bwr")
ax_cross.tick_params(labelsize=12)
# Add the color bar
cb_dbz = fig.colorbar(QCLOUD_contours, ax=ax_cross,ticks=v)
cb_dbz.ax.tick_params(labelsize=12)

# Fill in the mountain area
ht_fill = ax_cross.fill_between(xs, 0, to_np(ter_line),facecolor="saddlebrown")

# Set the x-ticks to use latitude and longitude labels
coord_pairs = to_np(z_cross.coords["xy_loc"])
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
ax_cross.set_ylim(0, 16000, {"fontsize":14})

# Add a title
ax_cross.set_title("flux (ug m^-2 s-1) ", {"fontsize": 16})
#ax_cross.set_title("QCLOUD (g/kg) ", {"fontsize" : 14})
pyplot.savefig("flux_time="+str(time)+"coord="+str(cross_start)+"_"+str(cross_end)+".png",dpi=300,bbox_inches='tight')
pyplot.show()


#'''






