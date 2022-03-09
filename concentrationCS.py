import numpy as np
from numpy import transpose

from matplotlib import pyplot
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES)
#import sys
#np.set_printoptions(threshold=sys.maxsize)
#np.set_printoptions(threshold=np.inf)

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\February\wrfout_d01_2018-02-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\February\wrf_trad_fields_d01_2018-02-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\March\wrf_trad_fields_d01_2018-03-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\May\wrfout_d01_2018-05-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\May\wrf_trad_fields_d01_2018-05-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\June\wrfout_d01_2018-06-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\June\wrf_trad_fields_d01_2018-06-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\July\wrfout_d01_2018-07-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\July\wrf_trad_fields_d01_2018-07-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\August\wrfout_d01_2018-08-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\August\wrf_trad_fields_d01_2018-08-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\September\wrfout_d01_2018-09-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\September\wrf_trad_fields_d01_2018-09-10_00%3A00%3A00")

#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_No_NE_BC\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_No_NE_BC\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf_file =  Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_Only_NE_BC\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_Only_NE_BC\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

time= ALL_TIMES
#Get the WRF variables
ht = getvar(wrf_file, "z", timeidx=0)
ter = getvar(wrf_file, "ter", timeidx=0)
BIN1 = getvar(wrf_file, "bc_a01", timeidx=time)
BIN2 = getvar(wrf_file, "bc_a02", timeidx=time)
BIN3 = getvar(wrf_file, "bc_a03", timeidx=time)
BIN4 = getvar(wrf_file, "bc_a04", timeidx=time)
T = getvar(trad, "TEMPERATURE", timeidx=time)  ##sensible temperature,units = K
P = getvar(trad, "PRESSURE", timeidx=time)
X = BIN1+BIN2+BIN3+BIN4
BIN1=None
BIN2=None
BIN3=None
BIN4=None
R = 287                 ##Universal gas constant,unit=J*kg^-1*K^-1
if time == ALL_TIMES:
    X = X.mean("Time")
    T = T.mean("Time")
    P = P.mean("Time")

K = (R*T)/P
CONC = X/K

# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=91.0)
cross_end = CoordPair(lat=30.0, lon=91.0)
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
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
#pyplot.savefig(r'D:\PHD\My PhD Reports\IICAQM conference\apr_no ne 2xbc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
pyplot.savefig(r'F:\WRF-CHEM ANALYSIS Chap 2\concentration and fraction\Only_NE_BC'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=v)
pyplot.show()



# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=92.0)
cross_end = CoordPair(lat=30.0, lon=92.0)
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#X=None
#T=None
#P=None
#K=None
dbz_cross_filled = np.ma.copy(to_np(z_cross))
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
pyplot.savefig(r'F:\WRF-CHEM ANALYSIS Chap 2\concentration and fraction\Only_NE_bc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=np.linspace(0,3,21))
pyplot.show()






# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=93.0)
cross_end = CoordPair(lat=30.0, lon=93.0)
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#X=None
#T=None
#P=None
#K=None
dbz_cross_filled = np.ma.copy(to_np(z_cross))
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
pyplot.savefig(r'F:\WRF-CHEM ANALYSIS Chap 2\concentration and fraction\Only_NE_bc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=np.linspace(0,3,21))
pyplot.show()





# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=94.0)
cross_end = CoordPair(lat=30.0, lon=94.0)
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#X=None
#T=None
#P=None
#K=None
dbz_cross_filled = np.ma.copy(to_np(z_cross))
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
pyplot.savefig(r'F:\WRF-CHEM ANALYSIS Chap 2\concentration and fraction\Only_NE_bc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=np.linspace(0,3,21))
pyplot.show()






# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=95.0)
cross_end = CoordPair(lat=30.0, lon=95.0)
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#X=None
#T=None
#P=None
#K=None
dbz_cross_filled = np.ma.copy(to_np(z_cross))
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
pyplot.savefig(r'F:\WRF-CHEM ANALYSIS Chap 2\concentration and fraction\Only_NE_bc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=np.linspace(0,3,21))
pyplot.show()




'''

# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=93.0)
cross_end = CoordPair(lat=30.0, lon=93.0)
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#X=None
#T=None
#P=None
#K=None
dbz_cross_filled = np.ma.copy(to_np(z_cross))
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
pyplot.savefig(r'D:\PHD\My PhD Reports\IICAQM conference\apr_no ne 2xbc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=np.linspace(0,3,21))
pyplot.show()





# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=94.0)
cross_end = CoordPair(lat=30.0, lon=94.0)
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#X=None
#T=None
#P=None
#K=None
dbz_cross_filled = np.ma.copy(to_np(z_cross))
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
pyplot.savefig(r'D:\PHD\My PhD Reports\IICAQM conference\apr_no ne 2xbc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=np.linspace(0,3,21))
pyplot.show()




# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=95.0)
cross_end = CoordPair(lat=30.0, lon=95.0)
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#X=None
#T=None
#P=None
#K=None
dbz_cross_filled = np.ma.copy(to_np(z_cross))
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
pyplot.savefig(r'D:\PHD\My PhD Reports\IICAQM conference\apr_no ne 2xbc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=np.linspace(0,3,21))
pyplot.show()




# Define the cross section start and end points
cross_start = CoordPair(lat=22.5, lon=96.0)
cross_end = CoordPair(lat=30.0, lon=96.0)
z_cross = vertcross(CONC, ht, wrfin=wrf_file,start_point=cross_start,end_point=cross_end,latlon=True, meta=True)
#X=None
#T=None
#P=None
#K=None
dbz_cross_filled = np.ma.copy(to_np(z_cross))
for i in range(dbz_cross_filled.shape[-1]):
    column_vals = dbz_cross_filled[:, i]
    first_idx = int(np.transpose((column_vals > -10).nonzero())[0])
    dbz_cross_filled[0:first_idx, i] = dbz_cross_filled[first_idx, i]
# Get the terrain heights along the cross section line
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,end_point=cross_end)

# Get the lat/lon points
lats, lons = latlon_coords(CONC)
fig = pyplot.figure(figsize=(15, 8))
ax_cross = pyplot.axes()
#v = np.linspace(np.min(dbz_cross_filled), 30, 45, endpoint=True)
v = np.linspace(0, 2 ,21, endpoint=True)
xs = np.arange(0, z_cross.shape[-1], 1)
ys = to_np(z_cross.coords["vertical"])
QCLOUD_contours = ax_cross.contourf(xs, ys, to_np(dbz_cross_filled),v, cmap="gist_ncar",extend='max')
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
#ax_cross.set_title("BC_concentration(aerfeebback) (ug m^-3) ", {"fontsize": 16})
pyplot.savefig(r'D:\PHD\My PhD Reports\IICAQM conference\apr_no ne 2xbc_CONCENTRATION (UG M3)-1'+str(cross_start)+str(cross_end)+'.png',dpi=600,bbox_inches='tight')
#pyplot.colorbar(ticks=np.linspace(0,3,21))
pyplot.show()


X=None
T=None
P=None
K=None
'''

