from netCDF4 import Dataset
from xarray import open_mfdataset,DataArray
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
import pandas as pd
import wrf as w
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim

####1 DAY GPM DATA
#ds=open_mfdataset(r"I:\GPM\2018\July_new\*.nc4",combine='by_coords',concat_dim='time')  #https://neetinayak.medium.com/combine-many-netcdf-files-into-a-single-file-with-python-469ba476fc14
#precp = ds.variables["precipitationCal"][:,:,:]
#ds.to_netcdf('July_precip_combined.nc')

#*****************USE THESE FOR REGRIDDING***************
#ncpdq -a time,lat,lon precip_combined.nc precip_combined_reordered.nc
#cdo setgridtype,lonlat -seltimestep,1 -selvar,RAINC wrfout_d01_2018-04-10_00%3A00%3A00 RAINC.nc
#PLACE THE RAINC.nc in C:\Users\neeldip\PycharmProjects\WRF
#cdo -remapbil,RAINC.nc precip_combined_reordered.nc GPM.nc
#******************REGRIDDING ENDS***********************


ds = open_mfdataset('July_GPM_d03.nc')
lats = ds.variables["lat"][:]
lons = ds.variables["lon"][:]
ds1=ds.variables["precipitationCal"][:,:,:]
ds2 = ds1.sum(dim="time")
df =pd.DataFrame(ds2)
#print(ds2)

#wrf = Dataset("G:\WRF_Chem_Output\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf = Dataset("G:\WRF_Outputs\Monsoon\Hong\wrfout_d03.nc")
#RAINC= wrf.variables["RAINC"][744,:,:] - wrf.variables["RAINC"][24,:,:]
RAINC= wrf.variables["RAINC"][768,:,:] - wrf.variables["RAINC"][24,:,:]
#print(np.shape(RAINC))
#RAINNC= wrf.variables["RAINNC"][744,:,:] - wrf.variables["RAINNC"][24,:,:]
RAINNC= wrf.variables["RAINNC"][768,:,:] - wrf.variables["RAINNC"][24,:,:]
RAIN = RAINC +RAINNC
#print(RAIN)
RAIN1= RAIN - df.round()
#RAIN1 = np.asarray(RAIN1)
#print(RAIN1)
#v = np.linspace(np.min(RAIN1),np.max(RAIN1),30)


font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
RAINNC766_1 = w.getvar(wrf, "RAINNC", timeidx=0)
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


plt.contourf(lons,lats,RAIN1,30,cmap=plt.get_cmap("nipy_spectral"),antialiased=False,
             transform=cartopy.crs.PlateCarree())
plt.colorbar(ax=ax,shrink=0.81)
plt.title("HONG-GPM")
plt.savefig('HONG-GPM_JULY.png',dpi=100,bbox_inches='tight')
plt.show()
