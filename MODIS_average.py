from xarray import open_dataset
import numpy as np
from wrf import get_cartopy,getvar,cartopy_xlim,cartopy_ylim
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
ds= open_dataset(r"F:\AOD_ANALYSIS\REGRIDDING MYD_L3\jan.nc")

AOD = ds.variables["AOD"][:,:,:]

lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]

AOD = np.nanmean(AOD,axis=0)


#wrf = Dataset(r"G:\WRF_Outputs\Monsoon\ACM2\wrfout_d02_2018-07-01_00%3A00%3A00")
wrf = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
font = {'family': 'serif',
        'color':  'black',
        'weight': 'bold',
        'size': 12,
        }
RAINNC766_1 = getvar(wrf, "RAINNC", timeidx=0)
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
v = np.linspace(0,1,20)
#plt.contourf(lon,lat,precp,v,cmap="nipy_spectral",orientation="vertical",
#             transform=cartopy.crs.PlateCarree())
#plt.pcolormesh(lon,lat,precp,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
plt.pcolormesh(lon,lat,AOD,vmin=0,vmax=1,transform=cartopy.crs.PlateCarree(),cmap="jet")
plt.colorbar(ticks=v,shrink=0.93)
plt.title("January")  #<----------------
plt.savefig(r'F:\AOD_ANALYSIS\MODIS_MOD_L3\mean_jan.png',dpi=300,bbox_inches='tight')  #<----------------
plt.show()