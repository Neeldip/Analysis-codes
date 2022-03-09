
#1.FIRST DO REGRIDDING OF MERRA2 DATA USING MERRA2 REGRIDDING.PY
#2.CDO REMAPBIL,TEMPLATE.NC INFILE OUTFILE

from xarray import open_dataset,open_mfdataset
import numpy as np
from netCDF4 import Dataset
from wrf import get_cartopy,getvar, ALL_TIMES,latlon_coords
from xarray import open_mfdataset
import matplotlib.pyplot as plt
import cartopy as cartopy
from cartopy.feature import NaturalEarthFeature
from wrf import omp_set_num_threads,get_cartopy,cartopy_xlim,cartopy_ylim
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter
from goodness_of_fit import rmse, d,me,r_pearson


#wrf_file = open_dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
#trad = open_dataset(r"G:\WRF_Chem_Output\202\NOR\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf_file = open_dataset(r"H:\WRF_Chem_Output\201\February\wrfout_d01_2018-02-10_00%3A00%3A00",engine='netcdf4')
trad = open_dataset(r"H:\WRF_Chem_Output\201\February\wrf_trad_fields_d01_2018-02-10_00%3A00%3A00",engine='netcdf4')

lat= wrf_file.variables["XLAT"][0,:,0]
lon = wrf_file.variables["XLONG"][0,0,:]
T = trad.variables["TEMPERATURE"][0:240,0,:,:] ##sensible temperature,units = K
P = trad.variables["PRESSURE"][0:240,0,:,:]  ##units = Pascal


###################BC###########################
BIN1 = wrf_file.variables["bc_a01"][0:240,0,:,:]
BIN2 = wrf_file.variables["bc_a02"][0:240,0,:,:]
BIN3 = wrf_file.variables["bc_a03"][0:240,0,:,:]
BIN4 = wrf_file.variables["bc_a04"][0:240,0,:,:]
R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 + BIN2 + BIN3 + BIN4
#BIN1=None
#BIN2=None
#BIN3=None
#BIN4=None
K = (R * T) / P  ### K unit = m^3*kg^-1 dry air
CONC_BC = X / K  ##unit = ug*m^-3

###################OC###########################
BIN1 = wrf_file.variables["oc_a01"][0:240,0,:,:]
BIN2 = wrf_file.variables["oc_a02"][0:240,0,:,:]
BIN3 = wrf_file.variables["oc_a03"][0:240,0,:,:]
BIN4 = wrf_file.variables["oc_a04"][0:240,0,:,:]
R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 + BIN2 + BIN3 + BIN4
#BIN1=None
#BIN2=None
#BIN3=None
#BIN4=None
K = (R * T) / P  ### K unit = m^3*kg^-1 dry air
CONC_OC = X / K  ##unit = ug*m^-3


###################dust###########################
BIN1 = wrf_file.variables["oin_a01"][0:240,0,:,:]
BIN2 = wrf_file.variables["oin_a02"][0:240,0,:,:]
BIN3 = wrf_file.variables["oin_a03"][0:240,0,:,:]
BIN4 = wrf_file.variables["oin_a04"][0:240,0,:,:]
R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 + BIN2 + BIN3 + BIN4
#BIN1=None
#BIN2=None
#BIN3=None
#BIN4=None
K = (R * T) / P  ### K unit = m^3*kg^-1 dry air
CONC_DUST = X / K  ##unit = ug*m^-3


###################SO4###########################
BIN1 = wrf_file.variables["so4_a01"][0:240,0,:,:]
BIN2 = wrf_file.variables["so4_a02"][0:240,0,:,:]
BIN3 = wrf_file.variables["so4_a03"][0:240,0,:,:]
BIN4 = wrf_file.variables["so4_a04"][0:240,0,:,:]
R = 287 ##gas constant, unit=J*kg^-1*K^-1
X = BIN1 + BIN2 + BIN3 + BIN4
#BIN1=None
#BIN2=None
#BIN3=None
#BIN4=None
K = (R * T) / P  ### K unit = m^3*kg^-1 dry air
CONC_SO4 = X / K  ##unit = ug*m^-3

BIN1=None
BIN2=None
BIN3=None
BIN4=None

ds1 =open_mfdataset([r"I:\Merra 2 M2T1NXAER\regridded\feb\10.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\11.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\12.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\13.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\14.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\15.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\16.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\17.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\18.nc4",
                     r"I:\Merra 2 M2T1NXAER\regridded\feb\19.nc4",])

RMSE_BC=np.empty((231,299))
ME_BC=np.empty((231,299))
IOA_BC=np.empty((231,299))
CORR_BC=np.empty((231,299))

RMSE_OC=np.empty((231,299))
ME_OC=np.empty((231,299))
IOA_OC=np.empty((231,299))
CORR_OC=np.empty((231,299))

RMSE_DUST=np.empty((231,299))
ME_DUST=np.empty((231,299))
IOA_DUST=np.empty((231,299))
CORR_DUST=np.empty((231,299))

RMSE_SO4=np.empty((231,299))
ME_SO4=np.empty((231,299))
IOA_SO4=np.empty((231,299))
CORR_SO4=np.empty((231,299))


for i in range(0,231,1):
    for j in range(0,299,1):
        bcsmass = ds1.variables["BCSMASS"][:, i, j] * 10 ** 9  ##CONVERT FROM KG/M3 TO UG/M3
        dusmass = ds1.variables["DUSMASS"][:, i, j] * 10 ** 9
        ocsmass = ds1.variables["OCSMASS"][:, i, j] * 10 ** 9
        so4smass = ds1.variables["SO4SMASS"][:, i, j] * 10 ** 9
        if np.isnan(bcsmass).any() == True:
           RMSE_BC[i,j] = np.nan
           IOA_BC[i, j] = np.nan
           ME_BC[i, j] = np.nan
           CORR_BC[i, j] = np.nan

           RMSE_OC[i,j] = np.nan
           IOA_OC[i, j] = np.nan
           ME_OC[i, j] = np.nan
           CORR_OC[i, j] = np.nan

           RMSE_DUST[i,j] = np.nan
           IOA_DUST[i, j] = np.nan
           ME_DUST[i, j] = np.nan
           CORR_DUST[i, j] = np.nan

           RMSE_SO4[i,j] = np.nan
           IOA_SO4[i, j] = np.nan
           ME_SO4[i, j] = np.nan
           CORR_SO4[i, j] = np.nan
           print(i)
           print(j)
        else:
           RMSE_BC[i, j] = rmse(CONC_BC[:,i, j],bcsmass)
           IOA_BC[i, j] = d(CONC_BC[:,i, j],bcsmass)
           ME_BC[i, j] = me(CONC_BC[:,i, j],bcsmass)
           CORR_BC[i, j] = r_pearson(CONC_BC[:,i, j],bcsmass)

           RMSE_OC[i, j] = rmse(CONC_OC[:,i, j],ocsmass)
           IOA_OC[i, j] = d(CONC_OC[:,i, j],ocsmass)
           ME_OC[i, j] = me(CONC_OC[:,i, j],ocsmass)
           CORR_OC[i, j] = r_pearson(CONC_OC[:,i, j],ocsmass)

           RMSE_DUST[i, j] = rmse(CONC_DUST[:,i, j],dusmass)
           IOA_DUST[i, j] = d(CONC_DUST[:,i, j],dusmass)
           ME_DUST[i, j] = me(CONC_DUST[:,i, j],dusmass)
           CORR_DUST[i, j] = r_pearson(CONC_DUST[:,i, j],dusmass)
           
           RMSE_SO4[i, j] = rmse(CONC_SO4[:,i, j],so4smass)
           IOA_SO4[i, j] = d(CONC_SO4[:,i, j],so4smass)
           ME_SO4[i, j] = me(CONC_SO4[:,i, j],so4smass)
           CORR_SO4[i, j] = r_pearson(CONC_SO4[:,i, j],so4smass)
           print(i)
           print(j)

print("BC")
print(np.nanmean(RMSE_BC))
print(np.nanmean(IOA_BC))
print(np.nanmean(ME_BC))
print(np.nanmean(CORR_BC))

print("OC")
print(np.nanmean(RMSE_OC))
print(np.nanmean(IOA_OC))
print(np.nanmean(ME_OC))
print(np.nanmean(CORR_OC))

print("DUST")
print(np.nanmean(RMSE_DUST))
print(np.nanmean(IOA_DUST))
print(np.nanmean(ME_DUST))
print(np.nanmean(CORR_DUST))

print("SO4")
print(np.nanmean(RMSE_SO4))
print(np.nanmean(IOA_SO4))
print(np.nanmean(ME_SO4))
print(np.nanmean(CORR_SO4))


CONC_BC=None
CONC_OC=None
CONC_DUST=None
CONC_SO4=None










'''
#ds2 = Dataset(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
ds2 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
RAINC231 = getvar(ds2, "RAINC", timeidx=0)
cart_proj = get_cartopy(RAINC231)
#ax = plt.axes(projection=ccrs.PlateCarree())
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
ax.set_xlim(cartopy_xlim(RAINC231))
ax.set_ylim(cartopy_ylim(RAINC231))


##for gridlines
'''
#gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
#                  linewidth=1, color='black', alpha=0.1, linestyle='-')
#gl.xlabels_top = False
#gl.xlabels_bottom = True
#gl.ylabels_left = True
#gl.ylabels_right = False
#gl.xlines = True
#gl.ylines = True
#gl.xlocator = mticker.FixedLocator([88, 90, 92, 94, 96, 98,100])
#gl.ylocator = mticker.FixedLocator([20, 22, 24, 26, 28, 30])
#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER
#gl.xlabel_style = {'size': 15, 'color': 'gray'}
#gl.xlabel_style = {'color': 'black', 'weight': 'normal'}
'''
#ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())

ax.set_xticks([75, 80, 85, 90, 95, 100], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())


#v = np.linspace(np.min(precp),np.max(precp),21)
#v = [-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0]
#v = [0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10]
#v = [-1.96,0,1.96]
#v= np.linspace(-30.0,30,11)
#v= np.linspace(0,50,21)
#v= np.linspace(-50.0,50.0,21)
#v= np.linspace(-0.006,0.006,21)
v= np.linspace(0,100,21)
plt.contourf(lon,lat,dusmass,v,cmap="jet",transform=cartopy.crs.PlateCarree(),extend='max',antialiased='False')
#plt.pcolormesh(lon,lat,precp,vmin=-1.96,vmax=1.96,transform=cartopy.crs.PlateCarree(),cmap="nipy_spectral")
plt.colorbar(ax=ax,ticks=v,shrink=0.96)
#plt.title("(d) decrease 0-0.21 mm/hr")

plt.savefig(r'D:\PHD\My PhD Reports\Manuscript3\Revision 1\Figures\dust_merra2_10.png',dpi=600,bbox_inches='tight')
plt.show()


'''

