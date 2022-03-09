import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER,mticker,LongitudeFormatter,LatitudeFormatter
#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/

#***********************************************************************************************************
###########################################DIRECT EFFECT+ SEMI-DIRECT+ INDIRECT#############################
#***********************************************************************************************************

#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset(r"K:\WRF_Chem_Output\202\April\only_emiss_IGP\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")

time = ALL_TIMES

C = getvar(wrf_file, "RAINC", timeidx=1)
#######################################ALL SKY RADIATIVE FORCING AT SURFACE##################################
CLD = getvar(wrf_file, "LWDNB", timeidx=time)
CLD1 = getvar(wrf_file, "LWUPB", timeidx=time)
CLD2 = getvar(wrf_file1, "LWDNB", timeidx=time)
CLD3 = getvar(wrf_file1, "LWUPB", timeidx=time)

CLD_1 = getvar(wrf_file, "SWDNB", timeidx=time)
CLD1_1 = getvar(wrf_file, "SWUPB", timeidx=time)
CLD2_1 = getvar(wrf_file1, "SWDNB", timeidx=time)
CLD3_1 = getvar(wrf_file1, "SWUPB", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")

        CLD_1=CLD_1.mean("Time")
        CLD1_1 = CLD1_1.mean("Time")
        CLD2_1 = CLD2_1.mean("Time")
        CLD3_1 = CLD3_1.mean("Time")

RF_SURFACE_ALLSKY = ((CLD+CLD_1)-(CLD1+CLD1_1))-((CLD2+CLD2_1)-(CLD3+CLD3_1))
print("DIRECT EFFECT+ SEMI-DIRECT+ INDIRECT ALL SKY RF SURFACE")
precp1 = RF_SURFACE_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))
#######################################ALL SKY RADIATIVE FORCING AT TOA##################################
CLD = getvar(wrf_file, "LWDNT", timeidx=time)
CLD1 = getvar(wrf_file, "LWUPT", timeidx=time)
CLD2 = getvar(wrf_file1, "LWDNT", timeidx=time)
CLD3 = getvar(wrf_file1, "LWUPT", timeidx=time)

CLD_1 = getvar(wrf_file, "SWDNT", timeidx=time)
CLD1_1 = getvar(wrf_file, "SWUPT", timeidx=time)
CLD2_1 = getvar(wrf_file1, "SWDNT", timeidx=time)
CLD3_1 = getvar(wrf_file1, "SWUPT", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")

        CLD_1 = CLD_1.mean("Time")
        CLD1_1 = CLD1_1.mean("Time")
        CLD2_1 = CLD2_1.mean("Time")
        CLD3_1 = CLD3_1.mean("Time")

RF_TOA_ALLSKY = ((CLD + CLD_1) - (CLD1 + CLD1_1)) - ((CLD2 + CLD2_1) - (CLD3 + CLD3_1))
print("DIRECT EFFECT+ SEMI-DIRECT+ INDIRECT ALL SKY RF TOA")
precp1 = RF_TOA_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))
#######################################ALL SKY ATMOSPHERIC RADIATIVE FORCING##################################
RF_ATM_ALLSKY = RF_TOA_ALLSKY-RF_SURFACE_ALLSKY
print("DIRECT EFFECT+ SEMI-DIRECT+ INDIRECT ALL SKY RF ATM")
precp1 = RF_ATM_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))





lats,lons =w.latlon_coords(RF_ATM_ALLSKY)

font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
cart_proj = get_cartopy(C)
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

#ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())

ax.set_xticks([75, 80, 85, 90, 95,100], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())

#ax.set_xlim(cartopy_xlim(RAINNC766_1))
#ax.set_ylim(cartopy_ylim(RAINNC766_1))
#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),20) ##COLORBAR LIMITS
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF), 30)
v = np.linspace(-30.0,30.0,21)


plt.contourf(lons, lats, RF_SURFACE_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 6\radiative forcing\TOTAL RADIATIVE FORCING_SURFACE_ALLSKY_only_emiss_IGP.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_TOA_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 6\radiative forcing\TOTAL RADIATIVE FORCING_TOA_ALLSKY_only_emiss_IGP.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_ATM_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 6\radiative forcing\TOTAL RADIATIVE FORCING_ATM_ALLSKY_only_emiss_IGP.png",dpi=600,bbox_inches='tight')



#plt.show()
CLD=None
CLD1=None
CLD2=None
CLD3=None
CLD_1=None
CLD1_1=None
CLD2_1=None
CLD3_1=None

'''
#***********************************************************************************************************
###########################################INDIRECT EFFECT##################################################
#***********************************************************************************************************


wrf_file = Dataset(r"G:\WRF_Chem_Output\202\\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")

time = ALL_TIMES

C = getvar(wrf_file, "RAINC", timeidx=1)
#######################################ALL SKY RADIATIVE FORCING AT SURFACE##################################
CLD = getvar(wrf_file, "LWDNB", timeidx=time)
CLD1 = getvar(wrf_file, "LWUPB", timeidx=time)
CLD2 = getvar(wrf_file1, "LWDNB", timeidx=time)
CLD3 = getvar(wrf_file1, "LWUPB", timeidx=time)

CLD_1 = getvar(wrf_file, "SWDNB", timeidx=time)
CLD1_1 = getvar(wrf_file, "SWUPB", timeidx=time)
CLD2_1 = getvar(wrf_file1, "SWDNB", timeidx=time)
CLD3_1 = getvar(wrf_file1, "SWUPB", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")

        CLD_1=CLD_1.mean("Time")
        CLD1_1 = CLD1_1.mean("Time")
        CLD2_1 = CLD2_1.mean("Time")
        CLD3_1 = CLD3_1.mean("Time")

RF_SURFACE_ALLSKY = ((CLD+CLD_1)-(CLD1+CLD1_1))-((CLD2+CLD2_1)-(CLD3+CLD3_1))
print("INDIRECT ALL SKY RF SURFACE")
precp1 = RF_SURFACE_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))
#######################################ALL SKY RADIATIVE FORCING AT TOA##################################
CLD = getvar(wrf_file, "LWDNT", timeidx=time)
CLD1 = getvar(wrf_file, "LWUPT", timeidx=time)
CLD2 = getvar(wrf_file1, "LWDNT", timeidx=time)
CLD3 = getvar(wrf_file1, "LWUPT", timeidx=time)

CLD_1 = getvar(wrf_file, "SWDNT", timeidx=time)
CLD1_1 = getvar(wrf_file, "SWUPT", timeidx=time)
CLD2_1 = getvar(wrf_file1, "SWDNT", timeidx=time)
CLD3_1 = getvar(wrf_file1, "SWUPT", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")

        CLD_1 = CLD_1.mean("Time")
        CLD1_1 = CLD1_1.mean("Time")
        CLD2_1 = CLD2_1.mean("Time")
        CLD3_1 = CLD3_1.mean("Time")

RF_TOA_ALLSKY = ((CLD + CLD_1) - (CLD1 + CLD1_1)) - ((CLD2 + CLD2_1) - (CLD3 + CLD3_1))
print("INDIRECT ALL SKY RF TOA")
precp1 = RF_TOA_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))
#######################################ALL SKY ATMOSPHERIC RADIATIVE FORCING##################################
RF_ATM_ALLSKY = RF_TOA_ALLSKY-RF_SURFACE_ALLSKY
print("INDIRECT ALL SKY RF ATM")
precp1 = RF_ATM_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))



lats,lons =w.latlon_coords(RF_ATM_ALLSKY)

font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
cart_proj = get_cartopy(C)
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

#ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())

ax.set_xticks([75, 80, 85, 90, 95,100], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())

#ax.set_xlim(cartopy_xlim(RAINNC766_1))
#ax.set_ylim(cartopy_ylim(RAINNC766_1))
#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),20) ##COLORBAR LIMITS
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF), 30)
v = np.linspace(-30.0,30.0,21)

plt.contourf(lons, lats, RF_SURFACE_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\INDIRECT RADIATIVE FORCING_SURFACE_ALLSKY.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_TOA_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\INDIRECT RADIATIVE FORCING_TOA_ALLSKY.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_ATM_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\INDIRECT RADIATIVE FORCING_ATM_ALLSKY.png",dpi=600,bbox_inches='tight')



#plt.show()
CLD=None
CLD1=None
CLD2=None
CLD3=None
CLD_1=None
CLD1_1=None
CLD2_1=None
CLD3_1=None

#***********************************************************************************************************
###########################################DIRECT +SEMI-DIRECT EFFECT##################################################
#***********************************************************************************************************
wrf_file = Dataset(r"G:\WRF_Chem_Output\202\\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")

time = ALL_TIMES

C = getvar(wrf_file, "RAINC", timeidx=1)
#######################################ALL SKY RADIATIVE FORCING AT SURFACE##################################
CLD = getvar(wrf_file, "LWDNB", timeidx=time)
CLD1 = getvar(wrf_file, "LWUPB", timeidx=time)
CLD2 = getvar(wrf_file1, "LWDNB", timeidx=time)
CLD3 = getvar(wrf_file1, "LWUPB", timeidx=time)

CLD_1 = getvar(wrf_file, "SWDNB", timeidx=time)
CLD1_1 = getvar(wrf_file, "SWUPB", timeidx=time)
CLD2_1 = getvar(wrf_file1, "SWDNB", timeidx=time)
CLD3_1 = getvar(wrf_file1, "SWUPB", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")

        CLD_1=CLD_1.mean("Time")
        CLD1_1 = CLD1_1.mean("Time")
        CLD2_1 = CLD2_1.mean("Time")
        CLD3_1 = CLD3_1.mean("Time")

RF_SURFACE_ALLSKY = ((CLD+CLD_1)-(CLD1+CLD1_1))-((CLD2+CLD2_1)-(CLD3+CLD3_1))
print("DIRECT +SEMI-DIRECT ALL SKY RF SURFACE")
precp1 = RF_SURFACE_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))
#######################################ALL SKY RADIATIVE FORCING AT TOA##################################
CLD = getvar(wrf_file, "LWDNT", timeidx=time)
CLD1 = getvar(wrf_file, "LWUPT", timeidx=time)
CLD2 = getvar(wrf_file1, "LWDNT", timeidx=time)
CLD3 = getvar(wrf_file1, "LWUPT", timeidx=time)

CLD_1 = getvar(wrf_file, "SWDNT", timeidx=time)
CLD1_1 = getvar(wrf_file, "SWUPT", timeidx=time)
CLD2_1 = getvar(wrf_file1, "SWDNT", timeidx=time)
CLD3_1 = getvar(wrf_file1, "SWUPT", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")

        CLD_1 = CLD_1.mean("Time")
        CLD1_1 = CLD1_1.mean("Time")
        CLD2_1 = CLD2_1.mean("Time")
        CLD3_1 = CLD3_1.mean("Time")

RF_TOA_ALLSKY = ((CLD + CLD_1) - (CLD1 + CLD1_1)) - ((CLD2 + CLD2_1) - (CLD3 + CLD3_1))
print("DIRECT +SEMI-DIRECT ALL SKY RF TOA")
precp1 = RF_TOA_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_TOA_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))
#######################################ALL SKY ATMOSPHERIC RADIATIVE FORCING##################################
RF_ATM_ALLSKY = RF_TOA_ALLSKY-RF_SURFACE_ALLSKY
print("DIRECT +SEMI-DIRECT ALL SKY RF ATM")
precp1 = RF_ATM_ALLSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_ATM_ALLSKY[168:178, 174:216]
print(np.nanmean(precp1))

#***********************************************************************************************************
###########################################DIRECT EFFECT##################################################
#***********************************************************************************************************
wrf_file = Dataset(r"G:\WRF_Chem_Output\202\\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_WRF\wrfout_d01_2018-04-10_00%3A00%3A00")

time = ALL_TIMES

C = getvar(wrf_file, "RAINC", timeidx=1)
#######################################CLEAR SKY RADIATIVE FORCING AT SURFACE##################################
CLD = getvar(wrf_file, "LWDNBC", timeidx=time)
CLD1 = getvar(wrf_file, "LWUPBC", timeidx=time)
CLD2 = getvar(wrf_file1, "LWDNBC", timeidx=time)
CLD3 = getvar(wrf_file1, "LWUPBC", timeidx=time)
CLD_1 = getvar(wrf_file, "SWDNBC", timeidx=time)
CLD1_1 = getvar(wrf_file, "SWUPBC", timeidx=time)
CLD2_1 = getvar(wrf_file1, "SWDNBC", timeidx=time)
CLD3_1 = getvar(wrf_file1, "SWUPBC", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")

        CLD_1 =CLD_1.mean("Time")
        CLD1_1 = CLD1_1.mean("Time")
        CLD2_1 = CLD2_1.mean("Time")
        CLD3_1 = CLD3_1.mean("Time")

RF_SURFACE_CLEARSKY = ((CLD+CLD_1)-(CLD1+CLD1_1))-((CLD2+CLD2_1)-(CLD3+CLD3_1))
print("DIRECT CLEAR SKY RF SURFACE")
precp1 = RF_SURFACE_CLEARSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_CLEARSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_CLEARSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_CLEARSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_SURFACE_CLEARSKY[168:178, 174:216]
print(np.nanmean(precp1))
#######################################CLEAR SKY RADIATIVE FORCING AT TOA##################################
CLD = getvar(wrf_file, "LWDNTC", timeidx=time)
CLD1 = getvar(wrf_file, "LWUPTC", timeidx=time)
CLD2 = getvar(wrf_file1, "LWDNTC", timeidx=time)
CLD3 = getvar(wrf_file1, "LWUPTC", timeidx=time)

CLD_1 = getvar(wrf_file, "SWDNTC", timeidx=time)
CLD1_1 = getvar(wrf_file, "SWUPTC", timeidx=time)
CLD2_1 = getvar(wrf_file1, "SWDNTC", timeidx=time)
CLD3_1 = getvar(wrf_file1, "SWUPTC", timeidx=time)

if time == ALL_TIMES:
        CLD = CLD.mean("Time")
        CLD1 = CLD1.mean("Time")
        CLD2 = CLD2.mean("Time")
        CLD3 = CLD3.mean("Time")

        CLD_1=CLD_1.mean("Time")
        CLD1_1 = CLD1_1.mean("Time")
        CLD2_1 = CLD2_1.mean("Time")
        CLD3_1 = CLD3_1.mean("Time")

RF_TOA_CLEARSKY = ((CLD+CLD_1)-(CLD1+CLD1_1))-((CLD2+CLD2_1)-(CLD3+CLD3_1))
print("DIRECT CLEAR SKY RF TOA")
precp1 = RF_TOA_CLEARSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_TOA_CLEARSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_TOA_CLEARSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_TOA_CLEARSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_TOA_CLEARSKY[168:178, 174:216]
print(np.nanmean(precp1))
#######################################CLEAR SKY ATMOSPHERIC RADIATIVE FORCING##################################
RF_ATM_CLEARSKY = RF_TOA_CLEARSKY-RF_SURFACE_CLEARSKY
print("DIRECT CLEAR SKY RF ATM")
precp1 = RF_ATM_CLEARSKY[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_ATM_CLEARSKY[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_ATM_CLEARSKY[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_ATM_CLEARSKY[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_ATM_CLEARSKY[168:178, 174:216]
print(np.nanmean(precp1))

lats,lons =w.latlon_coords(RF_ATM_ALLSKY)

font = {'family': 'serif','color':  'black','weight': 'bold','size': 12,}
cart_proj = get_cartopy(C)
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

#ax.set_xticks([88, 90, 92, 94, 96, 98], crs=ccrs.PlateCarree())
#ax.xaxis.set_major_formatter(LongitudeFormatter())
#ax.set_yticks([22, 24, 26, 28, 30], crs=ccrs.PlateCarree())
#ax.yaxis.set_major_formatter(LatitudeFormatter())

ax.set_xticks([75, 80, 85, 90, 95,100], crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.set_yticks([15, 20, 25, 30], crs=ccrs.PlateCarree())
ax.yaxis.set_major_formatter(LatitudeFormatter())

#ax.set_xlim(cartopy_xlim(RAINNC766_1))
#ax.set_ylim(cartopy_ylim(RAINNC766_1))
#ax.gridlines(color="black", linestyle="dotted")
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF),20) ##COLORBAR LIMITS
#v = np.linspace(np.min(CLD_DIFF),np.max(CLD_DIFF), 30)
v = np.linspace(-30.0,30.0,21)


plt.contourf(lons, lats, RF_SURFACE_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\DIRECT+SEMI RADIATIVE FORCING_SURFACE_ALLSKY.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_TOA_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\DIRECT+SEMI RADIATIVE FORCING_TOA_ALLSKY.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_ATM_ALLSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\DIRECT+SEMI RADIATIVE FORCING_ATM_ALLSKY.png",dpi=600,bbox_inches='tight')


plt.contourf(lons, lats, RF_SURFACE_CLEARSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\DIRECT RADIATIVE FORCING_SURFACE_CLEARSKY.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_TOA_CLEARSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\DIRECT RADIATIVE FORCING_TOA_CLEARSKY.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_ATM_CLEARSKY,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\DIRECT RADIATIVE FORCING__ATM_CLEARSKY.png",dpi=600,bbox_inches='tight')

#############################SEMI DIRECT EFFECT###############################
RF_SEMI_SURFACE=RF_SURFACE_ALLSKY-RF_SURFACE_CLEARSKY
print("SEMI DIRECT RF NET SURFACE")
precp1 = RF_SEMI_SURFACE[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_SEMI_SURFACE[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SEMI_SURFACE[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SEMI_SURFACE[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_SEMI_SURFACE[168:178, 174:216]
print(np.nanmean(precp1))

RF_SEMI_TOA=RF_TOA_ALLSKY-RF_TOA_CLEARSKY
print("SEMI DIRECT RF NET TOA")
precp1 = RF_SEMI_TOA[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_SEMI_TOA[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SEMI_TOA[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SEMI_TOA[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_SEMI_TOA[168:178, 174:216]
print(np.nanmean(precp1))

RF_SEMI_ATM=RF_ATM_ALLSKY-RF_ATM_CLEARSKY
print("SEMI DIRECT RF NET ATM")
precp1 = RF_SEMI_ATM[157:186, 174:216]
print(np.nanmean(precp1))
precp1 = RF_SEMI_ATM[157:186, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SEMI_ATM[186:203, 216:252]
print(np.nanmean(precp1))
precp1 = RF_SEMI_ATM[123:157, 200:237]
print(np.nanmean(precp1))
precp1 = RF_SEMI_ATM[168:178, 174:216]
print(np.nanmean(precp1))

plt.contourf(lons, lats, RF_SEMI_SURFACE,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\SEMI-DIRECT NET RADIATIVE FORCING_SURFACE.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_SEMI_TOA,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\SEMI-DIRECT NET RADIATIVE FORCING_TOA.png",dpi=600,bbox_inches='tight')

plt.contourf(lons, lats, RF_SEMI_ATM,v, cmap=plt.get_cmap("bwr"),antialiased=False,transform=cartopy.crs.PlateCarree(),extend='both')
#plt.title("TSK difference",weight='bold',fontdict=font)
#plt.figure(figsize=(8,6))
#plt.colorbar(ax=ax, shrink=.81, ticks=v)
plt.savefig(r"F:\WRF-CHEM ANALYSIS Chap 5\radiative forcing\SEMI-DIRECT NET RADIATIVE FORCING_ATM.png",dpi=600,bbox_inches='tight')


#plt.show()
CLD=None
CLD1=None
CLD2=None
CLD3=None
CLD_1=None
CLD1_1=None
CLD2_1=None
CLD3_1=None
'''
