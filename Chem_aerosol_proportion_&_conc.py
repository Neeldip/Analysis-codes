import numpy as np
import netCDF4 as nc
from wrf import getvar,omp_set_num_threads,ALL_TIMES
omp_set_num_threads(4)

wrf =  nc.Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf =  [nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-10_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-11_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-12_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-13_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-14_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-15_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-16_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-17_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-18_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-19_00%3A00%3A00"),
#        nc.Dataset(r"H:\WRF_Chem_Output\201\March\wrfout_d01_2018-03-20_00%3A00%3A00")]

#'''
wrf1 = nc.Dataset(r"H:\WRF_Chem_Output\201\July\wrfout_d01_2018-07-10_00%3A00%3A00")

wrf_lat = getvar(wrf,"XLAT",timeidx=0)
wrf_lon = getvar(wrf,"XLONG",timeidx=0)

#####select region########
sel_lat = 10.64#21
sel_lon = 71.67#83

a = abs(wrf_lat - sel_lat) + abs(wrf_lon - sel_lon)
b, c = np.unravel_index(a.argmin(), a.shape)
#print(b)
#print(c)
sel_lat = 31.22#30
sel_lon = 100.43#97

d = abs(wrf_lat - sel_lat) + abs(wrf_lon - sel_lon)
e, f = np.unravel_index(d.argmin(), d.shape)
#print(e)
#print(f)

####specify model levels to process#########
model_levels = 44

wrf_lat = wrf1.variables["XLAT"][0,b:e,0]
wrf_lon = wrf1.variables["XLONG"][0,0,c:f]
#wrf_lat = wrf.variables["XLAT"][0,:,0]
#wrf_lon = wrf.variables["XLONG"][0,0,:]
#print(np.shape(wrf_lat))
#print(np.shape(wrf_lon))

############create netcdf################

fn = 'F:\WRF-CHEM ANALYSIS\April_BC_fraction_conc_MYNN3.nc'   #<----------------
ds = nc.Dataset(fn, 'w', format='NETCDF4')
ds.set_fill_on()

#lat = ds.createDimension('lat', 35)
#lon = ds.createDimension('lon', 44)
lat = ds.createDimension('lat', e-b)
lon = ds.createDimension('lon', f-c)
lev = ds.createDimension('lev', model_levels)
lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
level = ds.createVariable('lev','f4', ('lev',))
#BC_FRACTION = ds.createVariable('BC_FRACTION', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
#BC_FRACTION.units = 'BC_FRACTION(%) out of so4 no3 bc oc oin nh4'
'''
SO4_FRACTION = ds.createVariable('SO4_FRACTION', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
SO4_FRACTION.units = 'SO4_FRACTION(%) out of so4 no3 bc oc oin nh4'

NO3_FRACTION = ds.createVariable('NO3_FRACTION', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
NO3_FRACTION.units = 'NO3_FRACTION(%) out of so4 no3 bc oc oin nh4'

OC_FRACTION = ds.createVariable('OC_FRACTION', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
OC_FRACTION.units = 'OC_FRACTION(%) out of so4 no3 bc oc oin nh4'

OIN_FRACTION = ds.createVariable('OIN_FRACTION', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
OIN_FRACTION.units = 'OIN_FRACTION(%) out of so4 no3 bc oc oin nh4'

NH4_FRACTION = ds.createVariable('NH4_FRACTION', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
NH4_FRACTION.units = 'NH4_FRACTION(%) out of so4 no3 bc oc oin nh4'
'''
BC_conc = ds.createVariable('BC_conc', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
BC_conc.units = 'BC_conc ug/kg dry air'
'''
SO4_conc = ds.createVariable('SO4_conc', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
SO4_conc.units = 'SO4_conc ug/kg dry air'

NO3_conc = ds.createVariable('NO3_conc', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
NO3_conc.units = 'NO3_conc ug/kg dry air'

OC_conc = ds.createVariable('OC_conc', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
OC_conc.units = 'OC_conc ug/kg dry air'

OIN_conc = ds.createVariable('OIN_conc', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
OIN_conc.units = 'OIN_conc ug/kg dry air'

NH4_conc = ds.createVariable('NH4_conc', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
NH4_conc.units = 'NH4_conc ug/kg dry air'
'''
#WATER_FRACTION = ds.createVariable('WATER_FRACTION', 'f4', ('lev', 'lat','lon'),fill_value=np.nan)
#WATER_FRACTION.units = 'WATER_FRACTION(%) out of so4 no3 bc oc oin nh4 water'

lats.units = 'degrees north'
lons.units = 'degrees east'
level.units = 'model level number'

level[:] = np.arange(0,model_levels,1)
lons[:] = wrf_lon#np.arange(71.67989, 100.52737994, 0.09647994)
lats[:] = wrf_lat# np.arange(10.64795, 31.32473994,  0.08950991 )

#bc_a01 = wrf1.variables["bc_a01"][:, 0:model_levels, b:e, c:f]
#bc_a02 = wrf1.variables["bc_a02"][:, 0:model_levels, b:e, c:f]
#bc_a03 = wrf1.variables["bc_a03"][:, 0:model_levels, b:e, c:f]
#bc_a04 = wrf1.variables["bc_a04"][:, 0:model_levels, b:e, c:f]
#pm10 = wrf1.variables["PM10"][:, 0:model_levels, b:e, c:f]
#T = trad1.variables["TEMPERATURE"][:, 0:model_levels, b:e, c:f]
#P = trad1.variables["PRESSURE"][:, 0:model_levels, b:e, c:f]

bc_a01 = getvar(wrf,"bc_a01",timeidx=ALL_TIMES)
bc_a02 = getvar(wrf,"bc_a02",timeidx=ALL_TIMES)
bc_a03 = getvar(wrf,"bc_a03",timeidx=ALL_TIMES)
bc_a04 = getvar(wrf,"bc_a04",timeidx=ALL_TIMES)
BC = bc_a01+bc_a02+bc_a03+bc_a04
bc_a01 = None
bc_a02 = None
bc_a03 = None
bc_a04 = None
bc1 = BC.mean(axis=0)
BC = None
'''
so4_a01 = getvar(wrf,"so4_a01",timeidx=ALL_TIMES)
so4_a02 = getvar(wrf,"so4_a02",timeidx=ALL_TIMES)
so4_a03 = getvar(wrf,"so4_a03",timeidx=ALL_TIMES)
so4_a04 = getvar(wrf,"so4_a04",timeidx=ALL_TIMES)
so4 = so4_a01+so4_a02+so4_a03+so4_a04
so4_a01 = None
so4_a02 = None
so4_a03 = None
so4_a04 = None
so41 = so4.mean(axis=0)
so4 = None

no3_a01 = getvar(wrf,"no3_a01",timeidx=ALL_TIMES)
no3_a02 = getvar(wrf,"no3_a02",timeidx=ALL_TIMES)
no3_a03 = getvar(wrf,"no3_a03",timeidx=ALL_TIMES)
no3_a04 = getvar(wrf,"no3_a04",timeidx=ALL_TIMES)
no3 = no3_a01+no3_a02+no3_a03+no3_a04
no3_a01 = None
no3_a02 = None
no3_a03 = None
no3_a04 = None
no31 = no3.mean(axis=0)
no3 = None

oc_a01 = getvar(wrf,"oc_a01",timeidx=ALL_TIMES)
oc_a02 = getvar(wrf,"oc_a02",timeidx=ALL_TIMES)
oc_a03 = getvar(wrf,"oc_a03",timeidx=ALL_TIMES)
oc_a04 = getvar(wrf,"oc_a04",timeidx=ALL_TIMES)
oc = oc_a01+oc_a02+oc_a03+oc_a04
oc_a01 = None
oc_a02 = None
oc_a03 = None
oc_a04 = None
oc1 = oc.mean(axis=0)
oc = None

oin_a01 = getvar(wrf,"oin_a01",timeidx=ALL_TIMES)
oin_a02 = getvar(wrf,"oin_a02",timeidx=ALL_TIMES)
oin_a03 = getvar(wrf,"oin_a03",timeidx=ALL_TIMES)
oin_a04 = getvar(wrf,"oin_a04",timeidx=ALL_TIMES)
oin = oin_a01+oin_a02+oin_a03+oin_a04
oin_a01 = None
oin_a02 = None
oin_a03 = None
oin_a04 = None
oin1 = oin.mean(axis=0)
oin = None

nh4_a01 = getvar(wrf,"nh4_a01",timeidx=ALL_TIMES)
nh4_a02 = getvar(wrf,"nh4_a02",timeidx=ALL_TIMES)
nh4_a03 = getvar(wrf,"nh4_a03",timeidx=ALL_TIMES)
nh4_a04 = getvar(wrf,"nh4_a04",timeidx=ALL_TIMES)
nh4 = nh4_a01+nh4_a02+nh4_a03+nh4_a04
nh4_a01 = None
nh4_a02 = None
nh4_a03 = None
nh4_a04 = None
nh41 = nh4.mean(axis=0)
nh4 = None

#water_a01 = getvar(wrf,"water_a01",timeidx=ALL_TIMES)
#water_a02 = getvar(wrf,"water_a02",timeidx=ALL_TIMES)
#water_a03 = getvar(wrf,"water_a03",timeidx=ALL_TIMES)
#water_a04 = getvar(wrf,"water_a04",timeidx=ALL_TIMES)
#water = water_a01+water_a02+water_a03+water_a04
#water_a01 = None
#water_a02 = None
#water_a03 = None
#water_a04 = None
#water1 = water.mean(axis=0)
#water = None
#pm10 = getvar(wrf,"PM10",timeidx=ALL_TIMES)
#T = getvar(trad,"PRESSURE",timeidx=ALL_TIMES)
#P = getvar(trad,"TEMPERATURE",timeidx=ALL_TIMES)

BC_frac = (bc1/(bc1+so41+no31+oc1+oin1+nh41))*100
#SO4_frac = (so41/(bc1+so41+no31+oc1+oin1+nh41))*100
#NO3_frac = (no31/(bc1+so41+no31+oc1+oin1+nh41))*100
#OC_frac = (oc1/(bc1+so41+no31+oc1+oin1+nh41))*100
#OIN_frac = (oin1/(bc1+so41+no31+oc1+oin1+nh41))*100
#NH4_frac = (nh41/(bc1+so41+no31+oc1+oin1+nh41))*100
#WATER_frac = (water1/(bc1+so41+no31+oc1+oin1+nh41+water1))*100
'''
print("mean complete")
#########################################PROPORTION########################################
###########ADDS NETCDF VALUES#############
for i in range(0,model_levels,1):
    for j in range(0,e-b,1):
        for k in range(0,f-c,1):
            #BC_FRACTION[i,j,k] = BC_frac[i,j+b,k+c]
            #SO4_FRACTION[i, j, k] = SO4_frac[i, j + b, k + c]
            #NO3_FRACTION[i, j, k] = NO3_frac[i, j + b, k + c]
            #OC_FRACTION[i, j, k] = OC_frac[i, j + b, k + c]
            #OIN_FRACTION[i, j, k] = OIN_frac[i, j + b, k + c]
            #NH4_FRACTION[i, j, k] = NH4_frac[i, j + b, k + c]
            BC_conc[i, j, k] = bc1[i, j + b, k + c]
            #SO4_conc[i, j, k] = so41[i, j + b, k + c]
            #NO3_conc[i, j, k] = no31[i, j + b, k + c]
            #OC_conc[i, j, k] = oc1[i, j + b, k + c]
            #OIN_conc[i, j, k] = oin1[i, j + b, k + c]
            #NH4_conc[i, j, k] = nh41[i, j + b, k + c]
            #WATER_FRACTION[i, j, k] = WATER_frac[i, j + b, k + c]

ds.close()

BC_FRACTION=None
SO4_FRACTION=None
NH4_FRACTION=None
NO3_FRACTION=None
OIN_FRACTION=None
OC_FRACTION=None
bc1=None
so41=None
no31=None
oc1=None
oin1=None
nh41=None
