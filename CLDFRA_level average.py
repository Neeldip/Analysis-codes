import wrf as w
from wrf import latlon_coords,omp_set_num_threads,get_cartopy,getvar,ALL_TIMES,cartopy_ylim,cartopy_xlim
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from cartopy.feature import NaturalEarthFeature
import cartopy as cartopy
import cartopy.crs as ccrs
#from xarray import Dataset


#https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
#trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"G:\WRF_Chem_Output\202\NOR\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\NO_BC_ABS\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_2xBC\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_4xBC\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_BC_no_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
#wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_dust_bc_absorbtion\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file1 = Dataset(r"G:\WRF_Chem_Output\202\NOFEED\wrfout_d01_2018-04-10_00%3A00%3A00")
time = ALL_TIMES

#time = 60
#level = range(13, 18, 1)  # [5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,30]
level = range(0,44,1)
C = getvar(wrf_file, "V10", timeidx=1)
CLD = getvar(wrf_file, "QVAPOR", timeidx=time)
CLD1 = getvar(wrf_file1, "QVAPOR", timeidx=time)
#print(CLD1)

if time == ALL_TIMES:
    CLD = CLD.mean(axis=0)
    CLD1 = CLD1.mean(axis=0)

    # CLD = CLD.std("Time")
    # CLD1 = CLD1.std("Time")

'''
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z) #* 100  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    ne=CLD_DIFF[123:204, 182:263]
    #ne_nosouth = CLD_DIFF[157:204, 182:263]
    print(np.nanmean(ne))
    #print(np.nanmean(ne_nosouth))
'''

print("######")

#for y in level:
#    # lats, lons = w.latlon_coords(CLD)
#    x = CLD[y, :, :]
#    z = CLD1[y, :, :]
#    CLD_DIFF = (x - z) * 1000  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
#    #ne=CLD_DIFF[123:204, 182:263]
#    ne_nosouth = CLD_DIFF[157:204, 182:263]
#    #print(np.nanmean(ne))
#    print(np.nanmean(ne_nosouth))
'''
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z) #* 100  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    #ne=CLD_DIFF[123:204, 182:263]
    ne_brv = CLD_DIFF[157:180, 182:263]
    #print(np.nanmean(ne))
    print(np.nanmean(ne_brv))
print("######")
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z) #* 100  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    #ne=CLD_DIFF[123:204, 182:263]
    ne_arunachal = CLD_DIFF[181:204, 182:263]
    #print(np.nanmean(ne))
    print(np.nanmean(ne_arunachal))
'''
print('########')
#'''
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z)# * 1000  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    #ne=CLD_DIFF[123:204, 182:263]
    #ne_west = CLD_DIFF[157:180, 182:216]
    ne_west = CLD_DIFF[157:186, 174:216]
    print(np.nanmean(ne_west))
print('########')
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z)# * 1000  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    #ne=CLD_DIFF[123:204, 182:263]
    #ne_west = CLD_DIFF[157:180, 182:216]
    #ne_west = CLD_DIFF[157:186, 174:216]
    #print(np.nanmean(ne_west))
    ne_west = CLD_DIFF[157:186, 216:252]
    print(np.nanmean(ne_west))
    #ne_west = CLD_DIFF[186:203, 216:252]
    #print(np.nanmean(ne_west))
    #ne_west = CLD_DIFF[123:157, 200:237]
    #print(np.nanmean(ne_west))
print('########')
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z)# * 1000  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    #ne=CLD_DIFF[123:204, 182:263]
    #ne_west = CLD_DIFF[157:180, 182:216]
    #ne_west = CLD_DIFF[157:186, 174:216]
    #print(np.nanmean(ne_west))
    #ne_west = CLD_DIFF[157:186, 216:252]
    #print(np.nanmean(ne_west))
    ne_west = CLD_DIFF[186:203, 216:252]
    print(np.nanmean(ne_west))
    #ne_west = CLD_DIFF[123:157, 200:237]
    #print(np.nanmean(ne_west))
print('########')
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z)# * 1000  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    #ne=CLD_DIFF[123:204, 182:263]
    #ne_west = CLD_DIFF[157:180, 182:216]
    #ne_west = CLD_DIFF[157:186, 174:216]
    #print(np.nanmean(ne_west))
    #ne_west = CLD_DIFF[157:186, 216:252]
    #print(np.nanmean(ne_west))
    #ne_west = CLD_DIFF[186:203, 216:252]
    #print(np.nanmean(ne_west))
    ne_west = CLD_DIFF[123:157, 200:237]
    print(np.nanmean(ne_west))
print('########')
'''
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z)# * 1000  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    #ne=CLD_DIFF[123:204, 182:263]
    #ne_west = CLD_DIFF[157:180, 182:216]
    #ne_west = CLD_DIFF[181:182, 181:182]
    ne_west = CLD_DIFF[174:175, 197:198]
    print(np.nanmean(ne_west))
'''
for y in level:
    # lats, lons = w.latlon_coords(CLD)
    x = CLD[y, :, :]
    z = CLD1[y, :, :]
    CLD_DIFF = (x - z)# * 1000  # 0000#*1000000000 ###*100 FOR BETTER READABILITY
    #ne=CLD_DIFF[123:204, 182:263]
    #ne_west = CLD_DIFF[157:180, 182:216]
    ne_west = CLD_DIFF[168:178, 174:216]
    print(np.nanmean(ne_west))
#'''

CLD = None
CLD1 = None
