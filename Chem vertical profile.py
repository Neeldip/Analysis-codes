from xarray import open_dataset
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import os
#filename = "name.png"
#path = "/path/to/save/location"
#fullpath = os.path.join(path, filename)

#DEFINE DATASETS
dsFeb = open_dataset("F:\WRF-CHEM ANALYSIS\Feb_fraction_conc.nc")
dsMar = open_dataset("F:\WRF-CHEM ANALYSIS\March_fraction_conc.nc")
dsApr = open_dataset("F:\WRF-CHEM ANALYSIS\April_fraction_conc.nc")
dsMay = open_dataset("F:\WRF-CHEM ANALYSIS\May_fraction_conc.nc")
dsJune = open_dataset("F:\WRF-CHEM ANALYSIS\June_fraction_conc.nc")
dsJuly = open_dataset("F:\WRF-CHEM ANALYSIS\July_fraction_conc.nc")
dsAug = open_dataset("F:\WRF-CHEM ANALYSIS\Aug_fraction_conc.nc")
dsSep = open_dataset("F:\WRF-CHEM ANALYSIS\Sep_fraction_conc.nc")

lat = dsFeb.variables["lat"][:]
lon = dsFeb.variables["lon"][:]

#SPECIFY LOCATION OF PROFILE
#For region 1
#lats=[26.32,26.40,26.47,26.18,26.31,26.44,26.17,25.28,25.56,26.02]
#lons=[89.45,90.27,90.54,90.62,91.00,91.44,91.75,91.56,91.87,89.97]

#For region 2
#lats=[27.22,27.48,27.48,27.28,27.92,28.06,27.47,27.58,27.66]
#lons=[94.09,94.92,95.36,95.65,96.14,95.32,94.55,91.86,95.86]

#For region 3
#lats=[26.34,26.62,25.91,26.76,25.65,26.98]
#lons=[92.68,92.78,93.72,94.21,94.09,94.63]

#For region 4
#lats=[24.80,23.83,23.74,24.81,22.91]
#lons=[93.93,91.28,92.73,92.75,92.74]

lats = [15]
lons = [77.5]
#print(i)
#print(j)


#BC_prop = dsApr.variables.keys()
#print(BC_prop)
for sel_lat,sel_lon in zip(lats,lons):
    a = abs(lat - sel_lat) + abs(lon - sel_lon)
    i, j = np.unravel_index(a.argmin(), a.shape)
    BC_conc = dsApr.variables["BC_conc"][:,i,j]
    SO4_conc = dsApr.variables["SO4_conc"][:,i,j]
    OC_conc = dsApr.variables["OC_conc"][:,i,j]
    OIN_conc = dsApr.variables["OIN_conc"][:,i,j]
    NO3_conc = dsApr.variables["NO3_conc"][:,i,j]
    NH4_conc = dsApr.variables["NH4_conc"][:,i,j]
    level = np.arange(0,44,1)
    plt.plot(BC_conc,level,label="BC_conc")
    plt.plot(SO4_conc,level,label="SO4_conc")
    plt.plot(OC_conc,level,label="OC_conc")
    plt.plot(OIN_conc,level,label="OIN_conc")
    plt.plot(NO3_conc,level,label="NO3_conc")
    plt.plot(NH4_conc,level,label="NH4_conc")
    plt.legend()
    plt.title("Aerosol_concentration at "+str(sel_lat)+"_"+str(sel_lon))
    #plt.savefig("Aerosol_concentration in Sep (REGION1) " +str(sel_lat)+"_"+str(sel_lon)+".png",bbox_inches='tight',dpi=300)
    plt.show()

'''
for sel_lat, sel_lon in zip(lats, lons):
    b = abs(lat - sel_lat) + abs(lon - sel_lon)
    i, j = np.unravel_index(b.argmin(), b.shape)
    BC_prop = dsApr.variables["BC_FRACTION"][:,i,j]
    SO4_prop = dsApr.variables["SO4_FRACTION"][:,i,j]
    OC_prop = dsApr.variables["OC_FRACTION"][:,i,j]
    OIN_prop = dsApr.variables["OIN_FRACTION"][:,i,j]
    NO3_prop = dsApr.variables["NO3_FRACTION"][:,i,j]
    NH4_prop = dsApr.variables["NH4_FRACTION"][:,i,j]
    level = np.arange(0,44,1)
    plt.plot(BC_prop,level,label="BC_fraction")
    plt.plot(SO4_prop,level,label="SO4_fraction")
    plt.plot(OC_prop,level,label="OC_fraction")
    plt.plot(OIN_prop,level,label="OIN_fraction")
    plt.plot(NO3_prop,level,label="NO3_fraction")
    plt.plot(NH4_prop,level,label="NH4_fraction")
    plt.legend()
    plt.title("Aerosol_fraction at "+str(sel_lat)+"_"+str(sel_lon))
    plt.savefig("Aerosol_fraction in April at "+str(sel_lat)+"_"+str(sel_lon)+".png",bbox_inches='tight',dpi=300)
    #plt.show()
'''

