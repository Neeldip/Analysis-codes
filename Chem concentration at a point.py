#import glob
#list = glob.glob(r"I:\observation_data\APRIL\DATA\Rainfall\IMD_gridded_datasets\RAW data\NETCDF\*.nc")
from xarray import open_mfdataset
import numpy as np
from xarray import open_dataset
from wrf import ll_to_xy
#print(list)
import warnings
warnings.filterwarnings("ignore")
from netCDF4 import Dataset
#import sys
#np.set_printoptions(threshold=sys.maxsize)
#sys.stdout = open(r'F:\WRF-CHEM ANALYSIS\aerosols.txt', 'w')
lats=[26.48,28.62,23.18,16.51,26.68,23.54,23.68,22.54,26.83,17.71,13.16,22.96,26.90,28.64]
lons=[80.24,77.24,75.76,80.51,88.41,87.28,86.94,88.34,80.89,83.30,80.26,76.06,75.83,77.31]

level =[0]

for sel_lat,sel_lon in zip(lats,lons):
    ds = Dataset(r"H:\WRF_Chem_Output\201\April\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
    #ds = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
    #ds1 = open_dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
    ds1 = open_dataset(r"H:\WRF_Chem_Output\201\April\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
    wrf_lats = ds.variables['XLAT'][0, :, :]  ##231*299
    wrf_lons = ds.variables['XLONG'][0, :, :]
    a = abs(wrf_lats-sel_lat)+abs(wrf_lons-sel_lon)
    i,j = np.unravel_index(a.argmin(), a.shape)
    #i,j=ll_to_xy(ds,sel_lat,sel_lon)
    #print(i)
    #print(j)
    o3 = ds.variables["o3"][:, level, i, j] ###unit=ppmv
    no = ds.variables["no"][:, level, i, j] ###unit=ppmv
    no2 = ds.variables["no2"][:, level, i, j] ###unit=ppmv
    so2 = ds.variables["so2"][:, level, i, j] ###unit=ppmv
    co = ds.variables["co"][:, level, i, j] ###unit=ppmv
    nh3 = ds.variables["nh3"][:, level, i, j] ###unit=ppmv
    T = ds1.variables["TEMPERATURE"][:, level, i, j]
    PRESS = ds1.variables["PRESSURE"][:, level, i, j]
    K = (287*T)/PRESS

    C_o3 = (PRESS * o3 * 48) / (287 * T)
    C_no = (PRESS * no * 30 ) / (287 * T)
    C_no2 = (PRESS * no2 * 46) / (287 * T)
    C_so2 = (PRESS * so2 * 64.07) / (287 * T)
    C_co = (PRESS * co * 28 ) / (287 * T)
    C_nh3 = (PRESS * nh3 * 17) / (287 * T)

    pm25 = ds.variables["PM2_5_DRY"][:, level, i, j]  #unit = ug/m3
    pm10 = ds.variables["PM10"][:, level, i, j] #unit = ug/m3
    
    bc_a01 = ds.variables["bc_a01"][:, level, i, j]
    bc_a02 = ds.variables["bc_a02"][:, level, i, j]
    bc_a03 = ds.variables["bc_a03"][:, level, i, j]
    bc_a04 = ds.variables["bc_a04"][:, level, i, j]
    bc = bc_a01+bc_a02+bc_a03+bc_a04
    bc = bc/K

    oc_a01 = ds.variables["oc_a01"][:, level, i, j]
    oc_a02 = ds.variables["oc_a02"][:, level, i, j]
    oc_a03 = ds.variables["oc_a03"][:, level, i, j]
    oc_a04 = ds.variables["oc_a04"][:, level, i, j]
    oc = oc_a01+oc_a02+oc_a03+oc_a04
    oc = oc/K

    so4_a01 = ds.variables["so4_a01"][:, level, i, j]
    so4_a02 = ds.variables["so4_a02"][:, level, i, j]
    so4_a03 = ds.variables["so4_a03"][:, level, i, j]
    so4_a04 = ds.variables["so4_a04"][:, level, i, j]
    so4= so4_a01+so4_a02+so4_a03+so4_a04
    so4= so4/K

    oin_a01 = ds.variables["oin_a01"][:, level, i, j]
    oin_a02 = ds.variables["oin_a02"][:, level, i, j]
    oin_a03 = ds.variables["oin_a03"][:, level, i, j]
    oin_a04 = ds.variables["oin_a04"][:, level, i, j]
    oin= oin_a01+oin_a02+oin_a03+oin_a04
    oin=oin/K

    no3_a01 = ds.variables["no3_a01"][:, level, i, j]
    no3_a02 = ds.variables["no3_a02"][:, level, i, j]
    no3_a03 = ds.variables["no3_a03"][:, level, i, j]
    no3_a04 = ds.variables["no3_a04"][:, level, i, j]
    no3=no3_a01+no3_a02+no3_a03+no3_a04
    no3=no3/K
    nh4_a01 = ds.variables["nh4_a01"][:, level, i, j]
    nh4_a02 = ds.variables["nh4_a02"][:, level, i, j]
    nh4_a03 = ds.variables["nh4_a03"][:, level, i, j]
    nh4_a04 = ds.variables["nh4_a04"][:, level, i, j]
    nh4=nh4_a01+nh4_a02+nh4_a03+nh4_a04
    nh4=nh4/K
    PM25=[]
    for timestep in range(0,240,24):
        pm251=pm25[timestep:timestep+24]
        pm251=np.mean(pm251)
        PM25.append(pm251)
        if len(PM25)== 10:
            print("pm2.5")
            PM25= np.asarray(PM25)
            PM25= np.reshape(PM25,(10,1))
            print(PM25)
            PM25=None
    PM10=[]
    for timestep in range(0, 240, 24):
        pm101=pm10[timestep:timestep+24]
        pm102=np.mean(pm101)
        PM10.append(pm102)
        if len(PM10)== 10:
            print("pm10")
            PM10=np.asarray(PM10)
            PM10 = np.reshape(PM10, (10, 1))
            print(PM10)
            PM10=None
    BC=[]
    for timestep in range(0, 240, 24):
        bc1=bc[timestep:timestep+24]
        bc2=np.mean(bc1)
        BC.append(bc2)
        if len(BC)== 10:
            print("bc")
            BC = np.asarray(BC)
            BC = np.reshape(BC, (10, 1))
            print(BC)
            BC=None
    OC=[]
    for timestep in range(0, 240, 24):
        oc1=oc[timestep:timestep+24]
        oc2=np.mean(oc1)
        OC.append(oc2)
        if len(OC)== 10:
            print("oc")
            OC = np.asarray(OC)
            OC = np.reshape(OC, (10, 1))
            print(OC)
            OC=None
    SO4=[]
    for timestep in range(0, 240, 24):
        so41=so4[timestep:timestep+24]
        so42=np.mean(so41)
        SO4.append(so42)
        if len(SO4)== 10:
            print("so4")
            SO4 = np.asarray(SO4)
            SO4 = np.reshape(SO4, (10, 1))
            print(SO4)
            SO4=None
    OIN=[]
    for timestep in range(0, 240, 24):
        oin1=oin[timestep:timestep+24]
        oin2=np.mean(oin1)
        OIN.append(oin2)
        if len(OIN)== 10:
            print("oin")
            OIN = np.asarray(OIN)
            OIN = np.reshape(OIN, (10, 1))
            print(OIN)
            OIN=None
    NO3=[]
    for timestep in range(0, 240, 24):
        no31=no3[timestep:timestep+24]
        no32=np.mean(no31)
        NO3.append(no32)
        if len(NO3)== 10:
            print("no3")
            NO3 = np.asarray(NO3)
            NO3 = np.reshape(NO3, (10, 1))
            print(NO3)
            NO3=None
    NH4=[]
    for timestep in range(0, 240, 24):
        nh41=nh4[timestep:timestep+24]
        nh42=np.mean(nh41)
        NH4.append(nh42)
        if len(NH4)== 10:
            print("nh4")
            NH4 = np.asarray(NH4)
            NH4 = np.reshape(NH4, (10, 1))
            print(NH4)
            NH4=None

    O3=[]
    for timestep in range(0, 240, 24):
        o31=C_o3[timestep:timestep+24]
        o32=np.mean(o31)
        O3.append(o32)
        if len(O3)== 10:
            print("o3")
            O3 = np.asarray(O3)
            O3 = np.reshape(O3, (10, 1))
            print(O3)
            O3=None

    NO=[]
    for timestep in range(0, 240, 24):
        no1=C_no[timestep:timestep+24]
        no2=np.mean(no1)
        NO.append(no2)
        if len(NO)== 10:
            print("no")
            NO = np.asarray(NO)
            NO = np.reshape(NO, (10, 1))
            print(NO)
            NO=None

    NO2=[]
    for timestep in range(0, 240, 24):
        no21=C_no2[timestep:timestep+24]
        no22=np.mean(no21)
        NO2.append(no22)
        if len(NO2)== 10:
            print("no2")
            NO2 = np.asarray(NO2)
            NO2 = np.reshape(NO2, (10, 1))
            print(NO2)
            NO2=None

    SO2=[]
    for timestep in range(0, 240, 24):
        so21=C_so2[timestep:timestep+24]
        so22=np.mean(so21)
        SO2.append(so22)
        if len(SO2)== 10:
            print("so2")
            SO2 = np.asarray(SO2)
            SO2 = np.reshape(SO2, (10, 1))
            print(SO2)
            SO2=None

    CO=[]
    for timestep in range(0, 240, 24):
        co1=C_co[timestep:timestep+24]
        co2=np.mean(co1)
        CO.append(co2)
        if len(CO)== 10:
            print("CO")
            CO = np.asarray(CO)
            CO = np.reshape(CO, (10, 1))
            print(CO)
            CO=None

    NH3=[]
    for timestep in range(0, 240, 24):
        nh31=C_nh3[timestep:timestep+24]
        nh32=np.mean(nh31)
        NH3.append(nh32)
        if len(NH3)== 10:
            print("nh3")
            NH3 = np.asarray(NH3)
            NH3 = np.reshape(NH3, (10, 1))
            print(NH3)
            NH3=None


    print("location over")
    print("\n \n \n \n \n \n \n \n")
#np.savetxt(sys.stdout.buffer, precp,fmt="%.3f")
#sys.stdout.close()

