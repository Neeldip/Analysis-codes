import wrf as w
from wrf import *
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
from netCDF4 import *

# https://nordicesmhub.github.io/climate-data-tutorial/03-visualization-python/
# trad = Dataset("G:\WRF_Chem_Output\ACM2\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")
# wrf_file = Dataset("G:\WRF_Chem_Output\ACM2\wrfout_d01_2018-04-10_00%3A00%3A00")
# wrf_file = Dataset(r"G:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
wrf_file = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2\wrfout_d03.nc")
wrf_file1 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_HongShin\wrfout_d03.nc")
wrf_file2 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYJ\wrfout_d03.nc")
wrf_file3 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_MYNN3\wrfout_d03.nc")
wrf_file4 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_YSU\wrfout_d03.nc")
wrf_file5 = Dataset("G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_QNSE\wrfout_d03.nc")
list = [wrf_file, wrf_file1, wrf_file2, wrf_file3, wrf_file4, wrf_file5]
time = ALL_TIMES

for file in list:
    QCLOUD = getvar(file, "QCLOUD", timeidx=time)
    Z = QCLOUD
    X = QCLOUD
    if time == ALL_TIMES:
        X = X.mean("Time")
    lats, lons = w.latlon_coords(Z)
    #print(X)
    x = X.sum(dim='bottom_top') * 1000
    print(x)
    v = np.linspace(np.min(x), 1, 30)
    plt.contourf(lons, lats, x,50, cmap=plt.get_cmap("nipy_spectral"),
                 extend='max')  ###use if conc does not cahnege much between plots
    # plt.contourf(lons, lats, CONC, cmap=plt.get_cmap("nipy_spectral"), extend='max')
    font = {'family': 'serif', 'color': 'black', 'weight': 'bold', 'size': 12, }
    plt.xlabel("Longitude", fontdict=font)
    plt.ylabel("Latitude", fontdict=font)
    plt.title("QCLOUD (g/kg)", weight='bold', fontdict=font)

    ##SET TICK LABELS PROPERTIES####
    ax = plt.gca()
    ax.set_xticklabels(ax.get_xticks(), font)
    ax.set_yticklabels(ax.get_yticks(), font)
    # plt.colorbar(shrink=.90)
    pbl = ["ACM2","HONG","MYJ","MYNN3","YSU","QNSE"]
    cb = plt.colorbar(fraction=0.046, pad=0.04).set_label(label='', size=15, weight='bold')
    for name in pbl:
        plt.savefig("QCLOUD_" + str(name) + ".png", dpi=500)
    plt.show()
