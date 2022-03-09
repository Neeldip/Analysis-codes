from netCDF4 import Dataset
from wrf import ll_to_xy,xy_to_ll

ds = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
ll_corner= ll_to_xy(ds,26.75,88.5,timeidx=0)
print(ll_corner)
#ul_corner=ll_to_xy(ds,30.0,98.0,timeidx=0)
#print(ul_corner)

ll=xy_to_ll(ds,182,157)
#print(ll)