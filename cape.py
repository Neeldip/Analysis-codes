import wrf as w
from wrf import omp_set_num_threads, to_np
import numpy as np
omp_set_num_threads(4)
import matplotlib.pyplot as plt
import netCDF4 as nc
ds = nc.Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
bpres=w.getvar(ds,"P",timeidx=203)
ppress=w.getvar(ds,"PB",timeidx=203)
press=(bpres+ppress)/100
T=w.getvar(ds,"THM",timeidx=203)
tkel=T+300
qv=w.getvar(ds,"QVAPOR",timeidx=203)
PH=w.getvar(ds,"PH",timeidx=203)
PHB=w.getvar(ds,"PHB",timeidx=203)
PHB = w.destagger(PHB, stagger_dim=0, meta=False)
PH = w.destagger(PH, stagger_dim=0, meta=False)
height=(PH+PHB)/9.81
terrain=w.getvar(ds,"HGT",timeidx=203)
psfc=w.getvar(ds,"PSFC",timeidx=203)
psfc_hpa=psfc/100
cape=w.cape_2d(press,tkel,qv,height,terrain,psfc_hpa,ter_follow=True,meta=False)
print(np.shape(cape))