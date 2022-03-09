import numpy as np
from numpy import *

from matplotlib import pyplot,colors, axes,gridspec,figure
from netCDF4 import Dataset
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 interpline, CoordPair,ALL_TIMES, enable_pyngl)


wrf_file = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrfout_d01_2018-04-10_00%3A00%3A00")
trad = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")

wrf_file1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrfout_d01_2018-04-10_00%3A00%3A00")
trad1 = Dataset(r"H:\WRF_Chem_Output\201\April\MYNN3_no_aer_feedback\wrf_trad_fields_d01_2018-04-10_00%3A00%3A00")



time = ALL_TIMES
#time = 0

list = np.linspace(0,44,1)


#ua = getvar(wrf_file,"ua",timeidx=time)
#va = getvar(wrf_file,"va",timeidx=time)
#ua1 = getvar(wrf_file1,"ua",timeidx=time)
#va1 = getvar(wrf_file1,"va",timeidx=time)

R = 287  ##Universal gas constant,unit=J*kg^-1*K^-1
x = np.linspace(71.679886, 100.43091, 299)
y = np.linspace(10.64794, 31.228256, 231)
# x = np.linspace(87.65714,98.07542,348)
# y = np.linspace(21.82169,30.20221,312)





for level in range(0,20,1):
    X = wrf_file.variables["QVAPOR"][:, level, :, :]
    X1 = wrf_file1.variables["QVAPOR"][:, level, :, :]

    UA = trad.variables["UMET"][:, level, :, :]
    VA = trad.variables["VMET"][:, level, :, :]
    UA1 = trad1.variables["UMET"][:, level, :, :]
    VA1 = trad1.variables["VMET"][:, level, :, :]

    T = trad.variables["TEMPERATURE"][:,level,:,:]
    P = trad.variables["PRESSURE"][:,level,:,:]
    T1 = trad1.variables["TEMPERATURE"][:,level,:,:]
    P1 = trad1.variables["PRESSURE"][:,level,:,:]

    UA = UA.mean(axis=0)
    VA = VA.mean(axis=0)
    UA1 = UA1.mean(axis=0)
    VA1 = VA1.mean(axis=0)

    G = X.mean(axis=0)
    H = T.mean(axis=0)
    I = P.mean(axis=0)
    G1 = X1.mean(axis=0)
    H1 = T1.mean(axis=0)
    I1 = P1.mean(axis=0)

    #U = UA
    #V = VA
    #U1 = UA1
    #V1 = VA1

    #u_diff=U-U1
    #v_diff=V-V1
    #print(np.shape(G))
    #s = G#.squeeze(axis=0)
    #A = H#.squeeze(axis=0)
    #B = I#.squeeze(axis=0)
    #s1 = G1#.squeeze(axis=0)
    #A1 = H1#.squeeze(axis=0)
    #B1 = I1#.squeeze(axis=0)


    K = (R * H) / I
    K1 = (R * H1) / I1

    L = G/K
    L1 = G1/K1
    u = UA*L
    v = VA*L
    u1 = UA1*L1
    v1 = VA1*L1

    #u = u#.squeeze(axis=0)
    #v = v#.squeeze(axis=0)
    #u1 = u1#.squeeze(axis=0)
    #v1 = v1#.squeeze(axis=0)
    SPD = sqrt(u**2+v**2)
    SPD1 = sqrt(u1**2+v1**2)
    SPD = asarray(SPD)*1000 #g/m2s
    SPD1 = asarray(SPD1)*1000
    SPD_DIFF=((SPD-SPD1)/SPD1)*100
    #SPD = SPD.squeeze(axis=0)
    #print(np.shape(u))
    #print(np.shape(SPD))
    pyplot.figure(figsize=(11, 7))
    k = np.linspace(np.min(SPD_DIFF), 100, 21, endpoint=True)
    #k = np.linspace(np.min(SPD), 50, 30, endpoint=True)
    strm=pyplot.streamplot(x, y, UA, VA, density=6, arrowsize=0.8, integration_direction='forward', color=SPD_DIFF, cmap="jet",linewidth=1.5)
    #pyplot.streamplot(x,y,u,v,density=5,arrowsize=1,integration_direction='forward',color=(0.9, 0.2, 0.5, 0.9),cmap="autumn",linewidth=2)
    pyplot.colorbar(strm.lines,ticks=k)
    pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\WRF level analysis\flux\QVAPORflux_NOR_level="+str(level)+"time="+str(time)+".png",dpi=600)
    pyplot.show()

    #k = np.linspace(np.min(SPD), np.max(SPD1), 21, endpoint=True)
    #strm = pyplot.streamplot(x, y, UA1, VA1, density=7, arrowsize=0.8, integration_direction='forward', color=SPD1,cmap="jet", linewidth=1.5)
    #fig.colorbar(strm.lines, ticks=k)
    #pyplot.savefig(r"F:\WRF-CHEM ANALYSIS\WRF level analysis\flux\QVAPORflux_NOFEED_level=" + str(level) + "time=" + str(time) + ".png",dpi=600)
    #pyplot.show()