import netCDF4 as nc
import wrf as w
import numpy as np
wrf = nc.Dataset("G:\WRF_Outputs\Premonsoon2018_with_nudging_ACM2\wrfout_d03_2018-03-31_01%3A00%3A00.nc")

z = [130,150,230,280,330,380,430,480,530,580,630,680,730,780,830,880,930,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2600,2800,3000]
w.eth()

for timeidx in range(23,743,1):             ###(23,24,1)-->1.4.2018_00 (35,36,1)-->1.4.2018_12  (383,384,1)-->16.4.2018_00  (395,396,1)-->16.4.2018_12 (719,720,1)-->30.4.2018_00 (731,732,1)-->30.4.2018_12
    for y in z:                     ###use x1_00_23 for 1.4.2018_00 and so on. *23--> time step index
           #time = w.getvar(wrf, 'Times', timeidx)
           #print(time)
           U = w.getvar(wrf,'ua', timeidx)
           V = w.getvar(wrf,'va', timeidx)
           T = w.getvar(wrf,'T', timeidx)
           #T156 = T[1,156,131]
           P = w.getvar(wrf,'P', timeidx)
           PB = w.getvar(wrf,'PB', timeidx)
           press= P+PB
           PSFC = w.getvar(wrf, 'PSFC', timeidx)
           PHB = w.getvar(wrf,"PHB", timeidx)
           PH = w.getvar(wrf,"PH", timeidx)
           #T = w.getvar(wrf, 'T', timeidx)
           #QV = w.getvar(wrf,'QVAPOR', timeidx)
           PHB = w.destagger(PHB, stagger_dim=0, meta=False)
           PH = w.destagger(PH, stagger_dim=0, meta=False)
           GPH = (PH+PHB)/9.81

           print(GPH[:,56,85])
           #U = w.interplevel(field3d=U, vert=GPH, desiredlev=y)
           #V = w.interplevel(field3d=V, vert=GPH, desiredlev=y)
           #press = w.interplevel(field3d=press, vert=GPH, desiredlev=y)
           #T = w.interplevel(field3d=T, vert=GPH, desiredlev=y)
           #T156 = T[156, 131]
           #TC = ((press[156,131]/PSFC[156,131])**0.286)*(T156+300)
           #SPD = np.sqrt(U[:, :] ** 2 + V[:, :] ** 2)
           #U157 = U[156,131]
           #V157 = V[156,131]
           #TC157 = TC[156, 131]
           #DIR = (270-np.rad2deg(np.arctan2(V157,U157)))%360

           #P = w.interplevel(field3d=P, vert=GPH, desiredlev=y)
           #PB = w.interplevel(field3d=PB, vert=GPH, desiredlev=y)
           #Press = P+PB

           #T = w.interplevel(field3d=T, vert=GPH, desiredlev=y)
           #T = T+300

           #QV = w.interplevel(field3d=QV, vert=GPH, desiredlev=y)

           #SPD = SPD.rename({'south_north':'latitude','west_east':"longitude"})
           #print(SPD[156,131],DIR,Press[156,131],T[156,131],QV[156,131])
           #print(U157,V157)
           #print(TC)