import netCDF4 as nc
import wrf as w
import numpy as np
wrf = nc.Dataset("G:\WRF_Outputs\Premonsoon2018_with_nudging_ACM2\wrfout_d03_2018-03-31_01%3A00%3A00.nc")

z = [80]
x1_00_23 =[86,173,252,297,314,513,596,763,782,895,981,1490,1530,1961,2003,2087,2454,2935,3109,3144,3250,3345,3783,3869,4187,4344,4569,5036,5753,5921,6335,6744,7413,7449,7780,8087,8585,9076,9433,9803,10663,11007,11486,11604,11874,12093,12158,12806,13133,13247,13757,13923,14048,14396,14712,15145,15819,16393,16563,16631,16762,17271,17650,17799,18110,18523,18565,18696,18937,19086,19153,19449,19835]
x1_12_35 =[74,288,593,745,897,1248,1378,1479,1776,1859,2092,2178,3041,3112,3652,3776,3940,4157,4196,4353,4433,4538,5634,5773,5819,6622,6722,7083,7443,7461,7625,8003,8350,9079,9463,9787,10604,10683,11087,11310,12123,12251,12619,13943,14501,14683,15409,15780,16413,17110,17315,17528,18063,18165,18375,18543,18991,19151,19180,19681,19936]
x16_00_383=[73,134,294,585,594,659,753,894,1272,1442,1483,1696,1778,2091,2273,2527,2991,3015,3085,3108,3564,4159,4198,4356,5072,5303,5450,5554,5721,5752,5783,7025,7445,7463,7517,8530,8674,8907,8972,9513,9537,9656,10636,10803,11637,11983,12313,12983,14143,14516,15338,15387,16381,16553,17404,17552,17976,18132,18530,18603,18769,19207,19498]
x16_12_395=[287,593,726,745,896,1479,1611,1683,1745,1901,2006,2091,2340,2696,2982,3051,3110,3311,3432,3740,4031,4186,4357,4571,4944,5346,5524,5783,6374,7111,7146,7473,7996,8601,8924,8991,9257,9533,9722,10348,10532,10803,10941,11831,11988,12149,12313,13621,14153,15071,15698,16573,16602,17305,17945,18047,18444,18613,18779,18881,19163,19218,19407,20122]
x30_00_719=[75,	146,190,296,431,596,735,829,895,1401,1431,1461,1601,1785,1952,2088,2489,2747,3082,3611,3636,3924,4194,4326,4554,4786,5170,5478,5753,5799,6434,7443,7841,8682,9326,9444,9513,9942,10089,10291,10602,10763,10817,11409,12017,12049,12146,12243,12308,12887,13548,13705,14073,14583,15541,16523,16847,16945,17192,17263,17538,17740,17747,18036,18501,18509,18593,18878,19034,19222,19515,19771,20041]
x30_12_731=[74,294,502,594,659,734,894,1194,1214,1462,1561,2087,2925,3091,3607,3984,4112,4191,4336,5139,5330,5588,5773,5835,6041,6416,6500,6703,6982,7417,7473,7777,8309,8861,9469,9563,10823,11015,11321,12303,12693,14072,14113,14194,14706,15668,16543,16584,16842,16910,17608,17699,18066,18165,18434,18540,18623,19546,19982,20168]

for timeidx in range(33,34,1):             ###(23,24,1)-->1.4.2018_00 (35,36,1)-->1.4.2018_12  (383,384,1)-->16.4.2018_00  (395,396,1)-->16.4.2018_12 (719,720,1)-->30.4.2018_00 (731,732,1)-->30.4.2018_12
    for y in x1_00_23:                     ###use x1_00_23 for 1.4.2018_00 and so on. *23--> time step index
           #time = w.getvar(wrf, 'Times', timeidx)
           #print(time)
           U = w.getvar(wrf,'ua', timeidx)
           V = w.getvar(wrf,'va', timeidx)
           T = w.getvar(wrf,'T', timeidx)
           #P = w.getvar(wrf,'P', timeidx)
           #PB = w.getvar(wrf,'PB', timeidx)
           PHB = w.getvar(wrf,"PHB", timeidx)
           PH = w.getvar(wrf,"PH", timeidx)
           #T = w.getvar(wrf, 'T', timeidx)
           #QV = w.getvar(wrf,'QVAPOR', timeidx)
           PHB = w.destagger(PHB, stagger_dim=0, meta=False)
           PH = w.destagger(PH, stagger_dim=0, meta=False)
           GPH = (PH+PHB)/9.81
           #U = w.interplevel(field3d=U, vert=GPH, desiredlev=y)
           #V = w.interplevel(field3d=V, vert=GPH, desiredlev=y)
           T = w.interplevel(field3d=T, vert=GPH, desiredlev=y)
           T = T+300
           #SPD = np.sqrt(U[:, :] ** 2 + V[:, :] ** 2)
           #U157 = U[156,131]
           #V157 = V[156,131]
           T157 = T[208, 245]
           w.eth()
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
           print(T157)

print("##########################################1_00_23end#############################################################################################")

for timeidx in range(35,36,1):  ###(23,24,1)-->1.4.2018_00 (35,36,1)-->1.4.2018_12  (383,384,1)-->16.4.2018_00  (395,396,1)-->16.4.2018_12 (719,720,1)-->30.4.2018_00 (731,732,1)-->30.4.2018_12
       for y in x1_12_35:  ###use x1_00_23 for 1.4.2018_00 and so on. *23--> time step index
              # time = w.getvar(wrf, 'Times', timeidx)
              # print(time)
              U = w.getvar(wrf, 'ua', timeidx)
              V = w.getvar(wrf, 'va', timeidx)
              # P = w.getvar(wrf,'P', timeidx)
              # PB = w.getvar(wrf,'PB', timeidx)
              PHB = w.getvar(wrf, "PHB", timeidx)
              PH = w.getvar(wrf, "PH", timeidx)
              # T = w.getvar(wrf, 'T', timeidx)
              # QV = w.getvar(wrf,'QVAPOR', timeidx)
              PHB = w.destagger(PHB, stagger_dim=0, meta=False)
              PH = w.destagger(PH, stagger_dim=0, meta=False)
              GPH = (PH + PHB) / 9.81
              U = w.interplevel(field3d=U, vert=GPH, desiredlev=y)
              V = w.interplevel(field3d=V, vert=GPH, desiredlev=y)
              # SPD = np.sqrt(U[:, :] ** 2 + V[:, :] ** 2)
              U157 = U[156, 131]
              V157 = V[156, 131]
              # DIR = (270-np.rad2deg(np.arctan2(V157,U157)))%360

              # P = w.interplevel(field3d=P, vert=GPH, desiredlev=y)
              # PB = w.interplevel(field3d=PB, vert=GPH, desiredlev=y)
              # Press = P+PB

              # T = w.interplevel(field3d=T, vert=GPH, desiredlev=y)
              # T = T+300

              # QV = w.interplevel(field3d=QV, vert=GPH, desiredlev=y)

              # SPD = SPD.rename({'south_north':'latitude','west_east':"longitude"})
              # print(SPD[156,131],DIR,Press[156,131],T[156,131],QV[156,131])
              print(U157, V157)

print("##########################################1_12_35end#############################################################################################")

for timeidx in range(383,384,1):  ###(23,24,1)-->1.4.2018_00 (35,36,1)-->1.4.2018_12  (383,384,1)-->16.4.2018_00  (395,396,1)-->16.4.2018_12 (719,720,1)-->30.4.2018_00 (731,732,1)-->30.4.2018_12
       for y in x16_00_383:  ###use x1_00_23 for 1.4.2018_00 and so on. *23--> time step index
              # time = w.getvar(wrf, 'Times', timeidx)
              # print(time)
              U = w.getvar(wrf, 'ua', timeidx)
              V = w.getvar(wrf, 'va', timeidx)
              # P = w.getvar(wrf,'P', timeidx)
              # PB = w.getvar(wrf,'PB', timeidx)
              PHB = w.getvar(wrf, "PHB", timeidx)
              PH = w.getvar(wrf, "PH", timeidx)
              # T = w.getvar(wrf, 'T', timeidx)
              # QV = w.getvar(wrf,'QVAPOR', timeidx)
              PHB = w.destagger(PHB, stagger_dim=0, meta=False)
              PH = w.destagger(PH, stagger_dim=0, meta=False)
              GPH = (PH + PHB) / 9.81
              U = w.interplevel(field3d=U, vert=GPH, desiredlev=y)
              V = w.interplevel(field3d=V, vert=GPH, desiredlev=y)
              # SPD = np.sqrt(U[:, :] ** 2 + V[:, :] ** 2)
              U157 = U[156, 131]
              V157 = V[156, 131]
              # DIR = (270-np.rad2deg(np.arctan2(V157,U157)))%360

              # P = w.interplevel(field3d=P, vert=GPH, desiredlev=y)
              # PB = w.interplevel(field3d=PB, vert=GPH, desiredlev=y)
              # Press = P+PB

              # T = w.interplevel(field3d=T, vert=GPH, desiredlev=y)
              # T = T+300

              # QV = w.interplevel(field3d=QV, vert=GPH, desiredlev=y)

              # SPD = SPD.rename({'south_north':'latitude','west_east':"longitude"})
              # print(SPD[156,131],DIR,Press[156,131],T[156,131],QV[156,131])
              print(U157, V157)

print("##########################################16_00_383end#############################################################################################")

for timeidx in range(395,396,1):  ###(23,24,1)-->1.4.2018_00 (35,36,1)-->1.4.2018_12  (383,384,1)-->16.4.2018_00  (395,396,1)-->16.4.2018_12 (719,720,1)-->30.4.2018_00 (731,732,1)-->30.4.2018_12
       for y in x16_12_395:  ###use x1_00_23 for 1.4.2018_00 and so on. *23--> time step index
              # time = w.getvar(wrf, 'Times', timeidx)
              # print(time)
              U = w.getvar(wrf, 'ua', timeidx)
              V = w.getvar(wrf, 'va', timeidx)
              # P = w.getvar(wrf,'P', timeidx)
              # PB = w.getvar(wrf,'PB', timeidx)
              PHB = w.getvar(wrf, "PHB", timeidx)
              PH = w.getvar(wrf, "PH", timeidx)
              # T = w.getvar(wrf, 'T', timeidx)
              # QV = w.getvar(wrf,'QVAPOR', timeidx)
              PHB = w.destagger(PHB, stagger_dim=0, meta=False)
              PH = w.destagger(PH, stagger_dim=0, meta=False)
              GPH = (PH + PHB) / 9.81
              U = w.interplevel(field3d=U, vert=GPH, desiredlev=y)
              V = w.interplevel(field3d=V, vert=GPH, desiredlev=y)
              # SPD = np.sqrt(U[:, :] ** 2 + V[:, :] ** 2)
              U157 = U[156, 131]
              V157 = V[156, 131]
              # DIR = (270-np.rad2deg(np.arctan2(V157,U157)))%360

              # P = w.interplevel(field3d=P, vert=GPH, desiredlev=y)
              # PB = w.interplevel(field3d=PB, vert=GPH, desiredlev=y)
              # Press = P+PB

              # T = w.interplevel(field3d=T, vert=GPH, desiredlev=y)
              # T = T+300

              # QV = w.interplevel(field3d=QV, vert=GPH, desiredlev=y)

              # SPD = SPD.rename({'south_north':'latitude','west_east':"longitude"})
              # print(SPD[156,131],DIR,Press[156,131],T[156,131],QV[156,131])
              print(U157, V157)

print("##########################################16_12_395end#############################################################################################")

for timeidx in range(719,720,1):  ###(23,24,1)-->1.4.2018_00 (35,36,1)-->1.4.2018_12  (383,384,1)-->16.4.2018_00  (395,396,1)-->16.4.2018_12 (719,720,1)-->30.4.2018_00 (731,732,1)-->30.4.2018_12
       for y in x30_00_719:  ###use x1_00_23 for 1.4.2018_00 and so on. *23--> time step index
              # time = w.getvar(wrf, 'Times', timeidx)
              # print(time)
              U = w.getvar(wrf, 'ua', timeidx)
              V = w.getvar(wrf, 'va', timeidx)
              # P = w.getvar(wrf,'P', timeidx)
              # PB = w.getvar(wrf,'PB', timeidx)
              PHB = w.getvar(wrf, "PHB", timeidx)
              PH = w.getvar(wrf, "PH", timeidx)
              # T = w.getvar(wrf, 'T', timeidx)
              # QV = w.getvar(wrf,'QVAPOR', timeidx)
              PHB = w.destagger(PHB, stagger_dim=0, meta=False)
              PH = w.destagger(PH, stagger_dim=0, meta=False)
              GPH = (PH + PHB) / 9.81
              U = w.interplevel(field3d=U, vert=GPH, desiredlev=y)
              V = w.interplevel(field3d=V, vert=GPH, desiredlev=y)
              # SPD = np.sqrt(U[:, :] ** 2 + V[:, :] ** 2)
              U157 = U[156, 131]
              V157 = V[156, 131]
              # DIR = (270-np.rad2deg(np.arctan2(V157,U157)))%360

              # P = w.interplevel(field3d=P, vert=GPH, desiredlev=y)
              # PB = w.interplevel(field3d=PB, vert=GPH, desiredlev=y)
              # Press = P+PB

              # T = w.interplevel(field3d=T, vert=GPH, desiredlev=y)
              # T = T+300

              # QV = w.interplevel(field3d=QV, vert=GPH, desiredlev=y)

              # SPD = SPD.rename({'south_north':'latitude','west_east':"longitude"})
              # print(SPD[156,131],DIR,Press[156,131],T[156,131],QV[156,131])
              print(U157, V157)

print("##########################################30_00_719end#############################################################################################")

for timeidx in range(731,732, 1):  ###(23,24,1)-->1.4.2018_00 (35,36,1)-->1.4.2018_12  (383,384,1)-->16.4.2018_00  (395,396,1)-->16.4.2018_12 (719,720,1)-->30.4.2018_00 (731,732,1)-->30.4.2018_12
       for y in x30_12_731:  ###use x1_00_23 for 1.4.2018_00 and so on. *23--> time step index
              # time = w.getvar(wrf, 'Times', timeidx)
              # print(time)
              U = w.getvar(wrf, 'ua', timeidx)
              V = w.getvar(wrf, 'va', timeidx)
              # P = w.getvar(wrf,'P', timeidx)
              # PB = w.getvar(wrf,'PB', timeidx)
              PHB = w.getvar(wrf, "PHB", timeidx)
              PH = w.getvar(wrf, "PH", timeidx)
              # T = w.getvar(wrf, 'T', timeidx)
              # QV = w.getvar(wrf,'QVAPOR', timeidx)
              PHB = w.destagger(PHB, stagger_dim=0, meta=False)
              PH = w.destagger(PH, stagger_dim=0, meta=False)
              GPH = (PH + PHB) / 9.81
              U = w.interplevel(field3d=U, vert=GPH, desiredlev=y)
              V = w.interplevel(field3d=V, vert=GPH, desiredlev=y)
              # SPD = np.sqrt(U[:, :] ** 2 + V[:, :] ** 2)
              U157 = U[156, 131]
              V157 = V[156, 131]
              # DIR = (270-np.rad2deg(np.arctan2(V157,U157)))%360

              # P = w.interplevel(field3d=P, vert=GPH, desiredlev=y)
              # PB = w.interplevel(field3d=PB, vert=GPH, desiredlev=y)
              # Press = P+PB

              # T = w.interplevel(field3d=T, vert=GPH, desiredlev=y)
              # T = T+300

              # QV = w.interplevel(field3d=QV, vert=GPH, desiredlev=y)

              # SPD = SPD.rename({'south_north':'latitude','west_east':"longitude"})
              # print(SPD[156,131],DIR,Press[156,131],T[156,131],QV[156,131])
              print(U157, V157)

print("##########################################30_12_731end#############################################################################################")
