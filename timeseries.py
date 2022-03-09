import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

#T2 = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\JULY\MATPLOTLIB.xlsx", sheet_name="Sheet2", skiprows=1, nrows=8)
T2 = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\JULY\MATPLOTLIB.xlsx", sheet_name="Sheet3", skiprows=0, nrows=49)
# plt.savefig('HIST_T2_July.png',dpi=300)
T2.plot(colormap='tab10')
#plt.xticks([0, 1, 2, 3, 4, 5, 6, 7], ['00', '03', '06', '09', '12', '15', '18', '21'])
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47],
           ['00', '','','','','','03','','','','','', '06','','','','','',
           '09','','','','','', '12','','','','','', '15','','','','','', '18','','','','','', '21','','','','',''])

# plt.show()
plt.gcf().set_size_inches(10, 5)
#plt.title("(a) April ")
plt.title("(b) July ")
plt.xlabel("Time")
#plt.ylabel("Bias (%)")
plt.ylabel("Temperature (K) ")
plt.savefig('T2_July.png',dpi=300,bbox_inches='tight')
plt.show()



#T2 = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\JULY\MATPLOTLIB.xlsx", sheet_name="Sheet2", skiprows=12,
#                   nrows=8)
T2 = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\JULY\MATPLOTLIB.xlsx", sheet_name="Sheet3", skiprows=51,
                   nrows=49)


T2.plot(colormap='tab10')
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47],
           ['00', '','','','','','03','','','','','', '06','','','','','',
           '09','','','','','', '12','','','','','', '15','','','','','', '18','','','','','', '21','','','','',''])
#plt.xticks([0, 1, 2, 3, 4, 5, 6, 7], ['00', '03', '06', '09', '12', '15', '18', '21'])

# plt.show()
plt.gcf().set_size_inches(10, 5)
#plt.title("(c) April ")
plt.title("(d) July ")
plt.xlabel("Time")
#plt.ylabel("RMSE (%)")
plt.ylabel("Relative humidity (%)")

plt.savefig('RH_july.png',dpi=300,bbox_inches='tight')
plt.show()
