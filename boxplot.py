###BOXPLOT FOR TEMPERATURE
'''
import pandas as pd
import matplotlib.pyplot as plt

T2 = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\APRIL\MATPLOTLIB.xlsx", sheet_name="T2", nrows=1441,)
observed = T2["observed"]
QNSE = T2["QNSE"]
ACM2 = T2["ACM2"]
MYNN3 = T2["MYNN3"]
MYJ = T2["MYJ"]
YSU = T2["YSU"]
HONG = T2["Hong"]
plt.boxplot([observed, QNSE, ACM2, YSU, MYJ, MYNN3, HONG],showmeans=True,showfliers=False)
#plt.legend(loc='upper right')
#plt.xlabel("Temperature (K)")
plt.ylabel("Temperature (K)")
plt.ylim(285,310)
plt.title("(a) April")
plt.xticks([1,2,3,4,5,6,7],['Observed', 'QNSE', 'ACM2', 'YSU', 'MYJ', 'MYNN3', 'HONG'])
plt.savefig('BOX_T2_APRIL.png',dpi=300)
plt.show()
RH = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\APRIL\MATPLOTLIB.xlsx", sheet_name="RH", nrows=1441)
observed = RH["observed"]
QNSE = RH["QNSE"]
ACM2 = RH["ACM2"]
MYNN3 = RH["MYNN3"]
MYJ = RH["MYJ"]
YSU = RH["YSU"]
HONG = RH["Hong"]
plt.boxplot([observed, QNSE, ACM2, YSU, MYJ, MYNN3, HONG],showmeans=True,showfliers=False)
#plt.legend(loc='upper left')
#plt.xlabel("Relative Humidity (%)")
plt.xticks([1,2,3,4,5,6,7],['Observed', 'QNSE', 'ACM2', 'YSU', 'MYJ', 'MYNN3', 'HONG'])
plt.ylabel("Relative Humidity (%)")
plt.ylim(10,105)
plt.title("(c) April")
plt.savefig('BOX_RH_APRIL.png',dpi=300)
plt.show()
'''
'''
###BOXPLOT FOR WIND SPEED
import pandas as pd
import matplotlib.pyplot as plt

T2 = pd.read_excel(r"D:\PHD\My PhD Reports\Manuscripts\Speedndir_APRIL.xlsx", sheet_name="speed", nrows=375,)
#T = pd.read_excel(r"I:\asos\GHY\ghyapril2018.xlsx", sheet_name="Sheet3")
observed = T2["observed"]
print(observed)
QNSE = T2["QNSE"]
ACM2 = T2["ACM2"]
MYNN3 = T2["MYNN3"]
MYJ = T2["MYJ"]
YSU = T2["YSU"]
HONG = T2["Hong"]
plt.boxplot([observed, QNSE, ACM2, YSU, MYJ, MYNN3, HONG],showmeans=True,showfliers=False)
#plt.legend(loc='upper right')
#plt.xlabel("Temperature (K)")
plt.ylabel("Wind speed (m/s)")
plt.ylim(0,10)
plt.title("(a) April")
plt.xticks([1,2,3,4,5,6,7],['Observed', 'QNSE', 'ACM2', 'YSU', 'MYJ', 'MYNN3', 'HONG'])
plt.savefig('BOX_windSpeed_APRIL.png',dpi=300,bbox_inches='tight')
plt.show()

RH = pd.read_excel(r"D:\PHD\My PhD Reports\Manuscripts\Speedndir_JULY.xlsx", sheet_name="speed", nrows=381)
#T = pd.read_excel(r"I:\asos\GHY\ghyjuly2018.xlsx", sheet_name="Sheet1")
observed = RH["observed"]
QNSE = RH["QNSE"]
ACM2 = RH["ACM2"]
MYNN3 = RH["MYNN3"]
MYJ = RH["MYJ"]
YSU = RH["YSU"]
HONG = RH["Hong"]
plt.boxplot([observed, QNSE, ACM2, YSU, MYJ, MYNN3, HONG],showmeans=True,showfliers=False)
#plt.legend(loc='upper left')
#plt.xlabel("Relative Humidity (%)")
plt.xticks([1,2,3,4,5,6,7],['Observed', 'QNSE', 'ACM2', 'YSU', 'MYJ', 'MYNN3', 'HONG'])
plt.ylabel("Wind speed (m/s)")
plt.ylim(0,10)
plt.title("(b) July")
plt.savefig('BOX_windSpeed_JULY.png',dpi=300,bbox_inches='tight')
plt.show()
'''
###BOXPLOT FOR surface parameters
import pandas as pd
import matplotlib.pyplot as plt


RH = pd.read_excel(r"D:\PHD\My PhD Reports\Manuscript2\New.xlsx", sheet_name="SWDOWN", nrows=12,skiprows=0)
#T = pd.read_excel(r"I:\asos\GHY\ghyjuly2018.xlsx", sheet_name="Sheet1")
#observed = RH["observed"]
QNSE = RH["QNSE"]
ACM2 = RH["ACM2"]
MYNN3 = RH["MYNN3"]
MYJ = RH["MYJ"]
YSU = RH["YSU"]
HONG = RH["HONG"]
plt.boxplot([QNSE, ACM2, YSU, MYJ, MYNN3, HONG],showmeans=True,showfliers=False)
#plt.legend(loc='upper left')
#plt.xlabel("Relative Humidity (%)")
plt.xticks([1,2,3,4,5,6],['QNSE', 'ACM2', 'YSU', 'MYJ', 'MYNN3', 'HONG'])
plt.ylabel("Solar radiation (W/m^2) ")
plt.ylim(0,1000)
#plt.title("(a) Day")
#plt.title("(b) Night")
plt.savefig('D:\PHD\My PhD Reports\Manuscript2\SWDOWN.png',dpi=600,bbox_inches='tight')
plt.show()