import pandas as pd
import matplotlib.pyplot as plt

T2 = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\JULY\MATPLOTLIB.xlsx", sheet_name="T2", nrows=1441)
observed = T2["observed"]
QNSE = T2["QNSE"]
ACM2 = T2["ACM2"]
MYNN3 = T2["MYNN3"]
MYJ = T2["MYJ"]
YSU = T2["YSU"]
HONG = T2["Hong"]
plt.hist([observed, QNSE, ACM2, YSU, MYJ, MYNN3, HONG],
         label=["Observed", "QNSE", "ACM2", "YSU", "MYJ", "MYNN3", "HONG"])
plt.legend(loc='upper right')
plt.xlabel("Temperature (K)")
plt.ylabel("Frequency")
plt.title("(b) July")
plt.savefig('HIST_T2_July.png',dpi=300)
plt.show()
RH = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\JULY\MATPLOTLIB.xlsx", sheet_name="RH", nrows=1441)
observed = RH["observed"]
QNSE = RH["QNSE"]
ACM2 = RH["ACM2"]
MYNN3 = RH["MYNN3"]
MYJ = RH["MYJ"]
YSU = RH["YSU"]
HONG = RH["Hong"]
plt.hist([observed, QNSE, ACM2, YSU, MYJ, MYNN3, HONG],
         label=["Observed", "QNSE", "ACM2", "YSU", "MYJ", "MYNN3", "HONG"])
plt.legend(loc='upper left')
plt.xlabel("Relative Humidity (%)")
plt.ylabel("Frequency")
plt.title("(d) July")
plt.savefig('HIST_RH_JULY.png',dpi=300)
plt.show()
