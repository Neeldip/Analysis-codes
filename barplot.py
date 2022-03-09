import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import numpy as np

#T2 = pd.read_excel(r"F:\WRF\Timeseries\SURFACE\Guwahati\JULY\MATPLOTLIB.xlsx", sheet_name="Sheet2", skiprows=1, nrows=8)
T2 = pd.read_excel(r"D:\PHD\My PhD Reports\Manuscripts\New.xlsx", sheet_name="APRIL_RAIN", skiprows=0, nrows=15)
#print(T2['Location'])
#locations = ['AG','DB','GH','SI','TZ','AZ','IM','BG','DH','SH','KO','JO','PG','DK']
#ds = T2.groupby("Location")[locations]
QNSE = T2["QNSE"].to_list()
GPM = T2["GPM"].to_list()
MYNN3 = T2["MYNN3"].to_list()
MYJ = T2["MYJ"].to_list()
YSU = T2["YSU"].to_list()
ACM2 = T2["ACM2"].to_list()
HONG = T2["HONG"].to_list()
X = np.arange(len(QNSE))
bar_width =0.1
plt.figure(figsize=(15,5))
plt.bar(X,GPM,width=bar_width,color='blue',zorder=2)
plt.bar(X+bar_width*1,QNSE,width=bar_width,color='coral',zorder=2)
plt.bar(X+bar_width*2,MYNN3,width=bar_width,color='gold',zorder=2)
plt.bar(X+bar_width*3,MYJ,width=bar_width,color='chocolate',zorder=2)
plt.bar(X+bar_width*4,YSU,width=bar_width,color='aqua',zorder=2)
plt.bar(X+bar_width*5,ACM2,width=bar_width,color='yellow',zorder=2)
plt.bar(X+bar_width*6,HONG,width=bar_width,color='grey',zorder=2)

blue = mpatch.Patch(color="blue", label="GPM")
coral = mpatch.Patch(color="coral", label="QNSE")
gold = mpatch.Patch(color="gold", label="MYNN3")
chocolate = mpatch.Patch(color="chocolate", label="MYJ")
aqua = mpatch.Patch(color="aqua", label="YSU")
yellow = mpatch.Patch(color="yellow", label="ACM2")
grey = mpatch.Patch(color="grey", label="HONG")
plt.legend(handles=[blue,coral,gold,chocolate,aqua,yellow,grey])
plt.xticks(X+bar_width*3,['AG','DB','GH','SI','TZ','AZ','IM','BG','DH','SH','KO','JO','PG','DK'])
plt.xlabel("Location")
plt.ylabel("Rainfall (mm)")
plt.title("(a) April")
#plt.grid(axis='y')
plt.savefig('Rainfall_April.png',dpi=300,bbox_inches='tight')
plt.show()