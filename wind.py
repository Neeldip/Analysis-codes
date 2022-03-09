from windrose import WindroseAxes
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import numpy as np
import  pandas as pd
import windrose as w
#df = pd.read_excel(r"G:\WRF_Outputs\Premonsoon\Premonsoon2018_with_nudging_ACM2/speed.xlsx", sheet_name="Sheet2")
#df = pd.read_excel(r"F:\WRF\WIND\U_V\ACM2\ACM2.xlsx", sheet_name="Sheet1")
#df = pd.read_excel(r"I:\observation_data\APRIL\DATA\Hourly_surface_data\2018\DS3505\Guwahati\DS3505\INTERPOLATION.xlsx", sheet_name="time")
#df = pd.read_excel(r"I:\observation_data\JULY\halfhourlydata\2018\GUWAHATI\Interpolated_July.xlsx", sheet_name="Sheet1")
df = pd.read_excel("D:\PHD\My PhD Reports\Manuscripts\Speedndir_APRIL.xlsx",sheet_name="speed")
df1 = pd.read_excel("D:\PHD\My PhD Reports\Manuscripts\Speedndir_APRIL.xlsx",sheet_name="direction")
#df = pd.read_excel(r"G:\WRF_Outputs\Monsoon\ACM2\speed.xlsx", sheet_name="Sheet2")
#df = pd.read_excel(r"I:\asos\GHY\ghyapril2018.xlsx", sheet_name="Sheet3")
ws = df["observed"]
wd = df1["observed"]

ax = WindroseAxes.from_ax()
#help(ax.bar)

ax.bar(wd, ws,normed=True,bins=np.arange(0, 10, 2),lw=10, cmap=cm.get_cmap('gist_rainbow'),nsector=36,blowto=True)
#ax.bar(wd, ws,normed=True,bins=6,lw=10, cmap=cm.get_cmap('gist_rainbow'),nsector=36)
ax.set_legend(title_fontsize=18,borderaxespad=-0.10, loc='lower right',bbox_to_anchor=(1.3, 0.3),units='m/s') #(x, y, width, height)
plt.title("ACM2",y=1.05)
plt.show()


