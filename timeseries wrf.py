import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
ds= pd.read_excel("D:\PHD\My PhD Reports\APS2\matplotlib.xlsx",sheet_name='28,93.5',usecols='A:B',skiprows=0)
#ds1=pd.read_excel("D:\PHD\My PhD Reports\APS2\matplotlib.xlsx",sheet_name='26.5,92.2',usecols='N',skiprows=0)
#ds=ds.drop(ds.index[[242,243,244,245,246,247,248]])
#ds1=ds1.drop(ds1.index[[242,243,244,245,246,247,248]])

ds.plot(colormap='flag')
#ds1.plot(colormap='tab10')
plt.gcf().set_size_inches(15, 5)
plt.xlim(0,240)
#plt.xticks(np.linspace(0,240,1))
#plt.xlabel(np.linspace(1,240,24))
plt.ylim(0,30)
plt.title('(c) Location 6')
plt.ylabel('Rainfall (mm/hr)')
plt.savefig(r'D:\PHD\My PhD Reports\APS2\images\28_93.5.png',dpi=600,bbox_inches='tight')
plt.show()