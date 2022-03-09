import pandas as pd
import numpy as np
import matplotlib.pyplot as pt

df=pd.read_excel(r'F:\WRF\UPPERAIR\PBL DATA\QNSE\VV.xlsx',
                   sheet_name='Sheet2')
#print(df)
avg_1_00 = df[8518:8638]
d=np.mean(avg_1_00)
print(d)

avg_1_12 = df[12838:12958]
d=np.mean(avg_1_12)
print(d)

avg_16_00 = df[138118:138238]
d=np.mean(avg_16_00)
print(d)

avg_16_12 = df[142438:142558]
d=np.mean(avg_16_12)
print(d)

avg_30_00 = df[259078:259198]
d=np.mean(avg_30_00)
print(d)

avg_30_12 = df[263399:263519]
d=np.mean(avg_30_12)
print(d)