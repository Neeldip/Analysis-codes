import pandas as pd
import matplotlib.pyplot as pt

df = pd.read_excel('D:\File2.xlsx')
x = df['Date_Time']
#print(type(df['Date_Time']))
y = df['Hum','Temp','WindSpd']
pt.plot(x,y)
pt.show()
