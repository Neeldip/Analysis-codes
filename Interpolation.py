import pandas as pd
import numpy as np
import matplotlib.pyplot as pt

df = pd.read_excel('E:\DATA\Hourly_surface_data\Guwahati-2018\IITGHY\AP.xlsx',
                   sheet_name='Sheet2', parse_dates=[['Date', 'Time']],)
df = df.set_index(['Date_Time']).resample('5min').last().reset_index()
df.replace(['N','NNE','NE','ENE','E','ESE',
                       'SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW','---','------'],
                      ['0','22.5','45','67.5','90','112.5','135','157.5','180','202.5',
                       '225','247.5','270','292.5','315','337.5','',''], inplace=True)
df = df.set_index(['Date_Time'])
df.replace(np.NaN,'', inplace=True)
df['Temp'] =pd.to_numeric(df['Temp'])
df['Hum'] =pd.to_numeric(df['Hum'])
df['Dewpnt'] =pd.to_numeric(df['Dewpnt'])
df['WindDir'] =pd.to_numeric(df['WindDir'])
df['WindSpd'] =pd.to_numeric(df['WindSpd'])
df['SolarRAd'] =pd.to_numeric(df['SolarRAd'])
df = df.interpolate(method='time')
print(df)
df.to_excel('File2.xlsx', sheet_name='New')



