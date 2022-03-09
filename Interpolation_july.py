import pandas as pd

ts = pd.read_excel(r'I:\observation_data\WRF-Chem\original\patna.xlsx', sheet_name="patna")  # READ EXCEL FILE
#ts['HrMn'].replace([0, 30, 100, 130, 200, 230, 300, 330, 400, 430, 500, 530, 600, 630, 700, 730, 800, 830, 900, 930],
#                   # REPLACE VALUES
#                   ["0000", "0030", "0100", "0130", "0200", "0230", "0300", "0330", "0400", "0430", "0500", "0530",
#                    "0600", "0630", "0700", "0730", "0800", "0830", "0900", "0930"], inplace=True)
#ts["Time"] = ts["Date"].astype(str) + ts["HrMn"].astype(str)  # MERGE COLUMNS INTO "TIME"
ts['Time'] = pd.to_datetime(ts['Time'])  # CONVERT "TIME" STRING TO TIMESERIES
#ts['DIR'].replace(999, "", inplace=True)  # REPLACE "999' VALEUES IN DIR TO NaN
ts.set_index('Time', inplace=True)  # SET THE INDEX TO TIME
#print(ts[ts.index.duplicated()])
dt = pd.date_range("2018-04-10 00:00:00", "2018-04-19 23:30:00", freq='30min',
                   name='Time')  # CREATE NEW TIME SERIES AT 30MIN FREQUENCY
idx = pd.DatetimeIndex(dt)
ts = ts.reindex(idx)
#ts = ts.drop(columns=['Date', 'HrMn', 'SLP'])  # REMOVE THE COLUMNS DATE AND HrMn
ts = ts.interpolate(method='time')
print(ts)
with pd.ExcelWriter(r'I:\observation_data\WRF-Chem\interpolated\patna.xlsx',
                    ) as writer:
    ts.to_excel(writer, sheet_name="akima")  # WRITE TO EXCEL
