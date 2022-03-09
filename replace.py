#from pathlib2 import Path
#path = Path(r"F:\WRF\UPPERAIR\WRFinterp\New1.txt")
#text = path.read_text()
#text = text.replace("<xarray.DataArray ()>")
#path.write_text(text)
import re
f1 = open(r"F:\WRF\UPPERAIR\WRFinterp\ACM2\1.txt", 'r')
f2 = open(r'F:\WRF\UPPERAIR\WRFinterp\ACM2\1_replaced.txt', 'w')
checkWords = ("<xarray.DataArray ()>","array(",",","dtype=float32)","Coordinates:","XLONG    float32 91.59026","XLAT     float32 26.100395","XTIME    float32","Time     datetime64[ns] 2018-04-01T01:00:00",\
              "level    int32 86 <xarray.DataArray ()>",",","Attributes:","FieldType:      104","units:          kg kg-1",\
              "stagger:","coordinates:    XLONG XLAT XTIME","projection:     Mercator(stand_lon=92.50599670410156, moad_cen_lat=26.005...",\
              "missing_value:  9.969209968386869e+36","_FillValue:     9.969209968386869e+36","vert_units:     None","level    int32",")","<xarray.DataArray 'T_interp' (>",\
              "<xarray.DataArray 'QVAPOR_interp' (>","projection:     Mercator(stand_lon=92.50599670410156 moad_cen_lat=26.005...","23040.0","Time","datetime64[ns] 2018-04-01",r"C:\Users\neeldip\Anaconda3\python.exe C:/Users/neeldip/PycharmProjects/WRF/VintpWRF.py","Process finished with exit code 0","<xarray.DataArray 'ua_interp' (>",\
              "    units:          m s-1","<xarray.DataArray 'va_interp' (>","23760.0","T12:00:00","44640.0","datetime64[ns] 2018-04-16","45360.0","datetime64[ns] 2018-04-16T12:00:00","64800.0","datetime64[ns] 2018-04-30","65520.0","datetime64[ns] 2018-04-30T12:00:00")
repWords = ("","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","")

for line in f1:
    for check, rep in zip(checkWords, repWords):
        line = line.replace(check, rep)
    f2.write(line)
f1.close()
f2.close()