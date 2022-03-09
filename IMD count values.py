from xarray import open_dataset
import numpy as np
import pandas as pd

ds = open_dataset(r"F:\IMD\Trend analysis\mann kendall\slope\slope_dec.nc")

latitude = ds.variables["lat"][:]
longitude = ds.variables["lon"][:]

sel_lat = 21.50
sel_lon = 89.25
a = abs(latitude - sel_lat) + abs(longitude - sel_lon)
i, j = np.unravel_index(a.argmin(), a.shape)
# print(i)
# print(j)

sel_lat = 30.0
sel_lon = 99.0
a = abs(latitude - sel_lat) + abs(longitude - sel_lon)
k, l = np.unravel_index(a.argmin(), a.shape)

z = ds.variables["value"][i:k, j:l]
size = np.size(z)
# print(size)
# z = np.reshape(z,(1326,1))
# z= pd.DataFrame(z)
# print(z.values)
# print(len(np.isnan(z)))
nonzero = np.count_nonzero(~np.isnan(z)) ##no of non nan values
# print(nonzero)
#nan = np.count_nonzero(np.isnan(z))  ##no of nan values
# print(nan)
z= np.asarray(z)
#print(z)
#z1 = z[np.where(z>1.96)]
#z1 = np.size(z1)
#z2 = z[np.where( z >0)]
#z2 = np.size(z2)
#z3 = z[np.where(z<0)]
#z3 = np.size(z3)
#z4 = z[np.where(z<-1.96)]
#z4 = np.size(z4)
#z5 = z[np.where(z==0)]
#z5 = np.size(z5)
##Significant increasing trend
#sig_inrease=(z1/nonzero)*100
#print(sig_inrease)

##Non Significant increasing trend
#nonsig_inrease=((z2-z1)/nonzero)*100
#print(nonsig_inrease)

##Non Significant decreasing trend
#nonsig_decrease=((z3-z4)/nonzero)*100
#print(nonsig_decrease)

##Significant decreasing trend
#sig_decrease=(z4/nonzero)*100
#print(sig_decrease)




##increasing trend
#increase=(z2/nonzero)*100
#print(increase)

##decreasing trend
#decrease=(z3/nonzero)*100
#print(decrease)

#no trend
#notrend= (z5/nonzero)*100
#print(notrend)


#increasing slope
#inc_slope = (z2/nonzero)*100
#print(inc_slope)
#decreasing slope
#dec_slope = (z3/nonzero)*100
#print(dec_slope)
#0 slope
#zero_slope =(z5/nonzero)*100
#print(zero_slope)

z1 = z[np.where(z>0)]
z1 = np.size(z1)
z2 = z[np.where( z >5)]
z2 = np.size(z2)
z6 = z[np.where( z >10)]
z6 = np.size(z6)
z3 = z[np.where(z<0)]
z3 = np.size(z3)
z4 = z[np.where(z<-5)]
z4 = np.size(z4)
z7 = z[np.where( z <-10)]
z7 = np.size(z7)
z5 = z[np.where(z==0)]
z5 = np.size(z5)


#increasing slope
inc_slope_0_5 = ((z1-z2)/nonzero)*100
print(inc_slope_0_5)
inc_slope_5_10 = ((z2-z6)/nonzero)*100
print(inc_slope_5_10)
inc_slope_grt10 = (z6/nonzero)*100
print(inc_slope_grt10)

#decreasing slope
dec_slope_0_5 = ((z3-z4)/nonzero)*100
print(dec_slope_0_5)
dec_slope_5_10 = ((z4-z7)/nonzero)*100
print(dec_slope_5_10)
dec_slope_grt10 = (z7/nonzero)*100
print(dec_slope_grt10)

#no slope
no_slope=(z5/nonzero)*100
print(no_slope)