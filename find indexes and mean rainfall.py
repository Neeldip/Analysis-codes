
##############
#THIS CODE FIND THE INDEXES WITH +VE OR -VE RAINFALL AND THEN FIND THE MEAN OF THE +VE OR NEGATIVE RAINFALL INTENSITIES
import numpy as np

import pandas as pd

ds= pd.read_excel("F:\WRF-CHEM ANALYSIS\RAIN_decrease.xlsx",sheet_name='28,93.5',skiprows=1,usecols='O')
total_norm= pd.read_excel("F:\WRF-CHEM ANALYSIS\RAIN_decrease.xlsx",sheet_name='28,93.5',skiprows=1,usecols='F')
total_nofeedback= pd.read_excel("F:\WRF-CHEM ANALYSIS\RAIN_decrease.xlsx",sheet_name='28,93.5',skiprows=1,usecols='N')
#print(total_nofeedback)
#RAINC_norm= pd.read_excel("F:\WRF-CHEM ANALYSIS\RAIN_increase.xlsx",sheet_name='26.5,92.2',skiprows=1,usecols='C')
#RAINNC_norm= pd.read_excel("F:\WRF-CHEM ANALYSIS\RAIN_increase.xlsx",sheet_name='26.5,92.2',skiprows=1,usecols='E')
#RAINC_nofeed= pd.read_excel("F:\WRF-CHEM ANALYSIS\RAIN_increase.xlsx",sheet_name='26.5,92.2',skiprows=1,usecols='K')
#AINNC_nofeed= pd.read_excel("F:\WRF-CHEM ANALYSIS\RAIN_increase.xlsx",sheet_name='26.5,92.2',skiprows=1,usecols='M')

ds=ds.drop(ds.index[[242,243,244,245,246,247,248]])
total_norm=total_norm.drop(total_norm.index[[242,243,244,245,246,247,248]])
total_nofeedback=total_nofeedback.drop(total_nofeedback.index[[242,243,244,245,246,247,248]])

#RAINC_norm=RAINC_norm.drop(ds.index[[242,243,244,245,246,247,248]])
#RAINNC_norm=RAINNC_norm.drop(ds.index[[242,243,244,245,246,247,248]])
#RAINC_nofeed=RAINC_nofeed.drop(ds.index[[242,243,244,245,246,247,248]])
#RAINNC_nofeed=RAINNC_nofeed.drop(ds.index[[242,243,244,245,246,247,248]])

negative_indx=np.where(ds<0)
negative_indx=np.asarray(negative_indx)
negative_indx=negative_indx[0:1,:] ##list of indexes with negative values
negative_indx=np.squeeze(negative_indx)
#print(negative_indx)
positive_indx=np.where(ds>0)
positive_indx=np.asarray(positive_indx)
positive_indx=positive_indx[0:1,:] ##list of indexes with positive values
positive_indx=np.squeeze(positive_indx)
#print(positive_indx)
#print(np.shape(positive_indx))
total_norm=np.squeeze(np.asarray(total_norm))
total_nofeedback=np.squeeze(np.asarray(total_nofeedback))
positive_norm=[]
negative_norm=[]
positive_nofeed=[]
negative_nofeed=[]
for i,j in zip(positive_indx,negative_indx):
    positivenorm=total_norm[i]
    positive_norm.append(positivenorm)
    positivenofeed=total_nofeedback[i]
    positive_nofeed.append(positivenofeed)
    negativenorm = total_norm[j]
    negative_norm.append(negativenorm)
    negativenofeed = total_nofeedback[j]
    negative_nofeed.append(negativenofeed)

#print(np.shape(positive_norm))
#print(positive_norm)
positive_norm_mean=np.mean(positive_norm)
positive_nofeed_mean= np.mean(positive_nofeed)
negative_norm_mean=np.mean(negative_norm)
negative_nofeed_mean= np.mean(negative_nofeed)

print(positive_norm_mean)
print(positive_nofeed_mean)
print(negative_norm_mean)
print(negative_nofeed_mean)
