import numpy as np
import pandas as pd
import os
import pickle
pd.set_option('display.max_rows',20)

file = open("MOLog")
skips = []
for i,line in enumerate(file):
    if not (" C " in line or " N " in line or " H " in line):
        skips.append(i)
df = pd.read_csv("MOLog",sep="\t",skiprows=skips,header=None)


df = df.iloc[:,0].str.split(expand=True)


DF = df[0:109]
for i in range(1,78):
    new_df=df[i*109:i*109+109]
    new_df.reset_index(drop=True,inplace=True)
    DF=pd.concat([DF,new_df],axis=1).reindex(DF.index)



DF.drop([0,1,2,3], axis=1,inplace=True)
DF.columns = [np.arange(0,DF.shape[1])]
DF.drop([103,207,311],axis=1,inplace=True)


DF1=DF.iloc[:,:103]
DF2=DF.iloc[:,103:206]
DF2.columns = [np.arange(0,DF2.shape[1])]
DF3=DF.iloc[:,206:]
DF3.columns = [np.arange(0,DF3.shape[1])]

with open("States.pkl",'wb') as f:
    pickle.dump(DF3,f)



