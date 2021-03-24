import numpy as np
import pandas as pd
import os
import pickle
pd.set_option('display.max_rows',20)

file = open("D-MATRIX")
skips = []
for i,line in enumerate(file):
    if not (" C " in line or " N " in line or " H " in line):
        skips.append(i)
df = pd.read_csv("D-MATRIX",sep="\t",skiprows=skips,header=None)

df = df.iloc[:,0].str.split(expand=True)

DF = df[0:103]

for i in range(1,52):
    new_df=df[i*103:i*103+103]
    new_df.reset_index(drop=True,inplace=True)
    DF=pd.concat([DF,new_df],axis=1).reindex(DF.index)


DF.drop([0,1,2,3], axis=1,inplace=True)
DF.columns = [np.arange(0,DF.shape[1])]
DF.drop(103,axis=1,inplace=True)
with open("Dens_Mat.pkl",'wb') as f:
    pickle.dump(DF,f)