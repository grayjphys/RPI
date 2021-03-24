import numpy as np
import pandas as pd
import os
import pickle
from pathlib import Path
pd.set_option('display.max_rows',20)


p = Path('.')
files = np.sort([x for x in p.iterdir() if ("efield" in x.name and "_0" in x.name)])



file = open(files[0])
data = []
for i,line in enumerate(file):
    if i > 16:
        for x in line.split():
            data.append(np.float(x))
data = np.array(data)
L = round(data.shape[0]**(1/3))
data = data.reshape(L,L,L)

with open("EX.pkl",'wb') as f:
    pickle.dump(data,f)


#How to read file

#with open("EX.pkl",'rb') as f:
    #data = pickle.load(f)


file = open(files[1])
data = []
for i,line in enumerate(file):
    if i > 16:
        for x in line.split():
            data.append(np.float(x))
data = np.array(data)
L = round(data.shape[0]**(1/3))
data = data.reshape(L,L,L)

with open("EY.pkl",'wb') as f:
    pickle.dump(data,f)


file = open(files[2])
data = []
for i,line in enumerate(file):
    if i > 16:
        for x in line.split():
            data.append(np.float(x))
data = np.array(data)
L = round(data.shape[0]**(1/3))
data = data.reshape(L,L,L)

with open("EZ.pkl",'wb') as f:
    pickle.dump(data,f)