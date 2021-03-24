import pickle
import numpy as np
import pandas as pd
from pathlib import Path
import sys
np.set_printoptions(threshold=sys.maxsize)

p = Path('.')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

p_n = 0
xyz = 0

VXC = np.zeros((2,11,3,103,103))
PE = np.zeros((2,11,3,103,103))
VKS = np.zeros((2,11,3,103,103))
KE = np.zeros((2,11,3,103,103))
Vext = np.zeros((2,11,3,103,103))
first_vext = np.zeros((11,3,103,103))

for d in dirs:

    atom_num = int(d.name[19:21])

    if d.name[22] == "X":
        xyz = 0
    elif d.name[22] == "Y":
        xyz = 1
    elif d.name[22] == "Z":
        xyz = 2

    if d.name[24] == "p":
        p_n = 0
    elif d.name[24] == "n":
        p_n = 1   

    with open(d.name + "/PE_Mat.pkl",'rb') as f:
        PE[p_n][atom_num-1][xyz] = np.array(pickle.load(f),dtype=float)

    with open(d.name + "/VXC_Mat.pkl",'rb') as f:
        VXC[p_n][atom_num-1][xyz] = np.array(pickle.load(f),dtype=float)

    with open(d.name + "/KS_Mat.pkl",'rb') as f:
        VKS[p_n][atom_num-1][xyz] = np.array(pickle.load(f),dtype=float)

    with open(d.name + "/KE_Mat.pkl",'rb') as f:
        KE[p_n][atom_num-1][xyz] = np.array(pickle.load(f),dtype=float)

#VKS = VH + VXC + VEXT
#VH = VKS - KE - PE - VXC
 
Vext = KE + PE

first_vext = (Vext[0] - Vext[1])/(2*0.000001)

with open("First-Derivative-Vext.pkl",'wb') as f:
    pickle.dump(first_vext,f)