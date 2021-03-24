import pickle
import numpy as np
import pandas as pd
from pathlib import Path
import sys
np.set_printoptions(threshold=sys.maxsize)

p = Path('.')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

disp = 0.001/0.529177 #bohr
p_n = 0
xyz = 0

first_vext = np.zeros((11,3,103,103))

for d in dirs:

    atom_num = int(d.name[19:21])

    if d.name[22] == "X":
        xyz = 0
    elif d.name[22] == "Y":
        xyz = 1
    elif d.name[22] == "Z":
        xyz = 2

    if d.name[24] == "P":
        p_n = 1
    elif d.name[24] == "N":
        p_n = -1   

    with open(d.name + "/Core_Mat.pkl",'rb') as f:
        first_vext[atom_num-1][xyz] += p_n*np.array(pickle.load(f),dtype=float)/(2*disp)

#VKS = VH + VXC + VEXT
#VH = VKS - KE - PE - VXC
# print(first_vext)
with open("First-Derivative-Vext.pkl",'wb') as f:
    pickle.dump(first_vext,f)
