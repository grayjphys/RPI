import pickle
import numpy as np
import pandas as pd
from pathlib import Path
import sys
np.set_printoptions(threshold=sys.maxsize)

p = Path('.')
dirs = [x for x in p.iterdir() if x.is_dir()]

hval = 0.0001
p_n = 0
xyz = 0
p_n2 = 0
xyz2 = 0
dens = np.zeros((103,103))
first_dens = np.zeros((3,2,103,103))
second_dens = np.zeros((3,3,2,2,103,103))

for d in dirs:
    if "UNPERTURBED" in d.name:
        print(d.name)
        with open(d.name+"/Dens_Mat.pkl",'rb') as f:
            dens = np.array(pickle.load(f),dtype=float)
    if "CHANGE" in d.name:
        if d.name[16] == "X":
            xyz = 0
        elif d.name[16] == "Y":
            xyz = 1
        elif d.name[16] == "Z":
            xyz = 2
                    
        if d.name[24] == "P":
            p_n = 0
        elif d.name[24] == "N":
            p_n = 1
        
        with open(d.name + "/Dens_Mat.pkl",'rb') as f:
            first_dens[xyz][p_n] = np.array(pickle.load(f),dtype=float)

        ddirs = [x for x in d.iterdir() if x.is_dir()]
        
        for dd in ddirs:
            if dd.name[17] == "X":
                xyz2 = 0
            elif dd.name[17] == "Y":
                xyz2 = 1
            elif dd.name[17] == "Z":
                xyz2 = 2
                
            if dd.name[29] == "P":
                p_n2 = 0
            elif dd.name[29] == "N":
                p_n2 = 1
            
            with open(d.name + "/" + dd.name + "/Dens_Mat.pkl",'rb') as f:
                second_dens[xyz][xyz2][p_n][p_n2] = np.array(pickle.load(f),dtype=float)

dens_second_derivs = np.zeros((3,3,103,103))

for i in range(3):
    for j in range(3):
        if i == j:
            dens_second_derivs[i][j] = (first_dens[i][0] - 2*dens + first_dens[i][1])/(hval**2)
        else:
            dens_second_derivs[i][j] = (second_dens[i][j][0][0] - second_dens[i][j][0][1] - second_dens[i][j][1][0] + second_dens[i][j][1][1])/(4*hval**2)


with open("Second-Derivative-Dens.pkl",'wb') as f:
    pickle.dump(dens_second_derivs,f)

print(dens_second_derivs)
