import numpy as np 
from pathlib import Path 
import os

cwd = Path('.')
mol = 'pyridine'
Mol = 'PYRIDINE'
dirs = [x.name for x in cwd.iterdir() if Mol in x.name and x.is_dir()]
cell = float(os.popen('grep ABC ' + Mol + '-UNPERTURBED/GEO-OPT.inp').read().split()[1])

old_cell = 30

Carb = 12.011  # mass in amu
Nit = 14.007  # mass in amu
Hel = 1.008  # mass in amu

masses = [Carb,Carb,Carb,Carb,Carb,Nit,Hel,Hel,Hel,Hel,Hel]
names = ['C','C','C','C','C','N','H','H','H','H','H']

for d in dirs:
    pos = np.genfromtxt(d + "/pyridine.xyz",skip_header=2,usecols=(1,2,3)).T

    c_x_ref = np.sum([masses[i]*pos[0][i] for i in range(len(masses))])/np.sum(masses)
    c_y_ref = np.sum([masses[i]*pos[1][i] for i in range(len(masses))])/np.sum(masses)
    c_z_ref = np.sum([masses[i]*pos[2][i] for i in range(len(masses))])/np.sum(masses)

    c_x_ref_dist = cell/2 - c_x_ref
    c_y_ref_dist = cell/2 - c_y_ref
    c_z_ref_dist = cell/2 - c_z_ref

    ref_dist = np.array([c_x_ref_dist,c_y_ref_dist,c_z_ref_dist])

    pos = pos.T

    for i,p in enumerate(pos):
        pos[i] += ref_dist

    with open(d + "/pyridine.xyz",'w') as f:
        f.write(f'{len(masses)}\n')
        f.write(f'{mol}\n')
        for i,p in enumerate(pos):
            f.write(f'{names[i]} {p[0]} {p[1]} {p[2]}\n')
