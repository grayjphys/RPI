import pickle
import numpy as np
import pandas as pd
from pathlib import Path
import sys
np.set_printoptions(threshold=sys.maxsize)

p = Path('.')
dirs = [x for x in p.iterdir() if x.is_dir() and "ATOM" in x.name]

for d in dirs:
    lines = []
    change = []
    atom_num = int(d.name[19:21])
    xyz = d.name[22]
    p_n = d.name[-3]
    f  = open(d.name+"/pyridine.xyz",'r')
    for line in f:
        lines.append(line)
    atom = lines[atom_num+1].split()[0]
    vec = np.array([float(lines[atom_num+1].split()[1]),float(lines[atom_num+1].split()[2]),float(lines[atom_num+1].split()[3])])
    if xyz == "X":
        change = np.array([0.000001,0,0])
    elif xyz == "Y":
        change = np.array([0,0.000001,0])
    elif xyz == "Z":
        change = np.array([0,0,0.000001])
    if p_n == "p":
        vec += change
    elif p_n == "n":
        vec -= change
    new_line = atom + " " + str(vec[0]) + " " +str(vec[1]) + " " +str(vec[2]) + "\n"
    new_file = open(d.name+"/new_pyridine.xyz",'w+')
    for l, line in enumerate(lines):
        if l != atom_num + 1:
            new_file.write(line)
        else:
            new_file.write(new_line)

    
    
    