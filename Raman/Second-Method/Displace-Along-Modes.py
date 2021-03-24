import numpy as np
from pathlib import Path
import sys
import os
np.set_printoptions(threshold=sys.maxsize)

p = Path('.')
dirs = [x for x in p.iterdir() if x.is_dir() and "MODE" in x.name]
dirs = np.sort(dirs)

disp = 0.0001
bohr_to_angstrom = 1.8897259886

vec = []
modes_disp = []
mode = []

with open("MODE-DATA") as f:
    for i,line in enumerate(f):
        if i % 11 == 0 and i != 0:
            modes_disp.append(mode)
            mode = []
        vec = np.array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
        mode.append(vec)
    modes_disp.append(mode)
    modes_disp = np.array(modes_disp)

# print(modes_disp)
modes_disp *= disp*bohr_to_angstrom

for d in dirs:
    lines = []
    new_line = ""
    mode_num = int(d.name[-2:])
    P_N = d.name[9]
    pos_neg = 0

    if P_N == "P":
        pos_neg = 1
    else:
        pos_neg = -1
    
    print(mode_num,pos_neg)

    os.system("cp pyridine.xyz " + d.name + "/pyridine.xyz")

    with open("pyridine.xyz",'r') as f:
        new_line = ""
        for i,line in enumerate(f):
            if len(line.split()) == 4:
                atom = line.split()[0]
                change = pos_neg*modes_disp[mode_num-1][i-2]
                vec = np.array([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]) + change
                new_line += "  {0}        {1:.10f}       {2:.10f}       {3:.10f}\n".format(atom,vec[0],vec[1],vec[2])
            else:
                new_line += line
    with open(d.name + "/pyridine.xyz",'w') as f:
        f.write(new_line)
