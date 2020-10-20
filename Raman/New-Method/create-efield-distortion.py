import pickle
import numpy as np
import pandas as pd
from pathlib import Path
import sys
import os
np.set_printoptions(threshold=sys.maxsize)


p = Path('.')
dirs = [x for x in p.iterdir() if x.is_dir()]

I = 1E10
I_au = I/3.50944E16

E0 = np.array([np.sqrt(I_au),0,0])

intensity = np.array([0,0,0])
scale = 0.000001
p_n = 0
xyz = np.array([1,0,0])
p_n2 = 0
xyz2 = np.array([1,0,0])

for d in dirs:
    if "CHANGE" in d.name:
        if d.name[24] == "P":
            p_n = 1
        elif d.name[24] == "N":
            p_n = -1

        if d.name[16] == "X":
            xyz = p_n*scale*np.array([1,0,0])
        elif d.name[16] == "Y":
            xyz = p_n*scale*np.array([0,1,0])
        elif d.name[16] == "Z":
            xyz = p_n*scale*np.array([0,0,1])

        E = E0 + xyz
        norm = np.linalg.norm(E)
        #I = np.sqrt(np.linalg.norm(E)/3.50944E16)
        I_new = norm**2
        I_new = I_new*3.50944E16
        P = E/norm

        f = open(str(d) + "/GEO-OPT.inp",'r')
        f2 = open(str(d) + "/new-GEO-OPT.inp",'w')
        for line in f:
            if "INTENSITY" in line:
                f2.write("      INTENSITY " + str(f"{I_new:.16f}") + "\n")
            elif "POLARISATION" in line:
                f2.write("      POLARISATION " + str(f"{P[0]:.16f}") + " " + str(f"{P[1]:.16f}") + " " + str(f"{P[2]:.16f}") + "\n")

            else:
                f2.write(line)
        os.system("mv " + str(d) + "/new-GEO-OPT.inp " + str(d) + "/GEO-OPT.inp")
        ddirs = [x for x in d.iterdir() if x.is_dir()]
        
        for dd in ddirs:
            if "CHANGE" in dd.name:
                if dd.name[29] == "P":
                    p_n2 = 1
                elif dd.name[29] == "N":
                    p_n2 = -1

                if dd.name[17] == "X":
                    xyz2 = p_n2*scale*np.array([1,0,0])
                elif dd.name[17] == "Y":
                    xyz2 = p_n2*scale*np.array([0,1,0])
                elif dd.name[17] == "Z":
                    xyz2 = p_n2*scale*np.array([0,0,1])

                E = E0 + xyz + xyz2
                norm = np.linalg.norm(E)
                #I = np.sqrt(np.linalg.norm(E)/3.50944E16)
                I_new = norm**2
                I_new = I_new*3.50944E16
                P = E/norm

                f = open(str(dd) + "/GEO-OPT.inp",'r')
                f2 = open(str(dd) + "/new-GEO-OPT.inp",'w')
                for line in f:
                    if "INTENSITY" in line:
                        f2.write("      INTENSITY " + str(f"{I_new:.16f}") + "\n")
                    elif "POLARISATION" in line:
                        f2.write("      POLARISATION " + str(f"{P[0]:.16f}") + " " + str(f"{P[1]:.16f}") + " " + str(f"{P[2]:.16f}") + "\n")

                    else:
                        f2.write(line)
                os.system("mv " + str(dd) + "/new-GEO-OPT.inp " + str(dd) + "/GEO-OPT.inp")
