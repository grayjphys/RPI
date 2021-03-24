import pickle
from re import S
import numpy as np
from numpy import ma
from numpy.core.defchararray import mod
import pandas as pd
from pathlib import Path
import sys
import matplotlib.pyplot as plt
np.set_printoptions(threshold=sys.maxsize)

class CloseEvent(object):

    def __init__(self):
        self.first = True

    def __call__(self):
        if self.first:
            self.first = False
            return
        sys.exit(0)

p = Path('.')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

broadening = 20
gamma_sqr= (0.5*broadening)**2
atom_num = 11
mode_num = 27
modes = np.zeros((mode_num,atom_num,3)) # in bohr later converted to angstrom
freqs = np.zeros(mode_num) # in cm^-1
C = 12.011
H = 1.00797
N = 14.0067
masses = np.array([C,C,C,C,C,N,H,H,H,H,H]) # in g/mol = amu
line_read = 0
count = 0

with open("PYRIDINE-VIBRATIONS-1.mol",'r') as f:
    for i,line in enumerate(f):
        if atom_num+2 < i < atom_num+mode_num+3:
            freqs[i-3-atom_num] = float(line.split()[0])
        if i > 2*atom_num + mode_num + 4:
            if "vibration" not in line:
                modes[int(count//atom_num)][count%atom_num][0] = float(line.split()[0])
                modes[int(count//atom_num)][count%atom_num][1] = float(line.split()[1])
                modes[int(count//atom_num)][count%atom_num][2] = float(line.split()[2]) 
                count += 1
# print(freqs)
# print(modes)

def lorentzian(x,I):
    L = 0
    for i,v in enumerate(freqs):
        L += I[i]*gamma_sqr/((x-v)**2+gamma_sqr)
    return L

x = np.linspace(500,1700,1000)

with open("Second-Derivative-Dens.pkl",'rb') as f:
    second_dens = np.array(pickle.load(f),dtype=float)

with open("First-Derivative-Vext.pkl",'rb') as f:
    first_vext = np.array(pickle.load(f),dtype=float)

A_mat =  np.zeros((mode_num,3,3))
for m in range(mode_num):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for a in range(atom_num):
                    A_mat[m][i][j] += 2*np.trace(second_dens[i][j].dot(first_vext[a][k]))*modes[m][a][k]/np.sqrt(masses[a])



hc = 1.98644586 # in Jcm and dividing by 10^-23
kT = 414.1947 # in J and dividing by 10^-23

 # in J
I = np.zeros(mode_num)

for m in range(mode_num):
    n = 1/(np.exp(hc*freqs[m]/kT) - 1)
    omega = 2*np.pi*29979245800*freqs[m]
    I[m] = (np.linalg.norm(np.array([1,0,0]).dot(A_mat[m].dot(np.array([0,0,1]))))**2)*(1/omega)*(1+n)
    I[m] /= 1e10

fig = plt.figure()
timer = fig.canvas.new_timer(interval = 3000) #creating a timer object and setting an interval of 3000 milliseconds
timer.add_callback(CloseEvent())

plt.plot(x,lorentzian(x,I),label='Pol X Prop Z')

# print(maxvals)
# annotations = [i for i in range(mode_num) if maxvals[i] > 1.0E-26 and freqs[i] < freqs[-1]+10]
plt.legend(fontsize=15)
plt.xlabel(r'Wavenumber cm$^{-1}$',fontsize=25)
plt.ylabel('A.U.',fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
# for a in annotations:
#     plt.annotate(str(float(round(freqs[a],1))),(freqs[a],maxvals[a]/(1.6*broadening)))
plt.tight_layout()
timer.start()
plt.show()

