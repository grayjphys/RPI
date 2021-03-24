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
hc = 1.98644586 # in Jcm and dividing by 10^-23
kT = 414.1947 # in J and dividing by 10^-23

line_read = 0

count = 0

with open("Pyridine-VIBRATIONS-1.mol",'r') as f:
    for i,line in enumerate(f):
        if atom_num+2 < i < atom_num+mode_num+3:
            freqs[i-3-atom_num] = float(line.split()[0])
        if i > 2*atom_num + mode_num + 4:
            if "vibration" not in line:
                modes[int(count//atom_num)][count%atom_num][0] = float(line.split()[0])
                modes[int(count//atom_num)][count%atom_num][1] = float(line.split()[1])
                modes[int(count//atom_num)][count%atom_num][2] = float(line.split()[2]) 
                count += 1

def lorentzian(x,I):
    L = 0
    for i,v in enumerate(freqs):
        L += I[i]*gamma_sqr/((x-v)**2+gamma_sqr)
    return L

def R_x(ang_num,ang_div):
    return np.array([[1,0,0],
                     [0,np.cos(2*np.pi*ang_num/ang_div),-np.sin(2*np.pi*ang_num/ang_div)],
                     [0,np.sin(2*np.pi*ang_num/ang_div),np.cos(2*np.pi*ang_num/ang_div)]])
def R_y(ang_num,ang_div):
    return np.array([[np.cos(2*np.pi*ang_num/ang_div),0,np.sin(2*np.pi*ang_num/ang_div)],
                     [0,1,0],
                     [-np.sin(2*np.pi*ang_num/ang_div),0,np.cos(2*np.pi*ang_num/ang_div)]])
def R_z(ang_num,ang_div):
    return np.array([[np.cos(2*np.pi*ang_num/ang_div),-np.sin(2*np.pi*ang_num/ang_div),0],
                     [np.sin(2*np.pi*ang_num/ang_div),np.cos(2*np.pi*ang_num/ang_div),0],
                     [0,0,1]])

def Full_Rot_yz(y_,z_,ang_div,vec):
    return R_y(y_,ang_div).dot(R_z(z_,ang_div).dot(vec))

def Full_Rot_zx(z_,x_,ang_div,vec):
    return R_z(z_,ang_div).dot(R_x(x_,ang_div).dot(vec))

def Full_Rot_xy(x_,y_,ang_div,vec):
    return R_x(x_,ang_div).dot(R_y(y_,ang_div).dot(vec))

x = np.linspace(500,1700,1000)

p = Path('./POL-X-PROP-Z')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

with open(p.name+"/Second-Derivative-Dens.pkl",'rb') as f:
    second_dens = np.array(pickle.load(f),dtype=float)

with open(p.name+"/First-Derivative-Vext.pkl",'rb') as f:
    first_vext = np.array(pickle.load(f),dtype=float)

A_mat_xz =  np.zeros((mode_num,3,3))
for m in range(mode_num):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for a in range(atom_num):
                    A_mat_xz[m][i][j] += 2*np.trace(second_dens[i][j].dot(first_vext[a][k]))*modes[m][a][k]/np.sqrt(masses[a])

p = Path('./POL-Y-PROP-Z')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

with open(p.name+"/Second-Derivative-Dens.pkl",'rb') as f:
    second_dens = np.array(pickle.load(f),dtype=float)

with open(p.name+"/First-Derivative-Vext.pkl",'rb') as f:
    first_vext = np.array(pickle.load(f),dtype=float)

A_mat_yz =  np.zeros((mode_num,3,3))
for m in range(mode_num):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for a in range(atom_num):
                    A_mat_yz[m][i][j] += 2*np.trace(second_dens[i][j].dot(first_vext[a][k]))*modes[m][a][k]/np.sqrt(masses[a])

p = Path('./POL-X-PROP-Y')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

with open(p.name+"/Second-Derivative-Dens.pkl",'rb') as f:
    second_dens = np.array(pickle.load(f),dtype=float)

with open(p.name+"/First-Derivative-Vext.pkl",'rb') as f:
    first_vext = np.array(pickle.load(f),dtype=float)

A_mat_xy =  np.zeros((mode_num,3,3))
for m in range(mode_num):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for a in range(atom_num):
                    A_mat_xy[m][i][j] += 2*np.trace(second_dens[i][j].dot(first_vext[a][k]))*modes[m][a][k]/np.sqrt(masses[a])

p = Path('./POL-Y-PROP-X')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

with open(p.name+"/Second-Derivative-Dens.pkl",'rb') as f:
    second_dens = np.array(pickle.load(f),dtype=float)

with open(p.name+"/First-Derivative-Vext.pkl",'rb') as f:
    first_vext = np.array(pickle.load(f),dtype=float)

A_mat_yx =  np.zeros((mode_num,3,3))
for m in range(mode_num):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for a in range(atom_num):
                    A_mat_yx[m][i][j] += 2*np.trace(second_dens[i][j].dot(first_vext[a][k]))*modes[m][a][k]/np.sqrt(masses[a])

p = Path('./POL-Z-PROP-X')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

with open(p.name+"/Second-Derivative-Dens.pkl",'rb') as f:
    second_dens = np.array(pickle.load(f),dtype=float)

with open(p.name+"/First-Derivative-Vext.pkl",'rb') as f:
    first_vext = np.array(pickle.load(f),dtype=float)

A_mat_zx =  np.zeros((mode_num,3,3))
for m in range(mode_num):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for a in range(atom_num):
                    A_mat_zx[m][i][j] += 2*np.trace(second_dens[i][j].dot(first_vext[a][k]))*modes[m][a][k]/np.sqrt(masses[a])

p = Path('./POL-Z-PROP-Y')
dirs = [x for x in p.iterdir() if x.is_dir() and "MOVE" in x.name]

with open(p.name+"/Second-Derivative-Dens.pkl",'rb') as f:
    second_dens = np.array(pickle.load(f),dtype=float)

with open(p.name+"/First-Derivative-Vext.pkl",'rb') as f:
    first_vext = np.array(pickle.load(f),dtype=float)

A_mat_zy =  np.zeros((mode_num,3,3))
for m in range(mode_num):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for a in range(atom_num):
                    A_mat_zy[m][i][j] += 2*np.trace(second_dens[i][j].dot(first_vext[a][k]))*modes[m][a][k]/np.sqrt(masses[a])

A_mat = (A_mat_xz + A_mat_yz + A_mat_xy + A_mat_yx + A_mat_zx + A_mat_zy)/6
#A_mat = (A_mat_xz + A_mat_yz)/2

A_mat_x = (A_mat_xz+A_mat_xy)/2
A_mat_y = (A_mat_yz+A_mat_yx)/2
A_mat_z = (A_mat_zx+A_mat_zy)/2

I = np.zeros(mode_num)

ang_div = 12

counter_total = 0

for m in range(mode_num):
    n = 1/(np.exp(hc*freqs[m]/kT) - 1)
    omega = 2*np.pi*29979245800*freqs[m]
    counter = 0
    for y1 in range(ang_div):
        for z1 in range(ang_div):
            for y2 in range(ang_div):
                for z2 in range(ang_div):
                    counter += 1
                    counter_total += 1
                    if counter%(ang_div**2) == 0:
                        print(m,100*(counter/(ang_div**4)),100*(counter_total/(27*ang_div**4)))
                    I[m] += (np.linalg.norm(Full_Rot_yz(y1,z1,ang_div,np.array([1,0,0])).dot(A_mat_x[m].dot(Full_Rot_yz(y2,z2,ang_div,np.array([1,0,0])))))**2)*(1/omega)*(1+n)
                    I[m] += (np.linalg.norm(Full_Rot_yz(y1,z1,ang_div,np.array([1,0,0])).dot(A_mat_y[m].dot(Full_Rot_yz(y2,z2,ang_div,np.array([1,0,0])))))**2)*(1/omega)*(1+n)
                    I[m] += (np.linalg.norm(Full_Rot_yz(y1,z1,ang_div,np.array([1,0,0])).dot(A_mat_z[m].dot(Full_Rot_yz(y2,z2,ang_div,np.array([1,0,0])))))**2)*(1/omega)*(1+n)

    I[m] /= 3e10
    I[m] /= counter
#plt.plot(x,lorentzian(x,I),label='Average Incident, Average Scattered')

I = np.zeros(mode_num)

ang_div = 12

counter_total = 0

for m in range(mode_num):
    n = 1/(np.exp(hc*freqs[m]/kT) - 1)
    omega = 2*np.pi*29979245800*freqs[m]
    counter = 0
    for y1 in range(ang_div):
        for z1 in range(ang_div):
            counter += 1
            counter_total += 1
            if counter%(ang_div) == 0:
                print(m,100*(counter/(ang_div**2)),100*(counter_total/(27*ang_div**2)))
            I[m] += (np.linalg.norm(Full_Rot_yz(y1,z1,ang_div,np.array([1,0,0])).dot(A_mat_x[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)
            I[m] += (np.linalg.norm(Full_Rot_yz(y1,z1,ang_div,np.array([1,0,0])).dot(A_mat_y[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)
            I[m] += (np.linalg.norm(Full_Rot_yz(y1,z1,ang_div,np.array([1,0,0])).dot(A_mat_z[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)

    I[m] /= 3e10
    I[m] /= counter
#plt.plot(x,lorentzian(x,I),label='Average Incident, Scattered Pol (1,0,0)')

I = np.zeros(mode_num)

ang_div = 12

counter_total = 0

for m in range(mode_num):
    n = 1/(np.exp(hc*freqs[m]/kT) - 1)
    omega = 2*np.pi*29979245800*freqs[m]
    counter = 0
    for y2 in range(ang_div):
        for z2 in range(ang_div):
            counter += 1
            counter_total += 1
            if counter%(ang_div) == 0:
                print(m,100*(counter/(ang_div**2)),100*(counter_total/(27*ang_div**2)))
            I[m] += (np.linalg.norm(np.array([1,0,0]).dot(A_mat_x[m].dot(Full_Rot_yz(y2,z2,ang_div,np.array([1,0,0])))))**2)*(1/omega)*(1+n)
            I[m] += (np.linalg.norm(np.array([0,1,0]).dot(A_mat_y[m].dot(Full_Rot_yz(y2,z2,ang_div,np.array([1,0,0])))))**2)*(1/omega)*(1+n)
            I[m] += (np.linalg.norm(np.array([0,0,1]).dot(A_mat_z[m].dot(Full_Rot_yz(y2,z2,ang_div,np.array([1,0,0])))))**2)*(1/omega)*(1+n)

    I[m] /= 3e10
    I[m] /= counter
#plt.plot(x,lorentzian(x,I),label='Single Incident Polarization, Average Scattered')

I = np.zeros(mode_num)

ang_div = 12

counter_total = 0

for m in range(mode_num):
    n = 1/(np.exp(hc*freqs[m]/kT) - 1)
    omega = 2*np.pi*29979245800*freqs[m]
    I[m] += (np.linalg.norm(np.array([1,0,0]).dot(A_mat_x[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([0,1,0]).dot(A_mat_y[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([0,0,1]).dot(A_mat_z[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)

    I[m] /= 3e10
#plt.plot(x,lorentzian(x,I),label='Single Incident Polarization, Scattered (1,0,0)')


I = np.zeros(mode_num)

ang_div = 12

counter_total = 0

for m in range(mode_num):
    n = 1/(np.exp(hc*freqs[m]/kT) - 1)
    omega = 2*np.pi*29979245800*freqs[m]
    I[m] += (np.linalg.norm(np.array([1,0,0]).dot(A_mat_x[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([0,1,0]).dot(A_mat_y[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([0,0,1]).dot(A_mat_z[m].dot(np.array([1,0,0]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([1,0,0]).dot(A_mat_x[m].dot(np.array([0,1,0]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([0,1,0]).dot(A_mat_y[m].dot(np.array([0,1,0]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([0,0,1]).dot(A_mat_z[m].dot(np.array([0,1,0]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([1,0,0]).dot(A_mat_x[m].dot(np.array([0,0,1]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([0,1,0]).dot(A_mat_y[m].dot(np.array([0,0,1]))))**2)*(1/omega)*(1+n)
    I[m] += (np.linalg.norm(np.array([0,0,1]).dot(A_mat_z[m].dot(np.array([0,0,1]))))**2)*(1/omega)*(1+n)
    I[m] /= 9e10
plt.plot(x,lorentzian(x,I),label='Single Incident Polarization, Scattered (1,0,0)')


#I = np.zeros(mode_num)

#ang_div = 12

#counter_total = 0

#for m in range(mode_num):
#    n = 1/(np.exp(hc*freqs[m]/kT) - 1)
#    omega = 2*np.pi*29979245800*freqs[m]
#    counter = 0
#    for y1 in range(ang_div):
#        for z1 in range(ang_div):
#            for y2 in range(ang_div):
#                for z2 in range(ang_div):
#                    counter += 1/3
#                    counter_total += 1/3
#                    if counter%(ang_div**2) == 0:
#                        print(m,100*(counter/(ang_div**4)),100*(counter_total/(27*ang_div**4)))
#                    I[m] += (np.linalg.norm(np.array([1,0,0]).dot(0.5*(A_mat_xz[m]+A_mat_xy).dot(Full_Rot_yz(y2,z2,ang_div,np.array([1,0,0])))))**2)*(1/omega)*(1+n)
#    for x1 in range(ang_div):
#        for z1 in range(ang_div):
#            for x2 in range(ang_div):
#                for z2 in range(ang_div):
#                    counter += 1/3
#                    counter_total += 1/3
#                    if counter%(ang_div**2) == 0:
#                        print(m,100*(counter/(ang_div**4)),100*(counter_total/(27*ang_div**4)))
#                    I[m] += (np.linalg.norm(np.array([0,1,0]).dot(0.5*(A_mat_yz[m]+A_mat_yx).dot(Full_Rot_zx(z2,x2,ang_div,np.array([0,1,0])))))**2)*(1/omega)*(1+n)
#    for x1 in range(ang_div):
#        for y1 in range(ang_div):
#            for x2 in range(ang_div):
#                for y2 in range(ang_div):
#                    counter += 1/3
#                    counter_total += 1/3
#                    if counter%(ang_div**2) == 0:
#                        print(m,100*(counter/(ang_div**4)),100*(counter_total/(27*ang_div**4)))
#                    I[m] += (np.linalg.norm(np.array([0,0,1]).dot(0.5*(A_mat_zx[m]+A_mat_zy).dot(Full_Rot_xy(x2,y2,ang_div,np.array([0,0,1])))))**2)*(1/omega)*(1+n)

#    I[m] /= 1e10
#plt.plot(x,lorentzian(x,I),label='Average Incident, Average Scattered')

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
plt.show()

