import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pickle
from scipy.fftpack import fft, fftfreq
from scipy.interpolate import interp1d
import os

mol_name = "Pyridine"
mol_name_caps = "PYRIDINE"

Carb = 12.011  # mass in amu
Nit = 14.007  # mass in amu
Hel = 1.008  # mass in amu

masses = [Carb,Carb,Carb,Carb,Carb,Nit,Hel,Hel,Hel,Hel,Hel]

skip = []
# skip = [21,20]

h = 0.0001

N = 2000  # number of timesteps
T = 0.02418884326509  # time interval in femtoseconds
M = 27  # number of modes

gamma = 0.1

planck = 4.135667696E-15  # eV * s
c = 2.99792458E10  # cm / s
KbT = 0.25851999786  # eV
wavelength_inc = 514  # nm
wavenumber_inc = 1e7/wavelength_inc  # 1 / cm

#######################################################################################################
###                                    GET POLARIZABILITY TENSOR                                    ###
#######################################################################################################

p = Path('.')

dirs = [x for x in p.iterdir() if x.is_dir() and "POL" in x.name]
dirs = np.sort(dirs)

pol_plus = np.zeros((M,3,3),dtype=complex)
pol_minus = np.zeros((M,3,3),dtype=complex)
files_found = 0

for progress, d in enumerate(dirs):
    print(d.name)
    with open(d.name + "/POL_TENSOR_PLUS.pkl",'rb') as f:
        pol_plus += np.array(pickle.load(f),dtype=complex)
    with open(d.name + "/POL_TENSOR_MINUS.pkl",'rb') as f:
        pol_minus += np.array(pickle.load(f),dtype=complex)

    files_found += 1

pol_plus /= files_found
pol_minus /= files_found

#######################################################################################################
###                                        GET FREQUENCIES                                          ###
#######################################################################################################

os.system(f"grep -A {M} FREQ " + mol_name + f"-VIBRATIONS-1.mol | tail -{M} > FREQ.dat")

freqs = np.zeros(M)
wavenumbers = np.zeros(M)
bose_einstein_term = np.zeros(M)

with open("FREQ.dat") as f:
    for m,line in enumerate(f):
        wavenumbers[m] = float(line)
        freqs[m] = float(line) * 2.99792458E-5 #* 2 * np.pi
        bose_einstein_term[m] = 1 / (1 - np.exp((-planck * c * wavenumbers[m]) / KbT))

w_values = fftfreq(N,d=T)[:N // 2]

indices = []
for i,w in enumerate(w_values):
    if min(freqs) <= w <= max(freqs):
        indices.append(i)
indices.append(indices[-1]+1)
indices = [indices[0] - 1] + indices

########################################################################################################
###                                      CALCULATE DERIVATIVES                                       ###
########################################################################################################

alpha_p_squared = np.zeros(M)
gamma_p_squared = np.zeros(M)
derivative_term = np.zeros(M)

for m in range(M):
    if m not in skip:
        alpha_p_squared[m] = (1 / 3) * np.linalg.norm((pol_plus[m][0][0] - pol_minus[m][0][0] + 
                            pol_plus[m][1][1] - pol_minus[m][1][1] + 
                            pol_plus[m][2][2] - pol_minus[m][2][2]) / (2 * h)) ** 2

        gamma_p_squared[m] = 0.5 * (np.linalg.norm((pol_plus[m][0][0] - pol_minus[m][0][0] - pol_plus[m][1][1] + pol_minus[m][1][1]) / (2 * h)) ** 2 +
                            np.linalg.norm((pol_plus[m][0][0] - pol_minus[m][0][0] - pol_plus[m][2][2] + pol_minus[m][2][2]) / (2 * h)) ** 2 +
                            np.linalg.norm((pol_plus[m][1][1] - pol_minus[m][1][1] - pol_plus[m][2][2] + pol_minus[m][2][2]) / (2 * h)) ** 2 +
                            6 * (np.linalg.norm((pol_plus[m][0][0] - pol_minus[m][0][0]) / (2 * h)) ** 2 +
                                    np.linalg.norm((pol_plus[m][1][1] - pol_minus[m][1][1]) / (2 * h)) ** 2 +
                                    np.linalg.norm((pol_plus[m][2][2] - pol_minus[m][2][2]) / (2 * h)) ** 2))

        derivative_term[m] = (45 * alpha_p_squared[m] + 7 * gamma_p_squared[m]) / 45

########################################################################################################
###                                    CALCULATE RAMAN INTENSITY                                     ###
########################################################################################################

raman_intensity = np.zeros(M)

for m in range(M):
    raman_intensity[m] = (1/wavenumbers[m]) * ((wavenumber_inc - wavenumbers[m]) ** 4) * derivative_term[m] * bose_einstein_term[m] * 5.76929578E-54 # needed for correct units of cm^2

x = np.linspace(500,1700,1000)

def lorentzian(x,I):
    broadening = 20
    gamma_sqr= (0.5*broadening)**2
    L = 0
    for i,v in enumerate(wavenumbers):
        L += I[i]*gamma_sqr/((x-v)**2+gamma_sqr)
    return L


data = []
for val in x:
    data.append(lorentzian(val,raman_intensity))

plt.plot(x,data,label="Total Raman")
plt.legend(fontsize=15)
plt.xlabel(r'Wavenumber cm$^{-1}$',fontsize=25)
plt.ylabel('A.U.',fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
plt.show()
plt.savefig("Raman-800-11-11-11.png")
