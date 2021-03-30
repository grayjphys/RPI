import numpy as np
import pickle
from pathlib import Path
import sys
import os
from scipy.fftpack import fft, fftfreq
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

np.set_printoptions(threshold=sys.maxsize)

p = Path('.') # This program should be run in a POL-[]-PROP-[] directory as described in the README.md file
folder = os.getcwd()
data_name = folder[-12:]

dirs = [x for x in p.iterdir() if x.is_dir() and "MODE" in x.name] # Collects all directories of the form PYRIDINE-[]-MODE-[]
dirs = np.sort(dirs)

mol_name = "Pyridine"
mol_name_caps = "PYRIDINE"

Carb = 12.011  # mass of carbon in amu
Nit = 14.007  # mass of nitrogen  in amu
Hyd = 1.008  # mass of hydrogen  in amu

masses = [Carb,Carb,Carb,Carb,Carb,Nit,Hyd,Hyd,Hyd,Hyd,Hyd]

skip = []
# if you want to skip specific modes, the code will account for that. for example: skip = [21,20]

h = 0.0001 # the step size for derivatives

box_size = 11 # DFT box size in angstroms
box_size_bohr = box_size * 1.8897259886 # \\ in bohr

N = 2000  # number of timesteps
A = 11  # number of atoms in the molecule
T = 0.02418884326509  # time interval in femtoseconds for fourier transform
M = 27  # number of modes

gamma = 0.1 # required to give lifetime to electronic transitions

planck = 4.135667696E-15  # eV * s
c = 2.99792458E10  # cm / s
KbT = 0.25851999786  # eV
wavelength_inc = 514  # nm
wavenumber_inc = 1e7/wavelength_inc  # 1 / cm


########################################################################################################
###                        COLLECT ELECTRIC FIELD VALUES AT CENTER OF BOX                            ###
########################################################################################################

Efield_x_minus = []
Efield_x_plus = []
Efield_y_minus = []
Efield_y_plus = []
Efield_z_minus = []
Efield_z_plus = []

for progress, d in enumerate(dirs): # loop through vibrational modes, plus and minus

    ### Read electric field files and create new files. This makes it easier to find the center of mass for each vibrational mode. ###
    
    os.system("rm " + d.name + "/POS.dat " + d.name + "/VOX.dat")
    os.system("head -" + str(6+A) + " " + d.name + "/"+ mol_name +"-efield_x-1_1.cube | tail -" + str(A) + " | awk \'{print $3,$4,$5}\'>" + d.name + "/POS.dat")
    os.system("head -4 " + d.name + "/"+ mol_name +"-efield_x-1_1.cube | tail -1 | awk \'{print $1,$2}\'>>" + d.name+"/VOX.dat")
    os.system("head -5 " + d.name + "/"+ mol_name +"-efield_x-1_1.cube | tail -1 | awk \'{print $1,$3}\'>>" + d.name+"/VOX.dat")
    os.system("head -6 " + d.name + "/"+ mol_name +"-efield_x-1_1.cube | tail -1 | awk \'{print $1,$4}\'>>" + d.name+"/VOX.dat")

    vox = np.loadtxt(d.name + "/VOX.dat")
    vox_T = vox.T

    ### There is a slight discrepancy between DFT box size and the size obtained from the electric field files. The factors take this into account. ###
    factor_x = box_size_bohr/(vox_T[0][0]*vox_T[1][0])
    factor_y = box_size_bohr/(vox_T[0][1]*vox_T[1][1])
    factor_z = box_size_bohr/(vox_T[0][2]*vox_T[1][2])

    pos = np.loadtxt(d.name + "/POS.dat")
    pos_T = pos.T

    ### calculate center of mass for each mode. ###
    com = np.array([np.sum([masses[e]*i for e, i in enumerate(pos_T[0])]),
	  			    np.sum([masses[e]*i for e, i in enumerate(pos_T[1])]),
	  			    np.sum([masses[e]*i for e, i in enumerate(pos_T[2])])])/np.sum(masses)

#     print(os.getcwd())
#     print(d.name)
#     print(pos)
#     print(com)
    
    x_range = vox_T[1][0]*np.arange(vox_T[0][0])
    y_range = vox_T[1][1]*np.arange(vox_T[0][1])
    z_range = vox_T[1][2]*np.arange(vox_T[0][2])

    print("{0}%".format(50*progress/M))
    
    ### Find electric field values from the .cube files and interpolate the value at the center of mass. ###
    ### The electric field is a vector, so you need (x,y,z) components. There is a separate .cube file for each direction and time. ###
    
    efield_files_x = np.array([e for e in d.iterdir() if mol_name + "-efield_x-1_" in e.name])
    ex_nums = []
    order = []
    for e in efield_files_x:
        num = int(e.name[20:-5])
        ex_nums.append(num-1)
        order = np.argsort(ex_nums)
    efield_files_x = efield_files_x[order]
    efield_x_values_minus = []
	
    if "MINUS" in d.name:
        for file in efield_files_x:
            with open(file) as f:
                data = []
                for num, line in enumerate(f):
                    if num > A + 5:
                        for ef in line.split():
                            data.append(float(ef))
                data = np.array(data).reshape((int(vox_T[0][0]),int(vox_T[0][1]),int(vox_T[0][2])))
                interpolating_function = RegularGridInterpolator((x_range, y_range, z_range), data)
                efield_x_values_minus.append(interpolating_function(com)[0])
        Efield_x_minus.append(efield_x_values_minus)
    efield_x_values_plus = []

    if "PLUS" in d.name:
        for file in efield_files_x:
            with open(file) as f:
                data = []
                for num, line in enumerate(f):
                    if num > A + 5:
                        for ef in line.split():
                            data.append(float(ef))
                data = np.array(data).reshape((int(vox_T[0][0]),int(vox_T[0][1]),int(vox_T[0][2])))
                interpolating_function = RegularGridInterpolator((x_range, y_range, z_range), data)
                efield_x_values_plus.append(interpolating_function(com)[0])
        Efield_x_plus.append(efield_x_values_plus)

    efield_files_y = np.array([e for e in d.iterdir() if mol_name + "-efield_y-1_" in e.name])
    ey_nums = []
    order = []
    for e in efield_files_y:
        num = int(e.name[20:-5])
        ey_nums.append(num-1)
        order = np.argsort(ey_nums)
    efield_files_y = efield_files_y[order]
    efield_y_values_minus = []
	
    if "MINUS" in d.name:
        for file in efield_files_y:
            with open(file) as f:
                data = []
                for num, line in enumerate(f):
                    if num > A + 5:
                        for ef in line.split():
                            data.append(float(ef))
                data = np.array(data).reshape((int(vox_T[0][0]),int(vox_T[0][1]),int(vox_T[0][2])))
                interpolating_function = RegularGridInterpolator((x_range, y_range, z_range), data)
                efield_y_values_minus.append(interpolating_function(com)[0])
        Efield_y_minus.append(efield_y_values_minus)
    efield_y_values_plus = []

    if "PLUS" in d.name:
        for file in efield_files_y:
            with open(file) as f:
                data = []
                for num, line in enumerate(f):
                    if num > A + 5:
                        for ef in line.split():
                            data.append(float(ef))
                data = np.array(data).reshape((int(vox_T[0][0]),int(vox_T[0][1]),int(vox_T[0][2])))
                interpolating_function = RegularGridInterpolator((x_range, y_range, z_range), data)
                efield_y_values_plus.append(interpolating_function(com)[0])
        Efield_y_plus.append(efield_y_values_plus)

    efield_files_z = np.array([e for e in d.iterdir() if mol_name + "-efield_z-1_" in e.name])
    ez_nums = []
    order = []
    for e in efield_files_z:
        num = int(e.name[20:-5])
        ez_nums.append(num-1)
        order = np.argsort(ez_nums)
    efield_files_z = efield_files_z[order]
    efield_z_values_minus = []
	
    if "MINUS" in d.name:
        for file in efield_files_z:
            with open(file) as f:
                data = []
                for num, line in enumerate(f):
                    if num > A + 5:
                        for ef in line.split():
                            data.append(float(ef))
                data = np.array(data).reshape((int(vox_T[0][0]),int(vox_T[0][1]),int(vox_T[0][2])))
                interpolating_function = RegularGridInterpolator((x_range, y_range, z_range), data)
                efield_z_values_minus.append(interpolating_function(com)[0])
        Efield_z_minus.append(efield_z_values_minus)
    efield_z_values_plus = []

    if "PLUS" in d.name:
        for file in efield_files_z:
            with open(file) as f:
                data = []
                for num, line in enumerate(f):
                    if num > A + 5:
                        for ef in line.split():
                            data.append(float(ef))
                data = np.array(data).reshape((int(vox_T[0][0]),int(vox_T[0][1]),int(vox_T[0][2])))
                interpolating_function = RegularGridInterpolator((x_range, y_range, z_range), data)
                efield_z_values_plus.append(interpolating_function(com)[0])
        Efield_z_plus.append(efield_z_values_plus)

Efield_x_minus = np.array(Efield_x_minus)
Efield_x_plus =  np.array(Efield_x_plus)
Efield_y_minus = np.array(Efield_y_minus)
Efield_y_plus =  np.array(Efield_y_plus)
Efield_z_minus = np.array(Efield_z_minus)
Efield_z_plus =  np.array(Efield_z_plus)

Efield_x_minus_fft = np.zeros((M,N),dtype=complex)
Efield_x_plus_fft = np.zeros((M,N),dtype=complex)
Efield_y_minus_fft = np.zeros((M,N),dtype=complex)
Efield_y_plus_fft = np.zeros((M,N),dtype=complex)
Efield_z_minus_fft = np.zeros((M,N),dtype=complex)
Efield_z_plus_fft = np.zeros((M,N),dtype=complex)

### Fourier transform the fields to get the electric field in frequency space. ###
for i in range(M):
    data = Efield_x_minus[i]
    Efield_x_minus_fft[i] = fft(data)
    data = Efield_x_plus[i]
    Efield_x_plus_fft[i] = fft(data)
    data = Efield_y_minus[i]
    Efield_y_minus_fft[i] = fft(data)
    data = Efield_y_plus[i]
    Efield_y_plus_fft[i] = fft(data)
    data = Efield_z_minus[i]
    Efield_z_minus_fft[i] = fft(data)
    data = Efield_z_plus[i]
    Efield_z_plus_fft[i] = fft(data)

print("************************** EFIELD DONE ****************************")

### Dump the result into a pickle file. This is just in case something goes wrong later on. It will save a lot of time to just read the values in. ###
with open('EFIELD_X_MINUS_FFT.pkl','wb') as file:
    pickle.dump(Efield_x_minus_fft, file)

with open('EFIELD_X_PLUS_FFT.pkl','wb') as file:
    pickle.dump(Efield_x_plus_fft, file)

with open('EFIELD_Y_MINUS_FFT.pkl','wb') as file:
    pickle.dump(Efield_y_minus_fft, file)

with open('EFIELD_Y_PLUS_FFT.pkl','wb') as file:
    pickle.dump(Efield_y_plus_fft, file)

with open('EFIELD_Z_MINUS_FFT.pkl','wb') as file:
    pickle.dump(Efield_z_minus_fft, file)

with open('EFIELD_Z_PLUS_FFT.pkl','wb') as file:
    pickle.dump(Efield_z_plus_fft, file)

#######################################################################################################
###                        COLLECT DIPOLE MOMENT VALUES AT CENTER OF MASS                           ###
#######################################################################################################

Moment_x_minus = []
Moment_x_plus = []
Moment_y_minus = []
Moment_y_plus = []
Moment_z_minus = []
Moment_z_plus = []

for enum, d in enumerate(dirs): # loop over each vibrational mode, plus and minus
	
    print("{0}%".format(50*enum/M))

    moments_x = []
    moments_y = []
    moments_z = []
    
    ### Prepare a file with just the moment data from the .mol file. ###
    os.system("grep X= " + d.name+"/Vibration-" + mol_name +".out | awk \'{print $2,$4,$6}\'>" + d.name+"/MOM.dat")

    with open(d.name+"/MOM.dat") as f:
        for i,line in enumerate(f):
            t = i*T
            if i == 0:
                mom_x0 = float(line.split()[0])
                mom_y0 = float(line.split()[1])
                mom_z0 = float(line.split()[2])
            else:
                if len(moments_x) < N:
                    moments_x.append((float(line.split()[0]) - mom_x0) * np.exp(-gamma*t)) # gamma factor necessary to give a time decay for electronic transitions.
                    moments_y.append((float(line.split()[1]) - mom_y0) * np.exp(-gamma*t))
                    moments_z.append((float(line.split()[2]) - mom_z0) * np.exp(-gamma*t))
        # print(d.name,len(moments_x))
	
    if "MINUS" in d.name:
        Moment_x_minus.append(moments_x)
        Moment_y_minus.append(moments_y)
        Moment_z_minus.append(moments_z)
    elif "PLUS" in d.name:
        Moment_x_plus.append(moments_x)
        Moment_y_plus.append(moments_y)
        Moment_z_plus.append(moments_z)
    

Moment_x_minus = np.array(Moment_x_minus)
Moment_x_plus =  np.array(Moment_x_plus)
Moment_y_minus = np.array(Moment_y_minus)
Moment_y_plus =  np.array(Moment_y_plus)
Moment_z_minus = np.array(Moment_z_minus)
Moment_z_plus =  np.array(Moment_z_plus)

Moment_x_minus_fft = np.zeros((M,N),dtype=complex)
Moment_x_plus_fft = np.zeros((M,N),dtype=complex)
Moment_y_minus_fft = np.zeros((M,N),dtype=complex)
Moment_y_plus_fft = np.zeros((M,N),dtype=complex)
Moment_z_minus_fft = np.zeros((M,N),dtype=complex)
Moment_z_plus_fft = np.zeros((M,N),dtype=complex)

### Fourier transform dipole moments into frequency space. ###
for i in range(M):
    if i not in skip:
        data = Moment_x_minus[i]
        Moment_x_minus_fft[i] = fft(data)
        data = Moment_x_plus[i]
        Moment_x_plus_fft[i] = fft(data)
        data = Moment_y_minus[i]
        Moment_y_minus_fft[i] = fft(data)
        data = Moment_y_plus[i]
        Moment_y_plus_fft[i] = fft(data)
        data = Moment_z_minus[i]
        Moment_z_minus_fft[i] = fft(data)
        data = Moment_z_plus[i]
        Moment_z_plus_fft[i] = fft(data)

print("************************** DIPOLE DONE ****************************")

### Dump the result into a pickle file. This is just in case something goes wrong later on. It will save a lot of time to just read the values in. ###

with open('MOMENT_X_MINUS_FFT.pkl','wb') as file:
    pickle.dump(Moment_x_minus_fft, file)

with open('MOMENT_X_PLUS_FFT.pkl','wb') as file:
    pickle.dump(Moment_x_plus_fft, file)

with open('MOMENT_Y_MINUS_FFT.pkl','wb') as file:
    pickle.dump(Moment_y_minus_fft, file)

with open('MOMENT_Y_PLUS_FFT.pkl','wb') as file:
    pickle.dump(Moment_y_plus_fft, file)

with open('MOMENT_Z_MINUS_FFT.pkl','wb') as file:
    pickle.dump(Moment_z_minus_fft, file)

with open('MOMENT_Z_PLUS_FFT.pkl','wb') as file:
    pickle.dump(Moment_z_plus_fft, file)

 #######################################################################################################
 ###                                        GET FREQUENCIES                                          ###
 #######################################################################################################

### Find frequency (actually wavenumber) of vibrational modes from .mol file. units of cm^-1 ###

os.system("grep -A {0} FREQ ".format(M) + mol_name_caps + "-VIBRATIONS-1.mol | tail -{0} > FREQ.dat".format(M))

freqs = np.zeros(M)
wavenumbers = np.zeros(M)
bose_einstein_term = np.zeros(M)

with open("FREQ.dat") as f:
    for m,line in enumerate(f):
        wavenumbers[m] = float(line)
        freqs[m] = float(line) * 2.99792458E-5 * 2 * np.pi
        bose_einstein_term[m] = 1 / (1 - np.exp((-planck * c * wavenumbers[m]) / KbT)) # Phonons follow bose-einstein statistics.

w_values = fftfreq(N,T)[:N // 2]

indices = []
for i,w in enumerate(w_values):
    if min(freqs) <= w <= max(freqs):
        indices.append(i)
indices.append(indices[-1]+1)
indices = [indices[0] - 1] + indices

# print(w_values[indices])
# print(freqs)

 #######################################################################################################
 ###                               CALCULATE POLARIZABILITY TENSOR                                   ###
 #######################################################################################################

### Read in previously pickled efields and moments. Comment the pickle sections and edit appropriately if you feel the need. ###
with open('EFIELD_X_MINUS_FFT.pkl','rb') as file:
    Efield_x_minus_fft = pickle.load(file)

with open('EFIELD_X_PLUS_FFT.pkl','rb') as file:
    Efield_x_plus_fft = pickle.load(file)

with open('EFIELD_Y_MINUS_FFT.pkl','rb') as file:
    Efield_y_minus_fft = pickle.load(file)

with open('EFIELD_Y_PLUS_FFT.pkl','rb') as file:
    Efield_y_plus_fft = pickle.load(file)

with open('EFIELD_Z_MINUS_FFT.pkl','rb') as file:
    Efield_z_minus_fft = pickle.load(file)

with open('EFIELD_Z_PLUS_FFT.pkl','rb') as file:
    Efield_z_plus_fft = pickle.load(file)

with open('MOMENT_X_MINUS_FFT.pkl','rb') as file:
    Moment_x_minus_fft = pickle.load(file)

with open('MOMENT_X_PLUS_FFT.pkl','rb') as file:
    Moment_x_plus_fft = pickle.load(file)

with open('MOMENT_Y_MINUS_FFT.pkl','rb') as file:
    Moment_y_minus_fft = pickle.load(file)

with open('MOMENT_Y_PLUS_FFT.pkl','rb') as file:
    Moment_y_plus_fft = pickle.load(file)

with open('MOMENT_Z_MINUS_FFT.pkl','rb') as file:
    Moment_z_minus_fft = pickle.load(file)

with open('MOMENT_Z_PLUS_FFT.pkl','rb') as file:
    Moment_z_plus_fft = pickle.load(file)

polarizability_tensor_minus = np.zeros((M,9,N),dtype = complex) # 3x3 polarizability tensor for each mode, 
                                                                # the N is due to the fact that the FFT will give a range of frequencies
								# including the frequency of the mode. 

for m in range(M):
    if m not in skip: # remember the skip from earlier? :P 
        for w in indices:
            if np.linalg.norm(Efield_x_minus_fft[m][w]) > 0:
                polarizability_tensor_minus[m][0][w] = Moment_x_minus_fft[m][w] / Efield_x_minus_fft[m][w]
                polarizability_tensor_minus[m][3][w] = Moment_y_minus_fft[m][w] / Efield_x_minus_fft[m][w]
                polarizability_tensor_minus[m][6][w] = Moment_z_minus_fft[m][w] / Efield_x_minus_fft[m][w]
            if np.linalg.norm(Efield_y_minus_fft[m][w]) > 0:
                polarizability_tensor_minus[m][1][w] = Moment_x_minus_fft[m][w] / Efield_y_minus_fft[m][w]
                polarizability_tensor_minus[m][4][w] = Moment_y_minus_fft[m][w] / Efield_y_minus_fft[m][w]
                polarizability_tensor_minus[m][7][w] = Moment_z_minus_fft[m][w] / Efield_y_minus_fft[m][w]
            if np.linalg.norm(Efield_z_minus_fft[m][w]) > 0:
                polarizability_tensor_minus[m][2][w] = Moment_x_minus_fft[m][w] / Efield_z_minus_fft[m][w]
                polarizability_tensor_minus[m][5][w] = Moment_y_minus_fft[m][w] / Efield_z_minus_fft[m][w]
                polarizability_tensor_minus[m][8][w] = Moment_z_minus_fft[m][w] / Efield_z_minus_fft[m][w]

polarizability_tensor_minus_interp = np.zeros((M,3,3),dtype=complex) 

### interpolate the polarizability tensor at the frequency of the mode. ###
for m in range(M):
    if m not in skip:
        polarizability_tensor_minus_interp[m][0][0] = interp1d(w_values[indices],polarizability_tensor_minus[m][0][indices])(freqs[m])
        polarizability_tensor_minus_interp[m][0][1] = interp1d(w_values[indices],polarizability_tensor_minus[m][1][indices])(freqs[m])
        polarizability_tensor_minus_interp[m][0][2] = interp1d(w_values[indices],polarizability_tensor_minus[m][2][indices])(freqs[m])
        polarizability_tensor_minus_interp[m][1][0] = interp1d(w_values[indices],polarizability_tensor_minus[m][3][indices])(freqs[m])
        polarizability_tensor_minus_interp[m][1][1] = interp1d(w_values[indices],polarizability_tensor_minus[m][4][indices])(freqs[m])
        polarizability_tensor_minus_interp[m][1][2] = interp1d(w_values[indices],polarizability_tensor_minus[m][5][indices])(freqs[m])
        polarizability_tensor_minus_interp[m][2][0] = interp1d(w_values[indices],polarizability_tensor_minus[m][6][indices])(freqs[m])
        polarizability_tensor_minus_interp[m][2][1] = interp1d(w_values[indices],polarizability_tensor_minus[m][7][indices])(freqs[m])
        polarizability_tensor_minus_interp[m][2][2] = interp1d(w_values[indices],polarizability_tensor_minus[m][8][indices])(freqs[m])

with open('POL_TENSOR_MINUS.pkl','wb') as file:
    pickle.dump(polarizability_tensor_minus_interp, file)

polarizability_tensor_plus= np.zeros((M,9,N),dtype=complex)

for m in range(M):
    if m not in skip:
        for w in indices:
            if np.linalg.norm(Efield_x_plus_fft[m][w]) > 0:
                polarizability_tensor_plus[m][0][w] = Moment_x_plus_fft[m][w] / Efield_x_plus_fft[m][w]
                polarizability_tensor_plus[m][3][w] = Moment_y_plus_fft[m][w] / Efield_x_plus_fft[m][w]
                polarizability_tensor_plus[m][6][w] = Moment_z_plus_fft[m][w] / Efield_x_plus_fft[m][w]
            if np.linalg.norm(Efield_y_plus_fft[m][w]) > 0:
                polarizability_tensor_plus[m][1][w] = Moment_x_plus_fft[m][w] / Efield_y_plus_fft[m][w]
                polarizability_tensor_plus[m][4][w] = Moment_y_plus_fft[m][w] / Efield_y_plus_fft[m][w]
                polarizability_tensor_plus[m][7][w] = Moment_z_plus_fft[m][w] / Efield_y_plus_fft[m][w]
            if np.linalg.norm(Efield_z_plus_fft[m][w]) > 0:
                polarizability_tensor_plus[m][2][w] = Moment_x_plus_fft[m][w] / Efield_z_plus_fft[m][w]
                polarizability_tensor_plus[m][5][w] = Moment_y_plus_fft[m][w] / Efield_z_plus_fft[m][w]
                polarizability_tensor_plus[m][8][w] = Moment_z_plus_fft[m][w] / Efield_z_plus_fft[m][w]

polarizability_tensor_plus_interp = np.zeros((M,3,3),dtype=complex)

### interpolate the polarizability tensor at the frequency of the mode. ###
for m in range(M):
    if m not in skip:
        polarizability_tensor_plus_interp[m][0][0] = interp1d(w_values[indices],polarizability_tensor_plus[m][0][indices])(freqs[m])
        polarizability_tensor_plus_interp[m][0][1] = interp1d(w_values[indices],polarizability_tensor_plus[m][1][indices])(freqs[m])
        polarizability_tensor_plus_interp[m][0][2] = interp1d(w_values[indices],polarizability_tensor_plus[m][2][indices])(freqs[m])
        polarizability_tensor_plus_interp[m][1][0] = interp1d(w_values[indices],polarizability_tensor_plus[m][3][indices])(freqs[m])
        polarizability_tensor_plus_interp[m][1][1] = interp1d(w_values[indices],polarizability_tensor_plus[m][4][indices])(freqs[m])
        polarizability_tensor_plus_interp[m][1][2] = interp1d(w_values[indices],polarizability_tensor_plus[m][5][indices])(freqs[m])
        polarizability_tensor_plus_interp[m][2][0] = interp1d(w_values[indices],polarizability_tensor_plus[m][6][indices])(freqs[m])
        polarizability_tensor_plus_interp[m][2][1] = interp1d(w_values[indices],polarizability_tensor_plus[m][7][indices])(freqs[m])
        polarizability_tensor_plus_interp[m][2][2] = interp1d(w_values[indices],polarizability_tensor_plus[m][8][indices])(freqs[m])

with open('POL_TENSOR_PLUS.pkl','wb') as file:
    pickle.dump(polarizability_tensor_plus_interp, file)


########################################################################################################
###                 CALCULATE DERIVATIVES (See mathematical formulae section of README.md)           ###
########################################################################################################

with open('POL_TENSOR_MINUS.pkl' , 'rb') as file:
    polarizability_tensor_minus = pickle.load(file)

with open('POL_TENSOR_PLUS.pkl' , 'rb') as file:
    polarizability_tensor_plus = pickle.load(file)

alpha_p_squared = np.zeros(M,dtype=float)
gamma_p_squared = np.zeros(M,dtype=float)
derivative_term = np.zeros(M,dtype=float)

for m in range(M):
    if m not in skip:
        alpha_p_squared[m] = np.linalg.norm((1 / 3) * (polarizability_tensor_plus[m][0][0] - polarizability_tensor_minus[m][0][0] + 
                            polarizability_tensor_plus[m][1][1] - polarizability_tensor_minus[m][1][1] + 
                            polarizability_tensor_plus[m][2][2] - polarizability_tensor_minus[m][2][2]) / (2 * h)) ** 2

        gamma_p_squared[m] = 0.5 * (np.linalg.norm((polarizability_tensor_plus[m][0][0] - polarizability_tensor_minus[m][0][0] - polarizability_tensor_plus[m][1][1] + polarizability_tensor_minus[m][1][1]) / (2 * h)) ** 2 +
                            np.linalg.norm((polarizability_tensor_plus[m][0][0] - polarizability_tensor_minus[m][0][0] - polarizability_tensor_plus[m][2][2] + polarizability_tensor_minus[m][2][2]) / (2 * h)) ** 2 +
                            np.linalg.norm((polarizability_tensor_plus[m][1][1] - polarizability_tensor_minus[m][1][1] - polarizability_tensor_plus[m][2][2] + polarizability_tensor_minus[m][2][2]) / (2 * h)) ** 2 +
                            6 * (np.linalg.norm((polarizability_tensor_plus[m][0][0] - polarizability_tensor_minus[m][0][0]) / (2 * h)) ** 2 +
                                    np.linalg.norm((polarizability_tensor_plus[m][1][1] - polarizability_tensor_minus[m][1][1]) / (2 * h)) ** 2 +
                                    np.linalg.norm((polarizability_tensor_plus[m][2][2] - polarizability_tensor_minus[m][2][2]) / (2 * h)) ** 2))

        derivative_term[m] =  (45 * alpha_p_squared[m] + 7 * gamma_p_squared[m]) / 45

########################################################################################################
###                                    CALCULATE RAMAN INTENSITY                                     ###
########################################################################################################

raman_intensity = np.zeros(M)

for m in range(M):
    raman_intensity[m] = ((wavenumber_inc - wavenumbers[m]) ** 4)/ wavenumbers[m] * derivative_term[m] * bose_einstein_term[m] * 5.76929578E-54 # needed for correct units

x = np.linspace(500,1700,1000)

def lorentzian(x,I): # Create spectrum with a broadening factor. This gives a more realistic spectrum because of quantum fluctuations.
    broadening = 20
    gamma_sqr= (0.5*broadening)**2
    L = 0
    for i,v in enumerate(wavenumbers):
        L += I[i]*gamma_sqr/((x-v)**2+gamma_sqr)
    return L


data = []
for val in x:
    data.append(lorentzian(val,raman_intensity))

with open(folder + "/" + data_name + "_SPECTRUM.dat",'w') as f:
    for e,d in enumerate(data):
        f.write("{0}\t{1}\n".format(x[e],d))

x = []
data = []
with open(folder + "/" + data_name + "_SPECTRUM.dat",'r') as f:
    for line in f:
        x.append(float(line.split()[0]))
        data.append(float(line.split()[1]))

### Plot the spectrum ###
plt.plot(x,data,label=data_name)
plt.legend(fontsize=15)
plt.xlabel(r'Wavenumber cm$^{-1}$',fontsize=25)
plt.ylabel('A.U.',fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
plt.show(block=False)
plt.pause(10)
plt.savefig("Raman-800-11-11-11.png")

