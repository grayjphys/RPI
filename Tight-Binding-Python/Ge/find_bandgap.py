import numpy as np
import os

Ef = float(os.popen("grep E-fermi OUTCAR").read().split()[2])
# print("e-fermi = ", Ef)

kpoints= []
energies = []
energy = []
with open("EIGENVAL") as f:
	for i,line in enumerate(f):
		if i <6:
			if i == 5:
				k_num = int(line.split()[1])
				band_num = int(line.split()[2])
			continue

		if (i-7)%(band_num+2) == 0:
			kpoints.append(np.asarray([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])]))
			energy = []

		for j in range(1,band_num+1):
			if (i-7) % (band_num+2)== j:
				energy.append(float(line.split()[1]))

		if (i-7)%(band_num+2) == (band_num+1):
			energies.append(energy)

kpoints = np.asarray(kpoints)
energies.append(energy)	
energies = energies[1:]

VB = []
CB = []
for i in range(len(energies)):
	for j in range(band_num):
		if energies[i][j] < Ef:
			VB.append(energies[i][j])
		if energies[i][j] > Ef:
			CB.append(energies[i][j])

VBM = max(VB)
CBM = min(CB)
Gap = CBM-VBM

print("VBM: " + str(VBM) + 'eV')
print("at: k = (" + str(kpoints[VB.index(VBM)][0]) + ', ' + str(kpoints[VB.index(VBM)][1]) + ', ' + str(kpoints[VB.index(VBM)][2]) + ')')
print("CBM: " + str(CBM) + 'eV')
print("at: k = (" + str(kpoints[CB.index(CBM)][0]) + ', ' + str(kpoints[CB.index(CBM)][1]) + ', ' + str(kpoints[CB.index(CBM)][2]) + ')')
print("Gap: " + str(Gap) + 'eV')
