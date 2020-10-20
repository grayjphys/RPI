import numpy as np
import matplotlib.pyplot as plt

first_n_bands = 22
num_kpoints = 0
num_bands = 0
num_ions = 0
ecounter = 0
vcounter = 0
norm = 0
tmp_dft_path = []
tmp_dft_band = []
tmp_dft_bands = []
tmp_dft_eigenvector = []
tmp_dft_eigenvectors = []
f = open("DAT")
for counter, line in enumerate(f):
	if counter == 0:
		num_kpoints =  int(line.split()[0])
	elif counter == 1:
		num_bands =  int(line.split()[0])
	elif counter == 2:
		num_ions =  int(line.split()[0])
	elif abs(counter-3) < num_kpoints:
		tmp_dft_path.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
	elif abs(counter-3) >= num_kpoints and abs(counter-3) < num_kpoints*(num_bands+1):
		if ecounter%num_bands == 0 and ecounter != 0:
			tmp_dft_bands.append(np.asarray(tmp_dft_band))
			tmp_dft_band = []
		en = float(line.split()[0])
		if ecounter%num_bands < first_n_bands:
			tmp_dft_band.append(en)
		ecounter+=1
	elif abs(counter-3) >= num_kpoints*(num_bands+1):
		if vcounter%num_bands == 0 and vcounter != 0:
			tmp_dft_eigenvectors.append(np.asarray(tmp_dft_eigenvector))
			tmp_dft_eigenvector = []
		s = float(line.split()[0])*0.5
		py = float(line.split()[1])*0.5
		pz = float(line.split()[2])*0.5
		px = float(line.split()[3])*0.5
		dxy = float(line.split()[4])*0.5
		dyz = float(line.split()[5])*0.5
		dz2 = float(line.split()[6])*0.5
		dzx = float(line.split()[7])*0.5
		dx2_y2 = float(line.split()[8])*0.5
		if vcounter%num_bands < first_n_bands:
			eigvec = np.array([s,s,px,py,pz,px,py,pz,dxy,dyz,dzx,dx2_y2,dz2,dxy,dyz,dzx,dx2_y2,dz2])
			norm = np.linalg.norm(eigvec)
			if norm == 0:
				norm = 1
			tmp_dft_eigenvector.append(eigvec/norm)
		vcounter+=1

tmp_dft_bands.append(np.asarray(tmp_dft_band))
tmp_dft_bands = np.asarray(tmp_dft_bands)
tmp_dft_eigenvectors.append(np.asarray(tmp_dft_eigenvector))
tmp_dft_bands = np.asarray(tmp_dft_bands)


###########################################################################################################
#
# Define region to print band structure
#
###########################################################################################################

dft_path = np.delete(tmp_dft_path,np.s_[:-1:100],0)
dft_path = np.insert(dft_path,0,tmp_dft_path[0],axis=0)
dft_bands= np.delete(tmp_dft_bands,np.s_[:-1:100],0)
dft_bands = np.insert(dft_bands,0,tmp_dft_bands[0],axis=0)
dft_eigenvectors= np.delete(tmp_dft_eigenvectors,np.s_[:-1:100],0)
dft_eigenvectors = np.insert(dft_eigenvectors,0,tmp_dft_eigenvectors[0],axis=0)

del tmp_dft_path
del tmp_dft_bands
del tmp_dft_eigenvectors
	
ticks = [0]
tick = 0
# print(0, dft_path[0], ticks[0])
for i in range(1,len(dft_path)):
	d = np.linalg.norm(dft_path[i]-dft_path[i-1])
	if d < 0.3:
		tick += d
	else:
		tick += 0.000001
	ticks.append(tick)
	# print(i, dft_path[i], ticks[i])
fig, ax = plt.subplots()

xcoords = [0,0.8660254037844393,1.5731321849709863, 1.8793184028195633,2.7885996831564928]
for xc in xcoords:
    plt.axvline(x=xc,c='grey')

for i,band in enumerate(dft_bands.T):
	ax.plot(ticks,band,c = 'k')
	if i == 0:
		l0, = ax.plot(ticks,band,c='k')

plt.xticks(xcoords,(r'$L$',r'$\Gamma$',r'$X$',r'$U,K$',r'$\Gamma$'),fontsize=20)
plt.yticks(fontsize=20)
plt.ylim((-7,15))
plt.ylabel('Energy (eV)',fontsize=30)
plt.title('Si Band Structure',fontsize=30)
plt.show()