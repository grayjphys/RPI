import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import minimize

a = 5.468721419 #Angstroms
dft_vbm = 5.4556495

def band_eigenvalues_and_eigenvectors(params,kx,ky,kz,dft_vbm):
	es = params[0]
	ep = params[1]
	ed = params[2]
	ed_ = params[3]
	vss = params[4]
	vxx = params[5]
	vxy = params[6]
	vsp = params[7]
	vsdS = params[8]
	vpdS = params[9]
	vpdP = params[10]
	vddS = params[11]
	vddP = params[12]
	vddD = params[13]
	v15 = vsdS/(3.0**0.5)
	v18 = 0 # Not sure if this is correct yet
	v19 = 0 # Not sure if this is correct yet
	v25 = ((3**0.5)*vpdS + vpdP)/(3.0*(3.0**0.5))
	v26 = ((3**0.5)*vpdS - 2.0*vpdP)/(3.0*(3.0**0.5))
	v27 = v26
	v28 = vpdP/(3.0**0.5)
	v29 = -1.0*vpdP/3.0
	v38 = -1.0*v28
	v39 = v29
	v48 = 0.0
	v49 = -2.0*v29
	v55 = (3.0*vddS+2.0*vddP+4.0*vddD)/9.0
	v56 = (3.0*vddS-1.0*vddP-2.0*vddD)/9.0
	v57 = v56
	v58 = 0.0
	v78 = (vddP-vddD)/3.0
	v59 = -2.0*v78/(3.0**0.5)
	v68 = -1.0*v78
	v79 = (vddP-vddD)/(3.0*(3.0**0.5))
	v69 = v79
	v88 = 2.0*(vddP/3.0)+(vddD/3.0)
	v89 = 0.0
	v99 = v88

	g = 0.25*np.array([
		complex(np.cos(0.5*np.pi*kx+0.5*np.pi*ky+0.5*np.pi*kz)+np.cos(0.5*np.pi*kx-0.5*np.pi*ky-0.5*np.pi*kz)+np.cos(-0.5*np.pi*kx+0.5*np.pi*ky-0.5*np.pi*kz)+np.cos(-0.5*np.pi*kx-0.5*np.pi*ky+0.5*np.pi*kz)
    ,np.sin(0.5*np.pi*kx+0.5*np.pi*ky+0.5*np.pi*kz)+np.sin(0.5*np.pi*kx-0.5*np.pi*ky-0.5*np.pi*kz)+np.sin(-0.5*np.pi*kx+0.5*np.pi*ky-0.5*np.pi*kz)+np.sin(-0.5*np.pi*kx-0.5*np.pi*ky+0.5*np.pi*kz)),
		
		complex(np.cos(0.5*np.pi*kx+0.5*np.pi*ky+0.5*np.pi*kz)+np.cos(0.5*np.pi*kx-0.5*np.pi*ky-0.5*np.pi*kz)-np.cos(-0.5*np.pi*kx+0.5*np.pi*ky-0.5*np.pi*kz)-np.cos(-0.5*np.pi*kx-0.5*np.pi*ky+0.5*np.pi*kz)
    ,np.sin(0.5*np.pi*kx+0.5*np.pi*ky+0.5*np.pi*kz)+np.sin(0.5*np.pi*kx-0.5*np.pi*ky-0.5*np.pi*kz)-np.sin(-0.5*np.pi*kx+0.5*np.pi*ky-0.5*np.pi*kz)-np.sin(-0.5*np.pi*kx-0.5*np.pi*ky+0.5*np.pi*kz)),
		
		complex(np.cos(0.5*np.pi*kx+0.5*np.pi*ky+0.5*np.pi*kz)-np.cos(0.5*np.pi*kx-0.5*np.pi*ky-0.5*np.pi*kz)+np.cos(-0.5*np.pi*kx+0.5*np.pi*ky-0.5*np.pi*kz)-np.cos(-0.5*np.pi*kx-0.5*np.pi*ky+0.5*np.pi*kz)
    ,np.sin(0.5*np.pi*kx+0.5*np.pi*ky+0.5*np.pi*kz)-np.sin(0.5*np.pi*kx-0.5*np.pi*ky-0.5*np.pi*kz)+np.sin(-0.5*np.pi*kx+0.5*np.pi*ky-0.5*np.pi*kz)-np.sin(-0.5*np.pi*kx-0.5*np.pi*ky+0.5*np.pi*kz)),
		
		complex(np.cos(0.5*np.pi*kx+0.5*np.pi*ky+0.5*np.pi*kz)-np.cos(0.5*np.pi*kx-0.5*np.pi*ky-0.5*np.pi*kz)-np.cos(-0.5*np.pi*kx+0.5*np.pi*ky-0.5*np.pi*kz)+np.cos(-0.5*np.pi*kx-0.5*np.pi*ky+0.5*np.pi*kz)
    ,np.sin(0.5*np.pi*kx+0.5*np.pi*ky+0.5*np.pi*kz)-np.sin(0.5*np.pi*kx-0.5*np.pi*ky-0.5*np.pi*kz)-np.sin(-0.5*np.pi*kx+0.5*np.pi*ky-0.5*np.pi*kz)+np.sin(-0.5*np.pi*kx-0.5*np.pi*ky+0.5*np.pi*kz))
		],dtype=complex)

	gc = np.conj(g)

	hamiltonian = np.array([
		[es,	     vss*g[0],  0.0,		0.0,		0.0,		vsp*g[1],  vsp*g[2],  vsp*g[3],  0.0,        0.0,        0.0,        0.0,        0.0,       v15*g[3], v15*g[1], v15*g[2],  v18*g[0], v19*g[0]],
		[vss*gc[0],  es,	   -vsp*gc[1], -vsp*gc[2], -vsp*gc[3],  0.0,	   0.0,		  0.0,       v15*gc[3],  v15*gc[1],  v15*gc[2],  v18*gc[0],  v19*gc[0], 0.0,      0.0,      0.0,       0.0,      0.0],
		[0.0,	    -vsp*g[1],  ep,		    0.0,		0.0,		vxx*g[0],  vxy*g[3],  vxy*g[2],  0.0,        0.0,        0.0,        0.0,        0.0,       v25*g[2], v26*g[0], v27*g[3],  v28*g[1], v29*g[1]],
		[0.0,       -vsp*g[2],  0.0,        ep,         0.0,        vxy*g[3],  vxx*g[0],  vxy*g[1],  0.0,        0.0,        0.0,        0.0,        0.0,       v27*g[1], v25*g[3], v26*g[0], -v28*g[2], v29*g[2]],
		[0.0,       -vsp*g[3],  0.0,        0.0,        ep,         vxy*g[2],  vxy*g[1],  vxx*g[0],  0.0,        0.0,        0.0,        0.0,        0.0,       v26*g[0], v27*g[2], v25*g[1],  v48*g[3], v49*g[3]],
		[gc[1]*vsp,  0.0,       vxx*gc[0],  vxy*gc[3],  vxy*gc[2],  ep,        0.0,       0.0,      -v25*gc[2], -v26*gc[0], -v27*gc[3], -v28*gc[1], -v29*gc[1], 0.0,      0.0,      0.0,       0.0,      0.0],
		[gc[2]*vsp,  0.0,       vxy*gc[3],  vxx*gc[0],  vxy*gc[1],  0.0,       ep,        0.0,      -v27*gc[1], -v25*gc[3], -v26*gc[0],  v28*gc[2], -v29*gc[2], 0.0,      0.0,      0.0,       0.0,      0.0],
		[gc[3]*vsp,  0.0,       vxy*gc[2],  vxy*gc[1],  vxx*gc[0],  0.0,       0.0,       ep,       -v26*gc[0], -v27*gc[2], -v25*gc[1], -v48*gc[3], -v49*gc[3], 0.0,      0.0,      0.0,       0.0,      0.0],
		[0.0,        v15*g[3],  0.0,        0.0,        0.0,       -v25*g[2], -v27*g[1], -v26*g[0],  ed,         0.0,        0.0,        0.0,        0.0,       v55*g[0], v56*g[2], v57*g[1],  v58*g[3], v59*g[3]],
		[0.0,        v15*g[1],  0.0,        0.0,        0.0,       -v26*g[0], -v25*g[3], -v27*g[2],  0.0,        ed,         0.0,        0.0,        0.0,       v56*g[2], v55*g[0], v56*g[3],  v68*g[1], v69*g[1]],
		[0.0,        v15*g[2],  0.0,        0.0,        0.0,       -v27*g[3], -v26*g[0], -v25*g[1],  0.0,        0.0,        ed,         0.0,        0.0,       v57*g[1], v57*g[3], v55*g[0],  v78*g[2], v79*g[2]],
		[0.0,        v18*g[0],  0.0,        0.0,        0.0,       -v28*g[1],  v28*g[2], -v48*g[3],  0.0,        0.0,        0.0,        ed_,        0.0,       v58*g[3], v68*g[1], v78*g[2],  v88*g[0], v89*g[0]],
		[0.0,        v19*g[0],  0.0,        0.0,        0.0,       -v29*g[1], -v29*g[2], -v49*g[3],  0.0,        0.0,        0.0,        0.0,        ed_,       v59*g[3], v69*g[1], v79*g[2],  v89*g[0], v99*g[0]],
		[v15*gc[3],  0.0,       v25*gc[2],  v27*gc[1],  v26*gc[0],  0.0,       0.0,       0.0,       v55*gc[0],  v56*gc[2],  v57*gc[1],  v58*gc[3],  v59*gc[3], ed,       0.0,      0.0,       0.0,      0.0],
		[v15*gc[1],  0.0,       v26*gc[0],  v25*gc[3],  v27*gc[2],  0.0,       0.0,       0.0,       v56*gc[2],  v55*gc[0],  v56*gc[3],  v68*gc[1],  v69*gc[1], 0.0,      ed,       0.0,       0.0,      0.0],
		[v15*gc[2],  0.0,       v27*gc[3],  v26*gc[0],  v25*gc[1],  0.0,       0.0,       0.0,       v57*gc[1],  v57*gc[3],  v55*gc[0],  v78*gc[2],  v79*gc[2], 0.0,      0.0,      ed,        0.0,      0.0],
		[v18*gc[0],  0.0,       v28*gc[1],  v28*gc[2],  v48*gc[3],  0.0,       0.0,       0.0,       v58*gc[3],  v68*gc[1],  v78*gc[2],  v88*gc[0],  v89*gc[0], 0.0,      0.0,      0.0,       ed_,      0.0],
		[v19*gc[0],  0.0,       v29*gc[1],  v29*gc[2],  v49*gc[3],  0.0,       0.0,       0.0,       v59*gc[3],  v69*gc[1],  v79*gc[2],  v89*gc[0],  v99*gc[0], 0.0,      0.0,      0.0,       0.0,      ed_]],dtype=complex)

	e, v = np.linalg.eigh(hamiltonian)
	VBM = len(e)*[dft_vbm]
	return e+VBM, v.T

def residual(params, dft_path, dft_bands, dft_eigenvectors):
	res = 0
	for p,k in enumerate(dft_path):
		tb_e, tb_v = band_eigenvalues_and_eigenvectors(params,k[0],k[1],k[2],dft_vbm)
		dft_e, dft_v = dft_bands[p],dft_eigenvectors[p]
		for i in range(len(tb_e)):
			for j in range(len(dft_e)):
				res += ((tb_e[i] - dft_e[j])**2)*np.linalg.norm(tb_v[i][::2].dot(dft_v[j]))
	return res


###########################################################################################################
#
# Get Data from DFT
#
###########################################################################################################

first_n_bands = 18
num_kpoints = 0
num_bands = 0
num_ions = 0
ecounter = 0
vcounter = 0
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
		s = float(line.split()[0])
		py = float(line.split()[1])
		pz = float(line.split()[2])
		px = float(line.split()[3])
		dxy = float(line.split()[4])
		dyz = float(line.split()[5])
		dz2 = float(line.split()[6])
		dzx = float(line.split()[7])
		dx2_y2 = float(line.split()[8])
		if vcounter%num_bands < first_n_bands:
			tmp_dft_eigenvector.append(np.array([s,px,py,pz,dxy,dyz,dzx,dx2_y2,dz2]))
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

dft_path = np.delete(tmp_dft_path,np.s_[::100],0)
dft_path = np.insert(dft_path,0,tmp_dft_path[0],axis=0)
dft_bands= np.delete(tmp_dft_bands,np.s_[::100],0)
dft_bands = np.insert(dft_bands,0,tmp_dft_bands[0],axis=0)
dft_eigenvectors= np.delete(tmp_dft_eigenvectors,np.s_[::100],0)
dft_eigenvectors = np.insert(dft_eigenvectors,0,tmp_dft_eigenvectors[0],axis=0)

del tmp_dft_path
del tmp_dft_bands
del tmp_dft_eigenvectors

dft_path = np.concatenate((dft_path[:201][:],dft_path[700:][:]),axis=0)
dft_bands = np.concatenate((dft_bands[:201][:],dft_bands[700:][:]),axis=0)
dft_eigenvectors = np.concatenate((dft_eigenvectors[:201][:][:],dft_eigenvectors[700:][:][:]),axis=0)
	
ticks = [0]
tick = 0
for i in range(1,len(dft_path)):
	d = np.linalg.norm(dft_path[i]-dft_path[i-1])
	if d < 1.0:
		tick += d
	else:
		tick += 0.000001
	ticks.append(tick)

fig, ax = plt.subplots()

xcoords = [0,0.8660254037844393,1.866025403784445,2.246606140472424,3.296553569714969]
for xc in xcoords:
    plt.axvline(x=xc,c='grey')

for i,band in enumerate(dft_bands.T):
	ax.plot(ticks,band,c = 'k')
	if i == 0:
		l0, = ax.plot(ticks,band,c='k')
# plt.show()

###########################################################################################################
#
# Get TB Bandstructure
#
###########################################################################################################

es = -4.1529
ep = 3.0471
ed = 6.99 
ed_ = 19.5033
vss = -8.23
vxx = 2.0471
vxy = 4.2553
vsp = 5.9134
vsdS = 0.8088
vpdS = -5.1862
vpdP = -11.395
vddS = -1.0481
vddP = 13.6586
vddD = -5.8832
params = np.array([es, ep, ed, ed_, vss, vxx, vxy, vsp, vsdS, vpdS, vpdP, vddS, vddP, vddD])

# tb_bands = []
# for k in dft_path:
# 	e, vs = band_eigenvalues_and_eigenvectors(params,k[0],k[1],k[2],dft_vbm)
# 	tb_bands.append(np.asarray(e)) 
# tb_bands = np.asarray(tb_bands)

# for i,band in enumerate(tb_bands.T):
# 	ax.plot(ticks,band,c='r')
# 	if i == 0:
# 		l1, = ax.plot(ticks,band,c='r')
# ax.legend((l0,l1),["SCAN-DFT","TB-sp3d5"], loc='upper right',shadow=True)
# plt.xticks(xcoords,(r'$L$',r'$\Gamma$',r'$X$',r'$U,K$',r'$\Gamma$'),fontsize=20)
# plt.yticks(fontsize=20)
# plt.ylabel('Energy (eV)',fontsize=30)
# plt.title('Si Band Structure',fontsize=30)

# plt.show()
###########################################################################################################
#
# Fitting
#
###########################################################################################################

path_length = 250
every = 99
res = minimize(residual, params, args=(dft_path[:path_length:every],dft_bands[:path_length:every],dft_eigenvectors[:path_length:every]),method='nelder-mead',
	options={'xtol': 1e-8, 'disp': True})

params = res.x

tb_bands = []
for k in dft_path:
	e, vs = band_eigenvalues_and_eigenvectors(params,k[0],k[1],k[2],dft_vbm)
	tb_bands.append(np.asarray(e)) 
tb_bands = np.asarray(tb_bands)

for i,band in enumerate(tb_bands.T):
	ax.plot(ticks,band,c='r')
	if i == 0:
		l1, = ax.plot(ticks,band,c='r')
ax.legend((l0,l1),["SCAN-DFT","TB-sp3d5"], loc='upper right',shadow=True)
plt.xticks(xcoords,(r'$L$',r'$\Gamma$',r'$X$',r'$U,K$',r'$\Gamma$'),fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('Energy (eV)',fontsize=30)
plt.title('Si Band Structure',fontsize=30)


plt.show()