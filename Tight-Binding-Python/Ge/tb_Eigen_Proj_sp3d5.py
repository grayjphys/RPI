import numpy as np
import cmath
import matplotlib.pyplot as plt 
from scipy.optimize import minimize
import scipy.linalg as LA

# np.__config__.show()
# scipy.__config__.show()

a = 5.468721419 #Angstroms exp 5.431
dft_vbm =5.67014

f = open("kpoint_integration_points.dat")
kspace = []
box = []
for i,line in enumerate(f):
	if i%28 ==0 and i !=0:
		kspace.append(np.array(box))
		box = []
	elif i%28 != 0:
		box.append(np.array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])]))
kspace.append(box[:])
kspace = np.array(kspace)

boxes = {}
s = []
py = []
pz = []
px = []
dxy = []
dyz = []
dz2 = []
dzx = []
dx2_y2 = []
iterate = 0
f = open("spherical_harmonics_at_kpoints.dat")
for i,line in enumerate(f):
	if i%(27*9) == 0:
		if i != 0:
			boxes.update({i//(27*9)-1 : np.asarray([s,py,pz,px,dxy,dyz,dz2,dzx,dx2_y2])})
			s = []
			py = []
			pz = []
			px = []
			dxy = []
			dyz = []
			dz2 = []
			dzx = []
			dx2_y2 = []
	if i%9 == 0:
		s.append(complex(float(line.split()[0]),float(line.split()[1])))
	if i%9 == 1:
		py.append(complex(float(line.split()[0]),float(line.split()[1])))
	if i%9 == 2:
		px.append(complex(float(line.split()[0]),float(line.split()[1])))
	if i%9 == 3:
		pz.append(complex(float(line.split()[0]),float(line.split()[1])))
	if i%9 == 4:
		dxy.append(complex(float(line.split()[0]),float(line.split()[1])))
	if i%9 == 5:
		dyz.append(complex(float(line.split()[0]),float(line.split()[1])))
	if i%9 == 6:
		dz2.append(complex(float(line.split()[0]),float(line.split()[1])))
	if i%9 == 7:
		dzx.append(complex(float(line.split()[0]),float(line.split()[1])))
	if i%9 == 8:
		dx2_y2.append(complex(float(line.split()[0]),float(line.split()[1])))
boxes.update({i//(27*9)-1 : np.asarray([s,py,pz,px,dxy,dyz,dz2,dzx,dx2_y2])})

weights = []
for i in [5/9,8/9,5/9]:
	for j in [5/9,8/9,5/9]:
		for k in [5/9,8/9,5/9]:
			weights.append(i*j*k)

#Have all boxes and kpoints in boxes with spherical harmonics evaluated at those points
#Now perform the gaussian integration for these k points and 

def project_wave_functions(params):
	# [s,px,py,pz,dxy,dyz,dzx,dx2_y2,dz2]
	proj_vec = 0
	for b,box in enumerate(boxes):
		ax = kspace[b][0][0]
		ay = kspace[b][0][1]
		az = kspace[b][0][2]
		bx = kspace[b][-1][0]
		by = kspace[b][-1][1]
		bz = kspace[b][-1][2]
		for w,k in enumerate(kspace[b]):
			es , vs = band_eigenvalues_and_eigenvectors(params,0.5*(bx-ax)*k[0]+0.5*(bx+ax),0.5*(by-ay)*k[1]+0.5*(by+ay),0.5*(bz-az)*k[2]+0.5*(bz+az),dft_vbm,[],1)
			vs*=0.5*(bx-ax)*0.5*(by-ay)*0.5*(bz-az)*weights[w]
			vs[0]*=boxes[b][0][w]
			vs[1]*=boxes[b][3][w]
			vs[2]*=boxes[b][1][w]
			vs[3]*=boxes[b][2][w]
			vs[4]*=boxes[b][4][w]
			vs[5]*=boxes[b][5][w]
			vs[6]*=boxes[b][7][w]
			vs[7]*=boxes[b][8][w]
			vs[8]*=boxes[b][6][w]
			vs[9]*=boxes[b][0][w]
			vs[10]*=boxes[b][3][w]
			vs[11]*=boxes[b][1][w]
			vs[12]*=boxes[b][2][w]
			vs[13]*=boxes[b][4][w]
			vs[14]*=boxes[b][5][w]
			vs[15]*=boxes[b][7][w]
			vs[16]*=boxes[b][8][w]
			vs[17]*=boxes[b][6][w]
			if w == 0 and b == 0:
				proj_vec = vs
			else:
				proj_vec += vs

	return proj_vec

def band_eigenvalues_and_eigenvectors(params,kx,ky,kz,dft_vbm,proj_vec,type_):
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
	v18 = 0
	v19 = 0
	v25 = ((3**0.5)*vpdS + vpdP)/(3.0*(3.0**0.5))
	v26 = ((3**0.5)*vpdS - (2.0*vpdP))/(3.0*(3.0**0.5))
	v27 = v25
	v28 = vpdP/(3.0**0.5)
	v29 = -1.0*vpdP/3.0
	v38 = -1.0*v28
	v39 = v29
	v48 = 0.0
	v49 = -2.0*v29
	v55 = ((3.0*vddS)+(2.0*vddP)+(4.0*vddD))/9.0
	v56 = ((3.0*vddS)-vddP-(2.0*vddD))/9.0
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

	g = 0.25*np.array([cmath.rect(1,(np.pi/2)*(kx+ky+kz)) +cmath.rect(1,(np.pi/2)*(kx-ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
					   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) +cmath.rect(1,(np.pi/2)*(kx-ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
					   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) -cmath.rect(1,(np.pi/2)*(kx-ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
					   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) -cmath.rect(1,(np.pi/2)*(kx-ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx-ky+kz))],dtype=complex)

	gc = np.conj(g)

	hamiltonian = np.array([
		[es, 	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,							 	vss*g[0], vsp*g[1], vsp*g[2], vsp*g[3], v15*g[3], v15*g[1], v15*g[2], v18*g[0], v19*g[0]],
		[0.0, 	ep,  	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,							   -vsp*g[1], vxx*g[0], vxy*g[3], vxy*g[2], v25*g[2], v26*g[0], v27*g[3], v28*g[1], v29*g[1]],
		[0.0,	0.0,	ep,	    0.0,	0.0,	0.0,	0.0,	0.0,	0.0,							   -vsp*g[2], vxy*g[3], vxx*g[0], vxy*g[1], v27*g[1], v25*g[3], v26*g[0],-v28*g[2], v29*g[2]],
		[0.0,	0.0,	0.0,	ep,	    0.0,	0.0,	0.0,	0.0,	0.0,							   -vsp*g[3], vxy*g[2], vxy*g[1], vxx*g[0], v26*g[0], v27*g[2], v25*g[1], v48*g[3], v49*g[3]],
		[0.0,	0.0,	0.0,	0.0,	ed,	    0.0,	0.0,	0.0,	0.0,	 							v15*g[3],-v25*g[2],-v27*g[1],-v26*g[0], v55*g[0], v56*g[2], v56*g[1], v58*g[3], v59*g[3]],
		[0.0,	0.0,	0.0,	0.0,	0.0,	ed,	    0.0,	0.0,	0.0,	 							v15*g[1],-v26*g[0],-v25*g[3],-v27*g[2], v56*g[2], v55*g[0], v56*g[3], v68*g[1], v69*g[1]],
		[0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	ed,   	0.0,	0.0,	 							v15*g[2],-v27*g[3],-v26*g[0],-v25*g[1], v56*g[1], v56*g[3], v55*g[0], v78*g[2], v79*g[2]],
		[0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	ed_,	0.0,								v18*g[0],-v28*g[1], v28*g[2],-v48*g[3], v58*g[3], v68*g[1], v78*g[2], v88*g[0], v89*g[0]],
		[0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	ed_,	 							v19*g[0],-v29*g[1],-v29*g[2],-v49*g[3], v59*g[3], v69*g[1], v79*g[2], v89*g[0], v99*g[0]],
		[vss*gc[0],-vsp*gc[1],-vsp*gc[2],-vsp*gc[3], v15*gc[3], v15*gc[1], v15*gc[2], v18*gc[0], v19*gc[0],	es, 	0.0, 	0.0,  	0.0,	0.0,	0.0,	0.0,	0.0,	0.0],
		[gc[1]*vsp, vxx*gc[0], vxy*gc[3], vxy*gc[2],-v25*gc[2],-v26*gc[0],-v27*gc[3],-v28*gc[1],-v29*gc[1],	0.0,	ep, 	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0],
		[gc[2]*vsp, vxy*gc[3], vxx*gc[0], vxy*gc[1],-v27*gc[1],-v25*gc[3],-v26*gc[0], v28*gc[2],-v29*gc[2],	0.0,	0.0,	ep,	    0.0,	0.0,	0.0,	0.0,	0.0,	0.0],
		[gc[3]*vsp, vxy*gc[2], vxy*gc[1], vxx*gc[0],-v26*gc[0],-v27*gc[2],-v25*gc[1],-v48*gc[3],-v49*gc[3],	0.0,	0.0,	0.0,	ep,	    0.0,	0.0,	0.0,	0.0,	0.0],
		[v15*gc[3], v25*gc[2], v27*gc[1], v26*gc[0], v55*gc[0], v56*gc[2], v56*gc[1], v58*gc[3], v59*gc[3],	0.0,	0.0,	0.0,	0.0,	ed,	    0.0,	0.0,	0.0,	0.0],
		[v15*gc[1], v26*gc[0], v25*gc[3], v27*gc[2], v56*gc[2],	v55*gc[0], v56*gc[3], v68*gc[1], v69*gc[1],	0.0,	0.0,	0.0,	0.0,	0.0,	ed,   	0.0,	0.0,	0.0],
		[v15*gc[2], v27*gc[3], v26*gc[0], v25*gc[1], v56*gc[1],	v56*gc[3], v55*gc[0], v78*gc[2], v79*gc[2],	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	ed, 	0.0,	0.0],
		[v18*gc[0], v28*gc[1],-v28*gc[2], v48*gc[3], v58*gc[3],	v68*gc[1], v78*gc[2], v88*gc[0], v89*gc[0],	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	ed_,	0.0],
		[v19*gc[0], v29*gc[1], v29*gc[2], v49*gc[3], v59*gc[3],	v69*gc[1], v79*gc[2], v89*gc[0], v99*gc[0],	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	0.0,	ed_]],dtype=complex)

	if type_ == 0:
		vss = proj_vec
		es = LA.eigvalsh(hamiltonian)
	else:
		es, vss = LA.eigh(hamiltonian)
	VBM = len(es)*[dft_vbm]	
	return es+VBM, vss

def residual(params, dft_path, dft_bands, dft_eigenvectors):
	res = 0
	proj_vec = project_wave_functions(params)
	for p,k in enumerate(dft_path):
		tb_e, tb_v = band_eigenvalues_and_eigenvectors(params,k[0],k[1],k[2],dft_vbm,proj_vec,0)
		dft_e, dft_v = dft_bands[p],dft_eigenvectors[p]
		for i in range(0,10):
			# weight = np.exp(-abs(i-2)/22)
			weight = 1
			for j in range(0,10):
				res += weight*((tb_e[i] - dft_e[j])**2)*(np.linalg.norm(tb_v[i].dot(dft_v[j])))
	return res


###########################################################################################################
#
# Get Data from DFT
#
###########################################################################################################

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
eigvec = []
sr = 0
pyr = 0
pzr = 0
pxr = 0
dxyr = 0
dyzr = 0
dz2r = 0
dzxr = 0
dx2_y2r = 0
tmp_dft_eigenvector = []
tmp_dft_eigenvectors = []
f = open("SCAN_DAT")
for counter, line in enumerate(f):
	if counter == 0:
		num_kpoints =  int(line.split()[0])
	elif counter == 1:
		num_bands =  int(line.split()[0])
	elif counter == 2:
		num_ions =  int(line.split()[0])
	elif abs(counter-3) < num_kpoints:
		kpoint = float(line.split()[0])*np.array([1,1,-1])+float(line.split()[1])*np.array([-1,1,1])+float(line.split()[2])*np.array([1,-1,1])
		tmp_dft_path.append(kpoint)
	elif abs(counter-3) >= num_kpoints and abs(counter-3) < num_kpoints*(num_bands+1):
		if ecounter%num_bands == 0 and ecounter != 0:
			tmp_dft_bands.append(np.asarray(tmp_dft_band))
			tmp_dft_band = []
		en = float(line.split()[0])
		if ecounter%num_bands < first_n_bands:
			tmp_dft_band.append(en)
		ecounter+=1
	elif abs(counter-3) >= num_kpoints*(num_bands+1):
		if vcounter%(num_bands*num_ions*2)%(num_ions*2) == 0:
			eigvec = np.array([])
		if vcounter%(num_bands*num_ions*2) == 0 and vcounter != 0:
			tmp_dft_eigenvectors.append(np.asarray(tmp_dft_eigenvector))
			tmp_dft_eigenvector = []
		if vcounter%(num_bands*num_ions*2)%(num_ions*2)%2 == 0:
			sr = float(line.split()[1])
			pyr = float(line.split()[2])
			pzr = float(line.split()[3])
			pxr = float(line.split()[4])
			dxyr = float(line.split()[5])
			dyzr = float(line.split()[6])
			dz2r = float(line.split()[7])
			dzxr = float(line.split()[8])
			dx2_y2r = float(line.split()[9])
		if vcounter%(num_bands*num_ions*2)%(num_ions*2)%2 == 1:
			si = float(line.split()[1])
			pyi = float(line.split()[2])
			pzi = float(line.split()[3])
			pxi = float(line.split()[4])
			dxyi = float(line.split()[5])
			dyzi = float(line.split()[6])
			dz2i = float(line.split()[7])
			dzxi = float(line.split()[8])
			dx2_y2i = float(line.split()[9])
			s = complex(sr,si)
			py = complex(pyr,pyi)
			pz = complex(pzr,pzi)
			px = complex(pxr,pxi)
			dxy = complex(dxyr,dxyi)
			dyz = complex(dyzr,dyzi)
			dz2 = complex(dz2r,dz2i)
			dzx = complex(dzxr, dzxi)
			dx2_y2 = complex(dx2_y2r, dx2_y2i)
			eigvec = np.append(eigvec,[s,px,py,pz,dxy,dyz,dzx,dx2_y2,dz2])
		if vcounter%(num_bands*num_ions*2)%num_bands  < first_n_bands and vcounter%(num_bands*num_ions*2)%(num_ions*2) == num_ions*2-1 and vcounter!=0:
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

xcoords = [0,0.8660254037844393,1.866025403784441, 2.2195787943794696,3.269526223622345]
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
# es = -2.55
# ep = 4.48593
# ed = 14.01053
# ed_ = 23.44607
# vss = -1.86600
# vxx = (4.08481 - 2*1.49207)*(4/3)
# vxy = (4.08481 + 1.49207)*(4/3)
# vsp = 2.91067
# vsdS = -2.23992
# vpdS = -1.66657
# vpdP = 2.39936
# vddS = -1.82945
# vddP = 3.08177
# vddD = -1.56676
params = np.array([es, ep, ed, ed_, vss, vxx, vxy, vsp, vsdS, vpdS, vpdP, vddS, vddP, vddD])

###########################################################################################################
#
# Fitting
#
###########################################################################################################

path_start = 0
path_end = -1
every = 1

print(dft_path[path_start:path_end:every])
res = residual(params,dft_path[path_start:path_end:every],dft_bands[path_start:path_end:every],dft_eigenvectors[path_start:path_end:every])
print(res)
print(params)
res = minimize(residual, params, args=(dft_path[path_start:path_end:every],dft_bands[path_start:path_end:every],dft_eigenvectors[path_start:path_end:every]),method='nelder-mead',
	options={'ftol': 1e-4, 'disp': True,'adaptive':True,'maxiter':10000})

params = res.x
print(res)


tb_bands = []
for k in dft_path:
	e, vs = band_eigenvalues_and_eigenvectors(params,k[0],k[1],k[2],dft_vbm,[],1)
	tb_bands.append(np.asarray(e)) 
tb_bands = np.asarray(tb_bands)

for i,band in enumerate(tb_bands.T):
	ax.plot(ticks,band,c='r')
	if i == 0:
		l1, = ax.plot(ticks,band,c='r')
ax.legend((l0,l1),["SCAN-DFT","TB-sp3d5"], loc='upper right',shadow=True,fontsize=20)
plt.xticks(xcoords,(r'$L$',r'$\Gamma$',r'$X$',r'$U,K$',r'$\Gamma$'),fontsize=20)
plt.yticks(fontsize=20)
plt.ylim(-8,15)
plt.ylabel('Energy (eV)',fontsize=30)
plt.title('Si Band Structure',fontsize=30)


plt.show()