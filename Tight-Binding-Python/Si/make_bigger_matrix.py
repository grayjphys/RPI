import numpy as np
import cmath

kx = 0.5
ky = 0.5
kz = 0.5

es = -4.1529
ep = 3.0471
ed = 6.99 
ed_ = 19.5033
vss = -8.23
vxx = 2.0471
vxy = 4.2553
vsp = 5.9134


g = 0.25*np.array([cmath.rect(1,(np.pi/2)*(kx+ky+kz)) +cmath.rect(1,(np.pi/2)*(kx-ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
				   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) +cmath.rect(1,(np.pi/2)*(kx-ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
				   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) -cmath.rect(1,(np.pi/2)*(kx-ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
				   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) -cmath.rect(1,(np.pi/2)*(kx-ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx-ky+kz))],dtype=complex)

gc = np.conj(g)

hamiltonian = np.array([
	[es, 	0.0,	0.0,	0.0,         vss*g[0], vsp*g[1], vsp*g[2], vsp*g[3]],
	[0.0, 	ep,  	0.0,	0.0,        -vsp*g[1], vxx*g[0], vxy*g[3], vxy*g[2]],
	[0.0,	0.0,	ep,	    0.0,        -vsp*g[2], vxy*g[3], vxx*g[0], vxy*g[1]],
	[0.0,	0.0,	0.0,	ep,         -vsp*g[3], vxy*g[2], vxy*g[1], vxx*g[0]],
	[vss*gc[0],-vsp*gc[1],-vsp*gc[2],-vsp*gc[3],	es, 	0.0, 	0.0,  	0.0],
	[gc[1]*vsp, vxx*gc[0], vxy*gc[3], vxy*gc[2],	0.0,	ep, 	0.0,	0.0],
	[gc[2]*vsp, vxy*gc[3], vxx*gc[0], vxy*gc[1],	0.0,	0.0,	ep,	    0.0],
	[gc[3]*vsp, vxy*gc[2], vxy*gc[1], vxx*gc[0],	0.0,	0.0,	0.0,	ep]],dtype=complex)

self_ham = np.array([[es, 	0.0, 	0.0,  	0.0],
					 [0.0,	ep, 	0.0,	0.0],
					 [0.0,	0.0,	ep,	    0.0],
					 [0.0,	0.0,	0.0,	ep]],dtype=complex)

interact_ham= np.array([[vss*g[0], vsp*g[1], vsp*g[2], vsp*g[3]],
					   [-vsp*g[1], vxx*g[0], vxy*g[3], vxy*g[2]],
					   [-vsp*g[2], vxy*g[3], vxx*g[0], vxy*g[1]],
					   [-vsp*g[3], vxy*g[2], vxy*g[1], vxx*g[0]]],dtype=complex)
conj_interact_ham = np.conj(interact_ham.T)

num_atoms = 3
new_hamiltonian = np.zeros((num_atoms,num_atoms,4,4),dtype=complex)
for i in range(num_atoms):
	for j in range(num_atoms):
		if i == j:
			new_hamiltonian[i][j] = self_ham
		if i == j + 1:
			new_hamiltonian[i][j] = interact_ham
		if j == i + 1:
			new_hamiltonian[i][j] = conj_interact_ham

new_hamiltonian = np.concatenate(new_hamiltonian,axis=2)
new_hamiltonian = np.concatenate(new_hamiltonian,axis=0)
print(new_hamiltonian)