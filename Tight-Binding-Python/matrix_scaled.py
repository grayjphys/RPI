import numpy as np
import cmath 

kx = 0.5
ky = 0.5
kz = 0.5 

g = 0.25*np.array([cmath.rect(1,(np.pi/2)*(kx+ky+kz)) +cmath.rect(1,(np.pi/2)*(kx-ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
					   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) +cmath.rect(1,(np.pi/2)*(kx-ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
					   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) -cmath.rect(1,(np.pi/2)*(kx-ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx-ky+kz)),
					   cmath.rect(1,(np.pi/2)*(kx+ky+kz)) -cmath.rect(1,(np.pi/2)*(kx-ky-kz)) -cmath.rect(1,(np.pi/2)*(-kx+ky-kz)) +cmath.rect(1,(np.pi/2)*(-kx-ky+kz))],dtype=complex)

gc = np.conj(g)

es = -4.1529
ep = 3.0471
vss = -8.23
vxx = 2.0471
vxy = 4.2553
vsp = 5.9134

hamiltonian = np.array([
		[es,	     vss*g[0],  0.0,		0.0,		0.0,		vsp*g[1],  vsp*g[2],  vsp*g[3]],
		[vss*gc[0],  es,	   -vsp*gc[1], -vsp*gc[2], -vsp*gc[3],  0.0,	   0.0,		  0.0     ],
		[0.0,	    -vsp*g[1],  ep,		    0.0,		0.0,		vxx*g[0],  vxy*g[3],  vxy*g[2]],
		[0.0,       -vsp*g[2],  0.0,        ep,         0.0,        vxy*g[3],  vxx*g[0],  vxy*g[1]],
		[0.0,       -vsp*g[3],  0.0,        0.0,        ep,         vxy*g[2],  vxy*g[1],  vxx*g[0]],
		[gc[1]*vsp,  0.0,       vxx*gc[0],  vxy*gc[3],  vxy*gc[2],  ep,        0.0,       0.0     ],
		[gc[2]*vsp,  0.0,       vxy*gc[3],  vxx*gc[0],  vxy*gc[1],  0.0,       ep,        0.0     ],
		[gc[3]*vsp,  0.0,       vxy*gc[2],  vxy*gc[1],  vxx*gc[0],  0.0,       0.0,       ep     ]],dtype=complex)


atom_number = 4

h = np.zeros((atom_number*2,atom_number*2))
for i in range(atom_number*2):
	for j in range(atom_number*2):
		

print(hamiltonian)
