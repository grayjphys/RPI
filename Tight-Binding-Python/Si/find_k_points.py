import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt 

# The number of intervals is equal to num-1 in one dimension
# The total number of boxes is (num-1)**3
num = 11

# Define the region of integration in units of 2pi/a
kxs = np.linspace(-1,1,num)
kys = np.linspace(-1,1,num)
kzs = np.linspace(-1,1,num)

#REMEMBER THE RECIPROCAL LATTICE IS A BCC LATTICE
# Points of the BCC lattice. First nearest and second nearest neighbors to the origin
xpos = np.array([[0,0,0],
				 [1,1,1],
				 [1,1,-1],
				 [1,-1,1],
				 [1,-1,-1],
				 [-1,1,1],
				 [-1,1,-1],
				 [-1,-1,1],
				 [-1,-1,-1],
				 [2,0,0],
				 [-2,0,0],
				 [0,2,0],
				 [0,-2,0],
				 [0,0,2],
				 [0,0,-2]])

# Gaussian 3 point quadrature points
points = [-0.6**0.5,0,0.6**0.5]
ks = []
for i in range(num-1):
	for j in range(num-1):
		for k in range(num-1):
			box = []
			for px in points:
				for py in points:
					for pz in points:
						# Define x, y, z coordinate within a box
						x = kxs[i] + 1/(num-1) + px/(num-1)
						y = kys[j] + 1/(num-1) + py/(num-1)
						z = kzs[k] + 1/(num-1) + pz/(num-1)
						# start with a large maximum
						d = 10
						nearest = []
						# Find the nearest atom to the point by comparing to the maximum distance
						for atom in xpos[1:9]:
							if ((atom[0]-x)**2 + (atom[1]-y)**2 + (atom[2]-z)**2)**0.5 < d:
								d = ((atom[0]-x)**2 + (atom[1]-y)**2 + (atom[2]-z)**2)**0.5
								nearest = atom
						p1 = np.array([x,y,z])
						# The unit vector of the nearest neighbor lattice position
						p2 = nearest/(np.linalg.norm(nearest))
						# The wigner seitz cell consists of the convex hull of bisections between 
						# nearest neighbors
						cutoff = np.linalg.norm(nearest)/2
						if abs(p1.dot(p2)) <= cutoff:
							box.append(p1)
			# If all of the points in a box are in the Wigner-Seitz Cell 
			# then the box is in the Wigner-Seitz cell
			if len(box) == 27:
				ks.append(np.array(box))

###################################################################
## Output the points to a file and plot the First Brillouin Zone ##
###################################################################
ks = np.array(ks)
print(len(ks))
f = open("kpoint_integration_points_280.dat",'w')
for b,box in enumerate(ks):
	box_num = str("box: %d"%b) + "\n"
	f.write(box_num)
	for k in box:
		k_string = str(round(k[0],2)) + "\t" + str(round(k[1],2)) + "\t" + str(round(k[2],2)) + "\n"
		f.write(k_string)
ks = ks.transpose(0,2,1)

xpos_in = xpos[:9].T
xpos_out = xpos[9:].T

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for box in ks:
	ax.scatter(box[0][27//2], box[1][27//2], box[2][27//2],s = 3,c='b')
ax.scatter(xpos_in[0], xpos_in[1], xpos_in[2],s = 20)
ax.scatter(xpos_out[0], xpos_out[1], xpos_out[2],s = 20)
ax.set_xlabel(r'$k_x$')
ax.set_ylabel(r'$k_y$')
ax.set_zlabel(r'$k_z$')

plt.show()
