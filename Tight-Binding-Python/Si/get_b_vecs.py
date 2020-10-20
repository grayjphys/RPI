import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

f = open("POSCAR")

a1 = []
a2 = []
a3 = []
for i,line in enumerate(f):
	tmp = line.split()
	if i == 2:
		a1 = np.array([float(tmp[0]),float(tmp[1]),float(tmp[2])])
	if i == 3:
		a2 = np.array([float(tmp[0]),float(tmp[1]),float(tmp[2])])
	if i == 4:
		a3 = np.array([float(tmp[0]),float(tmp[1]),float(tmp[2])])

volume = a1.dot(np.cross(a2,a3))
b1 = np.cross(a2,a3)/volume
b2 = np.cross(a3,a1)/volume
b3 = np.cross(a1,a2)/volume

print(b1,b2,b3)

positions = []
for i in range(-1,2):
	for j in range(-1,2):
		for k in range(-1,2):
			positions.append(i*b1+j*b2+k*b3)

positions = np.array(positions).T
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')



ax.scatter(positions[0], positions[1], positions[2],s = 20)
ax.scatter([0],[0],[0],s = 30, c = 'r')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()