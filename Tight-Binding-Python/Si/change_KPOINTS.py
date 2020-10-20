import numpy as np 

a = 5.430
c = 1

b1 = c*np.array([1,1,-1])
b2 = c*np.array([-1,1,1])
b3 = c*np.array([1,-1,1])

kpath = np.array([[0.5,0.5,0.5],[0.0,0.0,0.0],[0.5,0.0,0.5],[0.625,0.25,0.625],[0.375,0.375,0.75],[0.0,0.0,0.0]])

for k in kpath:
	new_k = k[0]*b1+k[1]*b2+k[2]*b3
	print(new_k[0],end="\t")
	print(new_k[1],end="\t")
	print(new_k[2],end="\n")