import numpy as np 
import random

random.seed(10)
N = 4
# A = np.zeros((N,N),dtype=complex)

# for i in range(N):
# 	for j in range(N):
# 		A[i][j] = complex(2 * random.random() - 1, 2 * random.random() - 1)

# A = (A + np.conj(A.T))/2



# k = 0
# b = np.zeros(N,dtype=complex)
# c = np.zeros(N-1,dtype=complex)
# gamma = 0
# phi = 0
# alpha = 0
# rho = complex(0,0)
# beta = complex(0,0)
# u = np.zeros(N,dtype=complex)
# t = np.zeros(N,dtype=complex)
# r = np.zeros(N,dtype=complex)

# for k in range(N-2):
# 	b[k] = A[k][k]
# 	gamma = 0
# 	for i in range(k+1,N):
# 		if np.linalg.norm(A[i][k]) > gamma:
# 			gamma = np.linalg.norm(A[i][k])
# 	phi = 0
# 	if gamma < 0.00000000001:
# 		c[k] = complex(0,0)
# 	else:
# 		alpha = 0
# 		if np.linalg.norm(A[k+1][k] < 0.0000000001):
# 			beta = complex(1,0)
# 		else:
# 			beta = A[k+1][k]/np.linalg.norm(A[k+1][k])
# 		for i in range(k+1,N):
# 			u[i] = A[i][k]/gamma
# 			alpha += np.linalg.norm(u[i])**2
# 		alpha = alpha**0.5
# 		c[k] = complex(-alpha,c[k].imag)
# 		phi = 1/(alpha*(alpha + np.linalg.norm(u[k+1])))
# 		u[k+1] += alpha*beta		
# 		rho = complex(0,0)

# 		for s in range(k+1,N):
# 			t[s] = complex(0,0)
# 			rho += np.conj(u[s])*t[s]

# 		for s in range(k+1,N):
# 			r[s] = rho*phi*u[s]/2

# 		for i in range(k+1,N):
# 			for j in range(k+1,i):
# 				A[i][j] += -t[i]*np.conj(u[j]) -np.conj(t[j])*u[j] + u[i]*np.conj(r[j]) + r[i]*np.conj(u[j])

# c[N-2] = A[N-1][N-2]
# b[N-2] = A[N-2][N-2]
# b[N-1] = A[N-1][N-1]

# Ak = np.zeros((N,N))

# for i in range(N):
# 	for j in range(N):
# 		if i == j:
# 			Ak[i][j] = b[i].real
# 		elif j == i + 1:
# 			if abs(c[i].imag) < 0.00000000001:
# 				if abs(c[i].real) < 0.00000000001:
# 					Ak[i][j] = 0
# 				else:
# 					Ak[i][j] = c[i].real
# 			else:
# 				Ak[i][j] = c[i].real
# 		elif i == j + 1:
# 			if abs(c[j].imag) < 0.00000000001:
# 				if abs(c[j].real) < 0.00000000001:
# 					Ak[i][j] = 0
# 				else:
# 					Ak[i][j] = c[j].real
# 			else:
# 				Ak[i][j] = c[j].real

Ak = np.asarray([[1.36075,-1.41243,0,0],
	 [-1.41243,-0.09041,-1.36329,0],
	 [0,-1.36329,-1.9348,0.01081],
	 [0,0,0.01081,0.05173]])

G = np.zeros((N-1,N,N))
row = 0
column = 0
temp = np.zeros((N,N))
temp2 = np.zeros((N,N))
sum_diags = 0
counter = 0
sign = 0
delta = 0
mu = 0
while(sum_diags < N):

	delta = (Ak[N-2][N-2]-Ak[N-1][N-1])/2
	if delta >= 0:
		sign = 1
	else:
		sign = -1
	mu = Ak[N-1][N-1] - sign*Ak[N-2][N-1]**2/(abs(delta) + (delta**2+Ak[N-2][N-1]**2)**0.5)

	for i in range(N):
		Ak[i][i] -= mu

	counter += 1
	for k in range(1,N):
		row = k
		column = k-1
		for i in range(N):
			for j in range(N):
				G[k-1][i][j] = 0
				if i == j and i != row - 1 and j != column and i != row  and j != column + 1:
					G[k-1][i][j] = 1
				elif (i == row - 1 and j == column) or (i == row and j == column + 1):
					G[k-1][i][j] = Ak[row-1][column]/(Ak[row-1][column]**2 + Ak[row][column]**2)**0.5
				elif (i == row and j == column) or (i == row - 1 and j == column + 1):
					G[k-1][i][j] = (j - i)*Ak[row][column]/(Ak[row-1][column]**2 + Ak[row][column]**2)**0.5

		Ak = G[k-1].dot(Ak)

	for k in range(1,N):
		Ak = Ak.dot(G[k-1].T)

	Ak += np.identity(N)*mu
	print("Ak:",k-1,counter)
	print(Ak)
	# print()
	sum_diags = 1
	for i in range(N):
		for j in range(N):
			if i == j + 1:
				if abs(Ak[i][j]) < 0.0000000001:
					sum_diags += 1
print([Ak[i][i] for i in range(N)])
print(counter)
