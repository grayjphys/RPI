import numpy as np
a = 5.43

## Read in k-points and build an array of points
f = open("kpoint_integration_points_280.dat")
kspace = []
box = []
for i,line in enumerate(f):
	if i%28 ==0 and i !=0:
		kspace.append(np.array(box))
		box = []
	elif i%28 != 0:
		box.append(np.array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])]))
kspace.append(box)
print(len(kspace))
kspace = np.array(kspace)

#Get the real part of the spherical harmonics times the Fourier phase
def sph_real(x,y,z,kx,ky,kz):
	real = np.cos(2*np.pi*(kx*x + ky*y + kz*z)/a)
	s = 0.5*np.sqrt(1/np.pi)
	py = 0.5*np.sqrt(3/np.pi)*(y/np.sqrt(x**2+y**2+z**2))
	pz = 0.5*np.sqrt(3/np.pi)*(z/np.sqrt(x**2+y**2+z**2))
	px = 0.5*np.sqrt(3/np.pi)*(x/np.sqrt(x**2+y**2+z**2))
	dxy = 0.5*np.sqrt(15/np.pi)*((x*y)/(x**2+y**2+z**2))
	dyz = 0.5*np.sqrt(15/np.pi)*((y*z)/(x**2+y**2+z**2))
	dz2 = 0.25*np.sqrt(5/np.pi)*((-x**2-y**2+2*z**2)/(x**2+y**2+z**2))
	dzx = 0.5*np.sqrt(15/np.pi)*((z*x)/(x**2+y**2+z**2))
	dx2_y2 = 0.25*np.sqrt(15/np.pi)*((x**2-y**2)/(x**2+y**2+z**2))
	return real*np.array([s,py,pz,px,dxy,dyz,dz2,dzx,dx2_y2])


#Get the imaginary part of the spherical harmonics times the Fourier phase
def sph_imag(x,y,z,kx,ky,kz):
	imag = np.sin(2*np.pi*(kx*x + ky*y + kz*z)/a)
	s = 0.5*np.sqrt(1/np.pi)
	py = 0.5*np.sqrt(3/np.pi)*(y/np.sqrt(x**2+y**2+z**2))
	pz = 0.5*np.sqrt(3/np.pi)*(z/np.sqrt(x**2+y**2+z**2))
	px = 0.5*np.sqrt(3/np.pi)*(x/np.sqrt(x**2+y**2+z**2))
	dxy = 0.5*np.sqrt(15/np.pi)*((x*y)/(x**2+y**2+z**2))
	dyz = 0.5*np.sqrt(15/np.pi)*((y*z)/(x**2+y**2+z**2))
	dz2 = 0.25*np.sqrt(5/np.pi)*((-x**2-y**2+2*z**2)/(x**2+y**2+z**2))
	dzx = 0.5*np.sqrt(15/np.pi)*((z*x)/(x**2+y**2+z**2))
	dx2_y2 = 0.25*np.sqrt(15/np.pi)*((x**2-y**2)/(x**2+y**2+z**2))
	return imag*np.array([s,py,pz,px,dxy,dyz,dz2,dzx,dx2_y2])

# Change the range of the function to -1 to 1. The typical range for Gaussian quadrature
def g_sph(f,x,y,z,kx,ky,kz,ax,bx,ay,by,az,bz):
	return f(0.5*(bx-ax)*x+0.5*(bx+ax),0.5*(by-ay)*y+0.5*(by+ay),0.5*(bz-az)*z+0.5*(bz+az),kx,ky,kz)*0.5*(bx-ax)*0.5*(by-ay)*0.5*(bz-az)

# Gaussian 3 point quadrature weights and points
weights = [5/9,8/9,5/9]
points = [-0.6**0.5,0,0.6**0.5]

# The position space integral is restricted to a sphere in cartesian coordinates
# Gaussian quadrature in spherical coordinates is not trivial
RWIGS = 2.480
num = 10
count = 0
boxes_real = []
boxes_imag = []
# Integrate every box and sum the results
for b,box in enumerate(kspace):
	print(b)
	integrals_real = []
	integrals_imag = []
	ax = kspace[b][0][0]
	ay = kspace[b][0][1]
	az = kspace[b][0][2]
	bx = kspace[b][-1][0]
	by = kspace[b][-1][1]
	bz = kspace[b][-1][2]
	for kpoint in box:
		integral_real = 0
		integral_imag = 0
		for ni in range(num+1):
			for nj in range(num+1):
				for nk in range(num+1):
					# Find the initial bounds of the box so that they can be fed into the g function to change the range to -1 to 1
					corner1 = ((-RWIGS+(2*RWIGS*ni/num))**2 + (-RWIGS+(2*RWIGS*nj/num))**2 + (-RWIGS+(2*RWIGS*nk/num))**2)**0.5
					corner2 = ((-RWIGS+(2*RWIGS*ni/num))**2 + (-RWIGS+(2*RWIGS*nj/num))**2 + (-RWIGS+(2*RWIGS*nk/num)+(2*RWIGS/num))**2)**0.5
					corner3 = ((-RWIGS+(2*RWIGS*ni/num))**2 + (-RWIGS+(2*RWIGS*nj/num)+(2*RWIGS/num))**2 + (-RWIGS+(2*RWIGS*nk/num))**2)**0.5
					corner4 = ((-RWIGS+(2*RWIGS*ni/num))**2 + (-RWIGS+(2*RWIGS*nj/num)+(2*RWIGS/num))**2 + (-RWIGS+(2*RWIGS*nk/num)+(2*RWIGS/num))**2)**0.5
					corner5 = ((-RWIGS+(2*RWIGS*ni/num)+(2*RWIGS/num))**2 + (-RWIGS+(2*RWIGS*nj/num))**2 + (-RWIGS+(2*RWIGS*nk/num))**2)**0.5
					corner6 = ((-RWIGS+(2*RWIGS*ni/num)+(2*RWIGS/num))**2 + (-RWIGS+(2*RWIGS*nj/num))**2 + (-RWIGS+(2*RWIGS*nk/num)+(2*RWIGS/num))**2)**0.5
					corner7 = ((-RWIGS+(2*RWIGS*ni/num)+(2*RWIGS/num))**2 + (-RWIGS+(2*RWIGS*nj/num)+(2*RWIGS/num))**2 + (-RWIGS+(2*RWIGS*nk/num))**2)**0.5
					corner8 = ((-RWIGS+(2*RWIGS*ni/num)+(2*RWIGS/num))**2 + (-RWIGS+(2*RWIGS*nj/num)+(2*RWIGS/num))**2 + (-RWIGS+(2*RWIGS*nk/num)+(2*RWIGS/num))**2)**0.5
					# If the box is inside the radius sum the integral
					if corner1 <= RWIGS and corner2 <= RWIGS and corner3 <= RWIGS and corner4 <= RWIGS and corner5 <= RWIGS and corner6 <= RWIGS and corner7 <= RWIGS and corner8 <= RWIGS: 
						for i in range(3):
							for j in range(3):
								for k in range(3):
									integral_real += g_sph(sph_real,points[i],points[j],points[k],0.5*(bx-ax)*kpoint[0]+0.5*(bx+ax),0.5*(by-ay)*kpoint[1]+0.5*(by+ay),0.5*(bz-az)*kpoint[2]+0.5*(bz+az),-RWIGS+(2*RWIGS*ni/num),-RWIGS+(2*RWIGS*ni/num)+(2*RWIGS/num),-RWIGS+(2*RWIGS*nj/num),-RWIGS+(2*RWIGS*nj/num)+(2*RWIGS/num),-RWIGS+(2*RWIGS*nk/num),-RWIGS+(2*RWIGS*nk/num)+(2*RWIGS/num))*weights[i]*weights[j]*weights[k]
									integral_imag += g_sph(sph_imag,points[i],points[j],points[k],0.5*(bx-ax)*kpoint[0]+0.5*(bx+ax),0.5*(by-ay)*kpoint[1]+0.5*(by+ay),0.5*(bz-az)*kpoint[2]+0.5*(bz+az),-RWIGS+(2*RWIGS*ni/num),-RWIGS+(2*RWIGS*ni/num)+(2*RWIGS/num),-RWIGS+(2*RWIGS*nj/num),-RWIGS+(2*RWIGS*nj/num)+(2*RWIGS/num),-RWIGS+(2*RWIGS*nk/num),-RWIGS+(2*RWIGS*nk/num)+(2*RWIGS/num))*weights[i]*weights[j]*weights[k]
						count += 1/(len(kspace)*27)

		integrals_real.append(integral_real)
		integrals_imag.append(integral_imag)
	boxes_real.append(np.array(integrals_real)*(27*count)**-1)
	boxes_imag.append(np.array(integrals_imag)*(27*count)**-1)
print(count)
boxes_real = np.array(boxes_real)
boxes_imag = np.array(boxes_imag)

# Write all of the data to a file
f = open("spherical_harmonics_at_kpoints_280.dat",'w')

for b in range(len(boxes_real)):
	for i in range(len(boxes_real[0])):
		for j in range(len(boxes_real[0][0])):
			line = str(boxes_real[b][i][j]) + "\t" + str(boxes_imag[b][i][j]) + "\n"
			f.write(line)
