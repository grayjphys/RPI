f = lambda x,y,z: (x*y*z) **2

points = [-0.90617984593866399279762687829939296512565191076253,
-0.53846931010568309103631442070020880496728660690556,
0,
0.53846931010568309103631442070020880496728660690556,
0.90617984593866399279762687829939296512565191076253]

weights = [0.23692688505618908751426404071991736264326000221241,
0.47862867049936646804129151483563819291229555334314,
0.5688888888888888888888888888888888888888888888889,
0.47862867049936646804129151483563819291229555334314,
0.23692688505618908751426404071991736264326000221241]

def g(n,x,y,z,ax,bx,ay,by,az,bz):
	diffx = (bx-ax)
	sumx = (bx+ax)/2
	diffy = (by-ay)
	sumy = (by+ay)/2
	diffz = (bz-az)
	sumz = (bz+az)/2
	return f(diffx*x+sumx,diffy*y+sumy,diffz*z+sumz)*diffx*diffy*diffz


integral = 0
n = 50
for x in range(n):
	for y in range(n):
		for z in range(n):
			for i in range(5):
				for j in range(5):
					for k in range(5):
						integral += g(n,points[i],points[j],points[k],(1/n)*x-1,(1/n)*x-1+(1/n),(1/n)*y-1,(1/n)*y-1+(1/n),(1/n)*z-1,(1/n)*z-1+(1/n))*weights[i]*weights[j]*weights[k]

print(integral)
