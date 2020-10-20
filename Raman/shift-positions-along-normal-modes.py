import numpy as np

molecule = "pyridine"
N = 11 # Number of atoms
M = 27 # Number of modes
dx_scale = 0.001 # The largest value of a vibrational mode displacement is about
                 # 0.25 angstroms, so 0.25*dx_scale is about 0.00025. With central differences, the
                 # error is (0.00025)^2 = 6.25e-9
bohr_to_ang = 0.529177
shift = np.array([10.0,10.0,10.0])
# Reads molden file and creates displaced coordinates along normal modes

f = open("PYRIDINE-VIBRATIONS-1.mol")
data = f.readlines()
positions = []
atom_list = []
pos = []
vibrational_modes = []
vibrational_mode = []
mod_atom_disp = []
collect = False
current_i = 0
for i,dat in enumerate(data):
	if 1<i <= N + 1:
		pos = dat.split()
		positions.append([float(pos[3])*bohr_to_ang,float(pos[4])*bohr_to_ang,float(pos[5])*bohr_to_ang])
		atom_list.append(pos[0])
	if N + 1 < i < (N+3)*M:
		if "vibration" in dat:
			vibrational_mode = []
			current_i = i
		if current_i < i <= current_i + N:
			mod_atom_disp =[float(dat.split()[0])*bohr_to_ang*dx_scale,float(dat.split()[1])*bohr_to_ang*dx_scale,float(dat.split()[2])*bohr_to_ang*dx_scale]
			vibrational_mode.append(mod_atom_disp)
		if i == current_i + N:
			vibrational_modes.append(vibrational_mode)

positions = np.array(positions)
# Required to align Nitrogen with Silver sphere in the case of Pyridene

center_of_ring = (1.0/6.0)*np.array([ sum(x) for x in zip(*positions[:6]) ])

#Rotation_Z = np.array([[0.86602540378, -0.5, 0.0],
#					   [0.5,  0.86602540378, 0.0],
#					   [0.0,			0.0, 1.0]])

#rotated_positions = np.array([Rotation_Z.dot(pos) for pos in positions - center_of_ring])



p = open("{}.xyz".format(molecule),'w')
p.write("{}\n".format(N))
p.write(molecule + "\n")

vibrational_modes_scaled = np.array(vibrational_modes)

positions_plus_modes = vibrational_modes_scaled + positions - center_of_ring + shift
positions_minus_modes = -1*(vibrational_modes_scaled - positions) - center_of_ring + shift

for a, atom in enumerate(positions-center_of_ring+shift):
        p.write("{} {} {} {}\n".format(atom_list[a],"{:.9f}".format(atom[0]),"{:.9f}".format(atom[1]),"{:.9f}".format(atom[2])))

#for m in range(M):
#	for n in range(N):
#		positions_plus_modes[m][n] = Rotation_Z.dot(positions_plus_modes[m][n] - center_of_ring)

#for m in range(M):
#	for n in range(N):
#		positions_minus_modes[m][n] = Rotation_Z.dot(positions_minus_modes[m][n] - center_of_ring)

modesscaled = open("{}.scaled-modes".format(molecule),'w')
modesscaled.write("Vibrational Modes Scaled (0.001)\n")

for i in range(M):
	modesscaled.write("Mode {}\n".format(i+1))
	vibrational_mode_scaled = vibrational_modes_scaled[i]

	for a, atom in enumerate(vibrational_mode_scaled):
		modesscaled.write("{} {} {} {}\n".format(atom_list[a],"{:.9f}".format(atom[0]),"{:.9f}".format(atom[1]),"{:.9f}".format(atom[2])))
	modesscaled.write("\n")

	modeplus = open("{}-plus-mode-{}.xyz".format(molecule,i+1),'w')
	modeplus.write("{}\n".format(N))
	modeplus.write("Pyridine Relaxed Plus Scaled (0.001) Mode {}\n".format(i+1))

	positions_plus_mode = positions_plus_modes[i]
	for a, atom in enumerate(positions_plus_mode):
		modeplus.write("{} {} {} {}\n".format(atom_list[a],"{:.9f}".format(atom[0]),"{:.9f}".format(atom[1]),"{:.9f}".format(atom[2])))

	modeminus = open("{}-minus-mode-{}.xyz".format(molecule,i+1),'w')
	modeminus.write("{}\n".format(N))
	modeminus.write("Pyridine Relaxed Minus Scaled (0.001) Mode {}\n".format(i+1))
	
	positions_minus_mode = positions_minus_modes[i]
	for a, atom in enumerate(positions_minus_mode):
		modeminus.write("{} {} {} {}\n".format(atom_list[a],"{:.9f}".format(atom[0]),"{:.9f}".format(atom[1]),"{:.9f}".format(atom[2])))









#IGNORE THIS


# out = subprocess.Popen(["grep", "-A", "11", "[Atoms]","PYRIDINE-VIBRATIONS-1.mol"], 
#            stdout=subprocess.PIPE, 
#            stderr=subprocess.STDOUT)


# stdout,stderr = out.communicate()

# data = [i.decode("utf-8") for i in stdout.split()]

# for i,dat in enumerate(data):
# 	if 3 < i < N*6 + 4 and not (i-3)%6 == 2  and not (i-3)%6 == 3:
# 		if (i-3)%6 == 1:
# 			pos = []
# 			pos.append(dat)
# 		if (i-3)%6 == 4:
# 			pos.append(float(dat))
# 		if (i-3)%6 == 5:
# 			pos.append(float(dat))
# 		if (i-3)%6 == 0:
# 			pos.append(float(dat))
# 			positions.append(pos)
# 	if N*6 + 4 < i and (i + 2) % (3*N + 3) > 2:
# 		if (i + 2) % (3*N + 3) == 3:
# 			vibrational_mode = []
# 		if ((i + 2) % (3*N + 3))%3 == 1:
# 			mod_atom_disp = []
# 			mod_atom_disp.append(float(dat))
# 		if ((i + 2) % (3*N + 3))%3 == 2:
# 			mod_atom_disp.append(float(dat))
# 		if ((i + 2) % (3*N + 3))%3 == 0:
# 			mod_atom_disp.append(float(dat))
# 			vibrational_mode.append(mod_atom_disp)
# 		if (i + 2) % (3*N + 3) == 3:
# 			vibrational_modes.append(vibrational_mode)

# for m,v in enumerate(vibrational_modes):
# 	print("vibration")
# 	print(m+1)
# 	for atom in v:
# 		print(atom)

