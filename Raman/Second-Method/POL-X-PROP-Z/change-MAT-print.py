import os

f = open("GEO-OPT.inp",'r')
f2 = open("GEO-OPT",'w+')
for line in f:
	if "         KOHN_SHAM_MATRIX .TRUE." in line:
		f2.write(line)
                f2.write("         KINETIC_ENERGY .TRUE.\n")
	else:
		f2.write(line)
os.system("mv GEO-OPT GEO-OPT.inp")
