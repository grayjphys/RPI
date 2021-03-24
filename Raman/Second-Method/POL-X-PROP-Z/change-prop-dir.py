import os

f = open("GEO-OPT.inp",'r')
f2 = open("GEO-OPT",'w+')
for line in f:
	if "MAT_EXP TAYLOR" in line:
		f2.write(line)
		f2.write("      DELTA_PULSE_DIRECTION 0 0 1\n")
	else:
		f2.write(line)
os.system("mv GEO-OPT GEO-OPT.inp")
