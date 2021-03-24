import os

path = os.getcwd()
dname = os.path.basename(path)
f = open("GEO-OPT.inp",'r')
f2 = open("GEO",'w+')
for line in f:
	if "&MO ON" in line or "EIGVECS" in line or "CARTESIAN" in line or "FILENAME cartesian-mos" in line or "&EACH" in line or "&END EACH" in line or "QS_SCF 0" in line or "&END MO" in line or "MOMENTS" in line or "        PERIODIC .FALSE." in line or "MAX_MOMENT 1" in line or "        ADD_LAST NUMERIC" in line or "COMMON_ITERATION_LEVELS 1" in line or "MD 1" in line:
		x = 1
        else:
		f2.write(line)
os.system("mv GEO GEO-OPT.inp")
