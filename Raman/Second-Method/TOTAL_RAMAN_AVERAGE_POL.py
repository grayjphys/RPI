import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

p = Path('.')

dirs = [x for x in p.iterdir() if x.is_dir() and "POL" in x.name]
dirs = np.sort(dirs)

num_wavenums = 1000
data = np.zeros((num_wavenums,2))
files_found = 0

for progress, d in enumerate(dirs):
    try:
        with open(d.name + "/" + d.name + "_SPECTRUM.dat",'r') as f:
            data += np.loadtxt(f,dtype=float)
        files_found += 1
    except:
        pass

plt.legend("Average")
plt.plot(data.T[0]/files_found,data.T[1]/files_found)
plt.show()
