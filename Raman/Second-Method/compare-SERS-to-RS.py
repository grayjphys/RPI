import numpy as np
import matplotlib.pyplot as plt

dataRS = np.loadtxt('POL-X-PROP-Z/spec.dat')
plt.plot(dataRS[:,0],dataRS[:,1])

dataRS = np.loadtxt('Efield-1nm/POL-X-PROP-Z/spec.dat')
plt.plot(dataRS[:,0],dataRS[:,1])

plt.show()