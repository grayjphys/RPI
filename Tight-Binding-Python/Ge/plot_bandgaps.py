import matplotlib.pyplot as plt

compositions = [0.0,   0.125, 0.25, 0.25,  0.375, 0.375, 0.5,   0.5,   0.5,  0.5,  0.5,   0.625, 0.625, 0.75, 0.75,  0.875, 1.0]
band_gaps =    [0.825, 0.627, 0.58, 0.629, 0.6,   0.746, 0.607, 0.617, 0.62, 0.62, 0.837, 0.48,  0.6,   0.32, 0.338, 0.115, 0.491]
space_groups = [227,215,115,166,35,215,160,51,91,95,216,215,35,166,115,215,227]

plt.title('Band Gap vs. Composition in ' + r'$Si_{1-x}Ge_x$',fontsize=30)
plt.yticks(fontsize=20)
plt.xticks(compositions,fontsize=20)
plt.ylabel('Band Gap (eV)',fontsize=20)
plt.xlabel('Fractional Composition (x)',fontsize=20)

for i in range(len(compositions)):

	plt.scatter(compositions[i],band_gaps[i])
	space_groups[i] = '(' +str(compositions[i]) + ', ' + str(space_groups[i])+')'
plt.legend(title='X, Space Groups',labels=space_groups,ncol=2,fontsize=20)
plt.show()