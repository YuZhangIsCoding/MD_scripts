import matplotlib.pyplot as plt
import numpy as np
import pdb

fig = plt.figure(figsize = (16,12), dpi = 1000)
fsize = 36
lwidth = 4.0
msize = 16.0

slit = np.loadtxt('ND_slit.txt')
ohslit = np.loadtxt('ND_ohslit.txt')

plt.plot(slit[:,0], slit[:,1], linewidth = lwidth, color = 'r', label = '$\mathregular{[EMIm^+]}$: Defunctionalized')
plt.plot(slit[:,0], slit[:,2], linewidth = lwidth, color = 'b', label = '$\mathregular{[TFSI^-]}$: Defunctionalized')
plt.plot(slit[:,0], ohslit[:,1], '--', linewidth = lwidth, color = 'r', label = '$\mathregular{[EMIm^+]}$: Oxidized')
plt.plot(slit[:,0], ohslit[:,2], '--', linewidth = lwidth, color = 'b', label = '$\mathregular{[TFSI^-]}$: Oxidized')
plt.legend(fontsize = 28, frameon = False)
plt.xticks(fontsize = fsize)
plt.yticks(fontsize = fsize)
plt.xlim([0.3, 0.8])
plt.xlabel('Distance (nm)', fontsize = fsize)
plt.ylabel('Number density ($\mathregular{\#/nm^3}$)', fontsize = fsize)
plt.gca().tick_params(labelsize = fsize, width = 2, length = 6)
for direction in ['top', 'bottom', 'left', 'right']:
    plt.gca().spines[direction].set_linewidth(2)

plt.tight_layout()
plt.savefig('ND_compare.tiff', format = 'tiff')
