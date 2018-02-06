#!/Users/yuzhang/anaconda/bin/python
import matplotlib.pyplot as plt
import numpy as np
import pdb

msd_x = np.loadtxt('msd_x.txt')
msd_y = np.loadtxt('msd_y.txt')

msd = msd_x+msd_y
msd[:, 0] /= 2

np.savetxt('msd_combined.txt', msd, fmt = ['%12.6f' for i in range(len(msd[0]))])

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

for i in range(1, len(msd[0])):
    plt.plot(msd[:, 0]/1000, msd[:, i], linewidth = lwidth, label = 'mol '+str(i))

plt.xlabel('t (ns)', fontsize = fsize)
plt.ylabel('MSD (nm^2)', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.legend(loc = 'best', fontsize = fsize)
plt.tight_layout()
plt.savefig('msd_combined.pdf')
