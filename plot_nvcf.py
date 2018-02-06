#!/Users/yuzhang/anaconda/bin/python

import matplotlib.pyplot as plt
import numpy as np
import pdb

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

diff = np.loadtxt('nvcf.txt')

lname = ['OMIm', 'TFSI']
cl = ['r', 'b', 'g']
for i in range(1, len(diff[0])):
    plt.plot(diff[:, 0], diff[:, i], color = cl[i-1], lw = lwidth, label = lname[i-1])

plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel('Time (ps)', fontsize = fsize)
plt.ylabel('Normalized Velocity Autocorrelation Functions', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    axes.tick_params(which = 'minor', length = 4, width = 1)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.grid()
plt.tight_layout()
plt.savefig('nvcf_out.pdf')

fig = plt.figure(figsize = (8,6), dpi = 300)
for i in range(1, len(diff[0])):
    plt.plot(diff[:, 0], diff[:, i], color = cl[i-1], lw = lwidth)
plt.xlim([0, 2])
plt.ylim([-0.3, 0.2])
plt.yticks(np.arange(-0.3, 0.21, 0.1))
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    axes.tick_params(which = 'minor', length = 4, width = 1)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.grid()
plt.tight_layout()
plt.savefig('nvcf_out_small.pdf')
