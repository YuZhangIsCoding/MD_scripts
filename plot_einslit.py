#!/Users/yuzhang/anaconda/bin/python
import matplotlib.pyplot as plt
import numpy as np
import pdb

fig, (ax2, ax1) = plt.subplots(2, figsize = (12,12), dpi = 1000)
fsize = 28
lwidth = 4.0
msize = 16.0
data = np.loadtxt('E_all.txt')

data[:, 1] -= 32.73
data[:, 3] += 32.73
data[:, 1] -= 30.03
data[:, 3] += 30.03

sep = 2000

fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (16, 16), dpi = 1000)
region = ['slit 1', 'bulk 1', 'slit 2', 'bulk 2']
for i in range(2):
    for j in range(2):
        c = 2*i+j
        axes[i, j].plot(data[:sep, 0], data[:sep, c+1], 'k-', linewidth = lwidth, label = 'Charge - '+region[c])
#axes[0, 0].plot(data[0, 0], data[0, 1]-32.73*2, 's', markeredgecolor = 'blue', markerfacecolor = 'white', markersize = msize, markeredgewidth = 2)
#axes[1, 0].plot(data[0, 0], data[0, 3]+32.73*2, 's', markeredgecolor = 'red', markerfacecolor = 'white', markersize = msize, markeredgewidth = 2)
axes[0, 0].plot(data[0, 0], data[0, 1]+30.03*2, 's', markeredgecolor = 'blue', markerfacecolor = 'white', markersize = msize, markeredgewidth = 2)
axes[1, 0].plot(data[0, 0], data[0, 3]-30.03*2, 's', markeredgecolor = 'red', markerfacecolor = 'white', markersize = msize, markeredgewidth = 2)
axes[0, 1].plot(data[0, 0], data[0, 2], 's', markeredgecolor = 'black', markerfacecolor = 'white', markersize = msize, markeredgewidth = 2)
axes[1, 1].plot(data[0, 0], data[0, 4], 's', markeredgecolor = 'black', markerfacecolor = 'white', markersize = msize, markeredgewidth = 2)
axes[1, 0].set_xlabel('Time (ns)', fontsize = fsize)
axes[1, 1].set_xlabel('Time (ns)', fontsize = fsize)
axes[0, 0].set_ylabel('Charge in the slit (e)', fontsize = fsize)
axes[1, 0].set_ylabel('Charge in the slit (e)', fontsize = fsize)
for axes in fig.axes:
    axes.set_xlim(-data[sep, 0]/50)
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
    axes.legend(loc = 'best', fontsize = 24)
    axes.grid()
plt.tight_layout()
plt.savefig('E_all_short.pdf')
