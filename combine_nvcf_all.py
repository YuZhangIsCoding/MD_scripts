#!/Users/yuzhang/anaconda/bin/python

import numpy as np
import os
import matplotlib.pyplot as plt

count = 0
for filename in os.listdir("./"):
    if 'nvcf_' in filename and filename.endswith('.txt') and '_avg' not in filename:
        count += 1
        if 'nvcf' in locals():
            nvcf[:, 1:] += np.loadtxt(filename)[:, 1:]
        else:
            nvcf = np.loadtxt(filename)
nvcf[:, 1:] /= count

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

colors = ['red', 'blue', 'green', 'pink']

for i in range(1, len(nvcf[0])):
    plt.plot(nvcf[:, 0], nvcf[:, i], linewidth = lwidth, color = colors[i-1], label = 'mol '+str(i))
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.xlabel('Time (ps)', fontsize = fsize)
plt.ylabel('Velocity Autocorrelation Function ($\mathregular{m^2\!/s^2\!}$)', fontsize = fsize)
plt.legend(loc = 'best', fontsize = fsize)
plt.tight_layout()
plt.savefig('nvcf_avg.pdf')
np.savetxt('nvcf_avg.txt', nvcf, fmt = ['%14.6E' for i in nvcf[0]])
