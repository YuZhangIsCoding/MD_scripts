#!/Users/yuzhang/anaconda/bin/python

import numpy as np
import os
import matplotlib.pyplot as plt

count = 0
for filename in os.listdir("./"):
    if 'msd_' in filename and filename.endswith('.txt') and '_avg.txt' not in filename:
        if 'msd' in locals():
            temp = np.loadtxt(filename)
            if len(temp[0]) == len(msd[0]):
                msd[:, 1:] += temp[:, 1:]
                count = [i+1 for i in count]
            elif len(temp[0]) > len(msd[0]):
                print 'Previous length:', len(msd[0]),
                print 'Current length:', len(temp[0])
                count = [i+1 for i in count]+[1 for _ in range(len(temp[0])-len(msd[0]))]
                msd[:, 1:] += temp[:, 1:len(msd[0])]
                msd = np.concatenate((msd, temp[:, len(msd[0]):]), axis = 1)
            else:
                print 'Previous length:', len(msd[0]),
                print 'Current length:', len(temp[0])
                msd[:, 1:len(temp[0])] += temp[:, 1:]
                for i in range(len(temp[0])-1):
                    count[i] += 1
        else:
            msd = np.loadtxt(filename)
            count = [1 for _ in range(len(msd[0])-1)]
msd[:, 1:] /= np.array(count)

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

colors = ['red', 'blue', 'green', 'pink']

pt_start = len(msd)/2
for i in range(1, len(msd[0])):
    plt.plot(msd[:, 0], msd[:, i], linewidth = lwidth, color = colors[i-1], label = 'mol '+str(i))
    coef = np.polyfit(msd[pt_start:, 0], msd[pt_start:, i], 1)
    print coef*10000/4, ' E-10 m^2/s'
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.xlabel('Time (ps)', fontsize = fsize)
plt.ylabel('MSD ($\mathregular{nm^2\!}$)', fontsize = fsize)
plt.legend(loc = 'best', fontsize = fsize)
plt.tight_layout()
plt.savefig('msd_avg.pdf')
np.savetxt('msd_avg.txt', msd, fmt = ['%12.6f' for i in msd[0]])
