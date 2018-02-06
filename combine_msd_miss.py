#!/Users/yuzhang/anaconda/bin/python
# Filename : combine_msd_miss.py
# Description:  this is a python script that combines several different msd_*.txt profile into one msd_avg.txt profile
#               However, there might be cases that some columns are missing in some files. This script combines msd files
#               according to the file lengths and the sequence of potential missing colomn from input.
# Date:     04-26-2017  Created

import numpy as np
import os, argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = 'Specify missing colomns')
parser.add_argument('-m', '--miss', dest = 'miss', type = int, nargs = '*', help = 'colums that might be missing, please input according to their possibilities of missing, from high to low')
parser.add_argument('-s', '--size', dest = 'size', type = int, help = 'output columns')
args = parser.parse_args()
if args.miss == []:
    miss = [int(_) for _ in raw_input('Please input the columns that might be missing:\n').split()]
else:
    miss = args.miss
if args.size == None:
    size = int(raw_input('Please input the output column number:\n'))
else:
    size = args.size
count = [0 for i in range(size)]
for filename in os.listdir("./"):
    if 'msd_' in filename and filename.endswith('.txt') and '_avg.txt' not in filename:
        temp = np.loadtxt(filename)
        diff = size-len(temp[0])+1
        if 'msd' not in locals():
            msd = np.zeros((len(temp), size+1))
            msd[:, 0] = temp[:, 0]
        if diff > 0:
            subcount = 0
            for i in range(1, size+1):
                if i not in miss[:diff]:
                    count[i-1] += 1
                    subcount += 1
                    msd[:, i] += temp[:, subcount]
        else:
            count = [i+1 for i in count]
            msd[:, 1:] += temp[:, 1:]
msd[:, 1:] /= np.array(count)
#import pdb
#pdb.set_trace()
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
