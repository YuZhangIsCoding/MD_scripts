#!/Users/yuzhang/anaconda/bin/python
# Filename:     combine_nvcf.py
# Description:  This is a python script to combine NVACFs and VACFs with different time steps
# Date:     04-05-2017  Created

import matplotlib.pyplot as plt
import numpy as np
import pdb, sys, os, argparse

fig = plt.figure(figsize = (16,12), dpi = 300)
fsize = 28
lwidth = 4.0
msize = 16.0

parser = argparse.ArgumentParser(description = 'User specified filenames')
parser.add_argument('-i', '--input', dest = 'filename', nargs = '*', help = 'Specify the NVACF files you want to combine')
args = parser.parse_args()

if args.filename == []:
    sys.exit('Exit: no file specified!')
else:
    raw = []
    for myfile in args.filename:
        if not os.path.isfile(myfile):
            print 'File', myfile, 'not found in current directory!'
            continue
        else:
            raw.append(np.loadtxt(myfile))
    if raw == []:
        sys.exit('Exit: no files found!')

mt = [item[-1, 0] for item in raw]
seq = [mt.index(i) for i in sorted(mt)]

data = raw[seq[0]]
for i in seq[1:]:
    for j in range(len(raw[i])):
        if raw[i][j, 0] > data[-1, 0]:
            data = np.concatenate((data, raw[i][j:]), axis = 0)
            break

np.savetxt('nvcf_combined.txt', data, fmt = ['%12.6f' for i in range(len(data[0]))])

mycmp = plt.cm.gist_ncar
plt.gca().set_color_cycle([mycmp(i) for i in np.linspace(0, 0.9, len(data[0])-1)])
for i in range(1, len(data[0])):
    plt.plot(data[:, 0], data[:, i], linewidth = lwidth, label = i)
plt.legend(loc = 'best', fontsize = fsize, frameon = True, numpoints = 1)
plt.xlabel('Time (ps)', fontsize = fsize)
plt.ylabel('NVACF', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    axes.tick_params(which = 'minor', length = 4, width = 1)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('nvcf_combined.pdf')
