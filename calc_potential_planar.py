#!/Users/yuzhang/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt
import argparse, sys, pdb

parser = argparse.ArgumentParser(description = 'User specified filenames, boundaries, etc.')
parser.add_argument('-b', '--begin', dest = 'begin', default = 0, type = float, help = 'input the first frame')
parser.add_argument('-d', '--dimension', dest = 'dim', nargs = '*', type = float, help = 'The dimension for the system')
parser.add_argument('-s', '--size', dest = 'size', type = float, nargs = '*', help = 'The stepsize for each bing')
parser.add_argument('-sc', '--surfcharge', dest = 'sc', type = float, help = 'Surface charge density')
parser.add_argument('-dist', '--distance', dest = 'dist', type = float, help = 'The distance between 2 plates')
args = parser.parse_args()

fsize = 28
lwidth = 4.0
msize = 16.0

dim = [int(np.ceil(args.dim[i]*int(1/args.size[i]))) for i in range(3)]

myfile = open('potential_surf_kr.dat', 'r')
chunk = 0

bulk_avg = []
flag = 0
for line in myfile:
    if line == '\n':
        count = 0
        if flag == 1:
            bulk_avg.append([np.mean(item) for item in bulk])
            flag = 0
        bulk = [[] for row in range(dim[2])]
    elif 'gfeng' in line:
        chunk += 1
    elif chunk >= args.begin:
        temp = [float(_) for _ in line.split()]
        flag = 1
        count += 1
        bulk[count%dim[2]].append(sum(temp))
print 'Total chunks:', chunk

bulk_avg = np.array(bulk_avg).T
bulk = np.mean(bulk_avg, axis = 1)

dens = np.loadtxt('Density.txt')
if args.sc == None:
    sc = float(raw_input('Please specify the surface charge density (C/m^2):\n'))
else:
    sc = args.sc
if args.dist == None:
    dist = float(raw_input('Please specify the distance between 2 charged planes (nm):\n'))
else:
    dist = args.dist

Q1 = -sc
Q2= sum(dens[:, 0]*dens[:, 2]*(dens[1, 0]-dens[0, 0]))/dist
const = 1000/8.85418784762/6.242/0.160217657

print Q1, Q2

for i, item in enumerate(bulk):
    bulk[i] = item*0.010364272+(Q1)*i*args.size[2]*const/dist

fig = plt.figure(figsize = (16, 12), dpi = 300)
plt.plot(range(len(bulk)), bulk*0.010364272, linewidth = lwidth, label = 'Potential in the bulk')
plt.legend(fontsize = fsize)

for ax in fig.axes:
    ax.set_ylabel('Potential (V)', fontsize = fsize)
    ax.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        ax.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('potential_fly.pdf')

np.savetxt('pd_fly.txt', np.column_stack((range(len(bulk)), bulk)),fmt=['%12.4f','%12.4f'])
