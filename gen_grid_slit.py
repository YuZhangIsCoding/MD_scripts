#!/Users/yuzhang/anaconda3/bin/python
# Filename:     gen_grid_slit.py
# Description:  Generate grids on the slit surfaces

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-x', nargs = '*', type = float, help = 'x range')
parser.add_argument('-y', nargs = '*', type = float, help = 'y range')
parser.add_argument('-z', nargs = '*', type = float, help = 'z range')
parser.add_argument('-bs', type = float, default = 0.1, help = 'bin size of the grid')
parser.add_argument('-switch', action = 'store_true', help = 'swith x and z')
parser.add_argument('-o', default = 'gridElectrode1.dat', help = 'output filename')
args = parser.parse_args()
gridfile = open(args.o, 'w')
for iind in range(len(args.x)//2):
    for i in np.arange(np.round(int(args.x[2*iind]/args.bs)*args.bs, decimals = 1), np.round(int(args.x[2*iind+1]/args.bs)*args.bs, decimals = 1), args.bs):
        for k in args.z:
            for j in np.arange(int(args.y[0]/args.bs)*args.bs, int(args.y[1]/args.bs)*args.bs, args.bs):
                if args.switch:
                    gridfile.write('%12.4f%12.4f%12.4f\n' %(k, j, i))
                else:
                    gridfile.write('%12.4f%12.4f%12.4f\n' %(i, j, k))
gridfile.close()
