#!/Users/yuzhang/anaconda/envs/py3/bin/python

import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'Input the dimensions')
parser.add_argument('-d', '--dimension', dest = 'dim', nargs = '*', type = float, help = 'The dimension for the system')
args = parser.parse_args()

size = 0.1
gridfile = open('gridElectrode1.dat', 'w')
gridfile.write('%5d\n' %(int(args.dim[1]*10+1)*(4*(int(args.dim[0]*5-110))+int(args.dim[0]*10))))
for seq in range(2):
    for i in np.arange(20+seq*args.dim[0]*5, args.dim[0]*5-70-20+seq*args.dim[0]*5):
        for j in np.arange(args.dim[1]*10):
            for k in [1.8*10, (args.dim[2]-1.8)*10]:
                gridfile.write('%12.4f%12.4f%12.4f\n' %(i*size, j*size, k*size))
for i in np.arange(args.dim[0]*10):
    for j in np.arange(args.dim[1]*10):
        gridfile.write('%12.4f%12.4f%12.4f\n' %(i*size, j*size, args.dim[2]/2))
gridfile.close()
