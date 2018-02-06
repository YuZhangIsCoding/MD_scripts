#!/Users/yuzhang/anaconda/bin/python
import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'Input the dimensions')
parser.add_argument('-d', '--dimension', dest = 'dim', nargs = '*', type = float, help = 'The dimension for the system')
parser.add_argument('-s', '--size', dest = 'size', nargs = '*', type = float, help = 'The stepsize of each bin')
args = parser.parse_args()

dim = [int(np.ceil(args.dim[i]*int(1/args.size[i]))) for i in range(3)]

gridfile = open('gridElectrode1.dat', 'w')
gridfile.write('%5d\n' %np.prod(dim))
for i in range(dim[0]):
    for j in range(dim[1]):
        for k in range(dim[2]):
            gridfile.write('%12.4f%12.4f%12.4f\n' %(i*args.size[0], j*args.size[1], k*args.size[2]))
gridfile.close()            
