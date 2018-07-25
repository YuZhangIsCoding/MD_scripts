#!/Users/yuzhang/anaconda3/bin/python
# Filename:     expand_gro.py
# Description:  This is a python script that used to expand the gro file
#               according to its periodic boundary conditions.
# Date: 06-14-2018  Created
#       07-04-2018  Added support for triclinic box

import argparse, itertools
import numpy as np
from gro_common import Gro

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help = 'input gro file')
parser.add_argument('-o', '--output', default = 'out.gro', help = 'output gro file')
parser.add_argument('-expand', nargs = 3, default = (1, 1, 1), 
                    type = int, help = 'expand rate in x and y direction')

args = parser.parse_args()

gro = Gro()
gro.read(args.input)

outfile = open(args.output, 'w')
outfile.write('Expand from %s\n' %args.input)
outfile.write('%5d\n' %(gro.num*np.prod(args.expand)))
for i, j, k in itertools.product(*(range(_) for _ in args.expand)):
    for item in gro.info:
        outfile.write('%s%8.3f%8.3f%8.3f\n' %(item[0], 
                item[1][0]+i*gro.box[0]+j*gro.box[5]+k*gro.box[7], 
                item[1][1]+j*gro.box[1]+i*gro.box[3]+k*gro.box[8], 
                item[1][2]+k*gro.box[2]+i*gro.box[4]+j*gro.box[6]))
outfile.write('%10.5f%10.5f%10.5f\n' %tuple([args.expand[i]*gro.box[i] for i in range(3)]))
outfile.close()
