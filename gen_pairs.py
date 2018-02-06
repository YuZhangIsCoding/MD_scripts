#!/Users/yuzhang/anaconda/bin/python
# Filename: gen_pairs.py
# Description: This is a python script to generate the cross LJ interaction parameters, knowing the parameters for the pure atoms
# Date: 02-24-2016 Created

import numpy as np
import argparse, os, sys

parser = argparse.ArgumentParser(description = 'User specificed filenames')
parser.add_argument('-i', '--input', dest = 'file_name', default = 'sig.txt', help = 'Filename for the nonbond parameter for pure atom')
args = parser.parse_args()

if not os.path.isfile(args.file_name):
    sys.exit('Exit: not file named %s found in current directory' %args.file_name)
file_in = open(args.file_name)

data = []
for line in file_in:
    if line[0] != ';':
        data.append(line.split())
file_in.close()
for i in range(len(data)):
    for j in range(i, len(data)):
        data[i][-1] = float(data[i][-1])
        data[j][-1] = float(data[j][-1])
        data[i][-2] = float(data[i][-2])
        data[j][-2] = float(data[j][-2])
        print '%6s%6s%5d%10.4f%10.4f' %(data[i][0], data[j][0], 1, (data[i][-2]+data[j][-2])/2, (data[i][-1]*data[j][-1])**0.5)
        #if data[i][0] == 'AU' or data[j][0] == 'AU':
        #    print '%6s%6s%5d%10.4f%10.4f' %(data[i][0], data[j][0], 1, (data[i][-2]+data[j][-2])/2, (data[i][-1]*data[j][-1])**0.5)
