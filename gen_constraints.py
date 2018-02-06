#!/usr/bin/python/
# Filename: gen_constraints.py
# Description: This is a python script to generate the constraints for gold atoms
# Date: 10-22-2015 Created

import sys
import numpy as np

if '-i' in sys.argv:
    temp = sys.argv.index('-i')+1
    groname = sys.argv[temp]
else:
    groname = raw_input('Please input the configuration file (.gro):\n')

grofile = open(groname, 'r')
coord = []
for i,line in enumerate(grofile):
    if i >= 2:
        entry = line.split()
        if len(entry) > 3:
            coord.append([float(entry[-3]), float(entry[-2]), float(entry[-1])])
grofile.close()

outfile = open('AU.itp','w')
outfile.write('[ moleculetype ]\nAU\t\t1\n[atoms]\n')
for i in range(len(coord)):
    outfile.write('%5d\tAU\t1\tAU\tAU%5d%10.4f%10.4f\n' %(i+1, i+1, 0.0, 196.967))

outfile.write('[ bonds ]\n')

dist_jdg = 0.288
dist0 = (dist_jdg+0.01)**2
coord = np.array(coord)
count = 0
for i in range(len(coord)):
    for j in range(i+1, len(coord)):
        dist = sum((coord[i]-coord[j])**2)
        if dist < dist0:
            count += 1
            #outfile.write('%5d%5d  1 %10.4f\n' %(i+1, j+1, dist_jdg))
            outfile.write('%5d%5d  1 %10.4f\t%10.4f\n' %(i+1, j+1, dist_jdg, 1.0e+07))
outfile.close()
print count
