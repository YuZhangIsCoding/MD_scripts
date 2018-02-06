#!/usr/bin/python
# Filename: gen_ndoped_quaternary.py
# This is a script rewrite to generate the .gro file and .itp file for ndoped graphene layers.
# The type of N here is quaternary, which meeans N substitude the C in the hexagonal ring.
# Date created: 06-10-2015

# dimensions
import pdb
import numpy as np

# Modify the percentage here
unit = [2, 4]
repeat = [7, 6]
l = unit[0]*repeat[0]
h = unit[1]*repeat[1]

# Some constants
bondlen = 0.142 
area = l*3/2*h/2*3**0.5*bondlen**2
e0 = 0.160217657
zgap = 4.53

# Inputs
sigma = input('Please input the surface charge(C/m^2):\n')

delta = sigma*area/e0/(unit[0]*unit[1]*repeat[0]*repeat[1])

# Define the functionalized spot
Nlist = [[0, 1]]
CNlist = [[0, 0], [0, 2], [1, 1]]

grofile = open('N_quaternary.gro', 'w')
itpfile = open('ngraphene.itp', 'w')

grofile.write('Quaternary N\n')
grofile.write(str(unit[0]*unit[1]*repeat[0]*repeat[1]*2)+'\n')

itpfile.write('[ moleculetype ]\n; molname\tnrexcl\n%s\t\t\t1\n' %('GPH'))
itpfile.write('[ atoms ]\n')
itpfile.write('; id type    res residu  at  cg  charge   mass\n')

natom = 0
for seq in [0, 1]:
    for i in range(0,l):
        for j in range(0,h):
            if i%2 == 0:
                coord = [1.5*i*bondlen+bondlen/2*(j%2), bondlen/2*3**0.5*j, seq*zgap]
            else:
                coord = [0.5*bondlen+1.5*i*bondlen-bondlen/2*(j%2), bondlen/2*3**0.5*j, seq*zgap]
            if [i%unit[0], j%unit[1]] in Nlist:
                atomname = 'NG'
                atomtype = 'grN'
                charge = delta*(-1)**seq-0.405
                molmass = 14.0067
            else:
                atomname = 'CG'
                atomtype = 'grC'
                molmass = 12.011
                if [i%unit[0], j%unit[1]] in CNlist:
                    charge = delta*(-1)**seq+0.135
                else:
                    charge = delta*(-1)**seq
            natom += 1
            grofile.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" %(1, 'GPH', atomname, natom,coord[0], coord[1], coord[2]))
            itpfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(natom, atomtype, 1, 'GPH', atomname, natom, charge, molmass ))
x=l*3/2*bondlen
y=h/2*3**0.5*bondlen
z=(x+y)/2

grofile.write('%10.5f%10.5f%10.5f' %(x, y, z))

grofile.close()
itpfile.close()
