#!/usr/bin/python
# Filename: convert_gro_itp.py
# This is a python script to generate .itp file for n-doped graphene systems

import sys

myfile = open(sys.argv[1], 'r')
mylines = myfile.readlines()
boxsize = mylines[-1].split()
area = float(boxsize[0])*float(boxsize[1])
sigma = input('Please input the surface charge density (C/m^2):\n')
e0 = 0.160217657
delta = sigma*area/e0/(len(mylines)-3)

itpfile = open('npyrrolic.itp','w')
itpfile.write('[ moleculetype ]\n; molname   nrexcl\nGPH         1\n[ atoms ]\n')
itpfile.write('; id type    res residu  at  cg  charge   mass\n')
for seq in [0, 1]:
    for count in range(len(mylines)):
        if count > 1 and count < len(mylines)-1: 
                n_atom = count-1+seq*(len(mylines)-3)
                n_judge = n_atom%34
                if n_judge == 1:
                    atomtype = 'HN'
                    atomname = 'HN'
                    charge = 0.35+delta*(-1)**seq
                    molmass = 1.008
                elif n_judge in [2, 3]:
                    atomtype = 'HC'
                    atomname = 'HC'
                    charge = 0.207+delta*(-1)**seq
                    molmass = 1.008
                elif n_judge == 0:
                    atomtype = 'grN'
                    atomname = 'NC'
                    charge = -0.53+delta*(-1)**seq
                    molmass = 14.0067
                elif n_judge in [17, 24]:
                    atomtype = 'grC'
                    atomname = 'CN'
                    charge = 0.09+delta*(-1)**seq
                    molmass = 12.0110
                elif n_judge in [23, 30]:
                    atomtype = 'grC'
                    atomname = 'CH'
                    charge = -0.207+delta*(-1)**seq
                    molmass = 12.011
                else:
                    atomtype = 'grC'
                    atomname = 'CG'
                    charge = delta*(-1)**seq
                    molmass = 12.011
                itpfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(n_atom, atomtype, 1, 'GPH', atomname, n_atom, charge, molmass))
myfile.close()
itpfile.close()
