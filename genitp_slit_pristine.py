#!/Users/yuzhang/anaconda/bin/python
# Filename: gen_slit_pore_functionalized.py
# This script is to generate functionlized slit pore systems. The outer-most layer has dimensions of 28*28*16,
# and the smallest unit cell is 4*4, so the percentage of the functional group (here just use hydroxyl group) is 6.25% or several times it.
# Date: 6-1-2015

import numpy as np
import pdb

############ Functions ###########

def gen_frame(l,h,stype):
    '''General frame to generate a graphene sheet'''
    gframe = [[] for row in range(l*h)]
    for i in range(0,l):
        for j in range(0,h):
            seq = i*h+j
            if stype == 1:
                if i%2 == 0:
                    gframe[seq] = np.array([1.5*i+(j%2)/2.0, j*3**0.5/2, 0])*bondlen
                else:
                    gframe[seq] = np.array([0.5+1.5*i-(j%2)/2.0, j*3**0.5/2, 0])*bondlen
            elif stype == 2:
                if i%2 == 0:
                    gframe[seq] = np.array([1.5*i-(j%2)/2.0, j*3**0.5/2, 0])*bondlen
                else:
                    gframe[seq] = np.array([-0.5+1.5*i+(j%2)/2.0, j*3**0.5/2, 0])*bondlen
    return gframe

def gen_functionalized(l, h):
    '''Generate frame for functionalized graphene'''
    bond_ccp = 0.151
    bond_coh = 0.1425
    bond_ohho = 0.096
    angle_coh = 109.5
    ccpdz = (bond_ccp**2-bondlen**2)**0.5
    ohdx = bond_ohho*np.sin(np.radians(angle_coh))
    ohdz = -bond_ohho*np.cos(np.radians(angle_coh))
    gframe = gen_frame(l, h, 1)
    fframe = []
    fspot = []
    for i in range(0, l):
        for j in range(0, h):
            seq = i*h+j
            if i%unit[0] == 1 and j%unit[0] == 1:
                fspot.append(seq)
                unit_num = i/unit[0]*h/unit[1]+j/unit[1]
                gframe[seq][2] = ccpdz*(-1)**unit_num
                fframe.append(np.array([gframe[seq][0], gframe[seq][1], (ccpdz+bond_coh)*(-1)**unit_num]))
                fframe.append(np.array([gframe[seq][0]+ohdx, gframe[seq][1], (ccpdz+bond_coh+ohdz)*(-1)**unit_num]))
    return gframe, fframe, fspot
    

def gen_circle(pstart, lx, ly, lz, myfile, gtype):
    '''Generate a circle of graphene sheets'''
    print 'Input circle size:',lx,ly,lz
    nx = int(np.floor(lx/bondlen/3*2))+1
    ny = int(np.ceil(ly/bondlen/3**0.5))*2
    nz = int(np.floor(lz/bondlen/3*2))
    # Button and top sheets
    lfspot = []
    vfspot = []
    while True:
        if gtype == 1 or gtype == 3:
            lframe = gen_frame(nx-2, ny, 1)
            if nx%2 == 0:
                x_real = bondlen+max(lframe[-1][0],lframe[-2][0])+bondlen
            else:
                x_real = bondlen+min(lframe[-1][0],lframe[-2][0])+bondlen*2
        elif gtype == 2:
            lframe, lfframe, lfspot = gen_functionalized(nx-2, ny)
            if nx%2 == 0:
                x_real = bondlen+max(lframe[-1][0],lframe[-2][0])+bondlen
            else:
                x_real = bondlen+min(lframe[-1][0],lframe[-2][0])+bondlen*2
        else:
            print 'Functionalization type unknown'
        dx = lx-x_real
        print '\tdx',dx
        if dx > 2*bondlen:
            nx += 1
        else:
            break
    # Left and right sheets
    while True:
        if gtype == 1 or gtype == 3:
            vframe = gen_frame(nz, ny, 1)
            lz_real = max(vframe[-1][0],vframe[-2][0])
        elif gtype == 2:
            vframe, vfframe, vfspot = gen_functionalized(nz, ny)
            lz_real = max(vframe[-1][0],vframe[-2][0])
        else:
            print 'Functionalization type unknown'
        dz = lz-lz_real
        print '\tdz',dz
        if dz > 2*bondlen:
            nz += 1
        else:
            break
    if gtype == 3:
        e0 = 0.160217657
        delta = surc*(ny*nz*bondlen**2*3*3**0.5/4)/e0/(ny*nz)
    elif gtype == 2:
        e0 = 0.160217657
        delta = surc*(ny*nz*bondlen**2*3*3**0.5/4)/e0/(ny*nz+2*len(vfspot))
    else:
        delta = 0
    if 'atom_num' not in globals():
        global atom_num
        atom_num = 0
    fspot = []
    for rept in range(2):
        for i, item in enumerate(lframe):
            atom_num += 1
            if i in lfspot:
                fspot.append(atom_num)
                myfile.write('%6d%6s%6d%6s%6s%6d%12.7f%8.3f\n' %(atom_num, 'COH', 1, 'GPH', 'COH', atom_num, delta+0.373, 12.011))
            else:
                myfile.write('%6d%6s%6d%6s%6s%6d%12.7f%8.3f\n' %(atom_num, 'CG', 1, 'GPH', 'CG', atom_num, delta, 12.011))
    for rept in range(2):
        for i, item in enumerate(vframe):
            atom_num += 1
            if i in vfspot:
                fspot.append(atom_num)
                myfile.write('%6d%6s%6d%6s%6s%6d%12.7f%8.3f\n' %(atom_num, 'COH', 1, 'GPH', 'COH', atom_num, delta+0.373, 12.011))
            else:
                myfile.write('%6d%6s%6d%6s%6s%6d%12.7f%8.3f\n' %(atom_num, 'CG', 1, 'GPH', 'CG', atom_num, delta, 12.011))
    if gtype == 2:
        fgroup = []
        for rept in range(2):
            for i, item in enumerate(lfframe):
                atom_num += 1
                fgroup.append(atom_num)
                if i%2 == 0:
                    myfile.write('%6d%6s%6d%6s%6s%6d%12.7f%8.3f\n' %(atom_num, 'OG', 1, 'GPH', 'OG', atom_num, delta-0.789, 15.9994))
                else:
                    myfile.write('%6d%6s%6d%6s%6s%6d%12.7f%8.3f\n' %(atom_num, 'HG', 1, 'GPH', 'HG', atom_num, delta+0.416, 1.008))
        for rept in range(2):
            for i, item in enumerate(vfframe):
                atom_num += 1
                fgroup.append(atom_num)
                if i%2 == 0:
                    myfile.write('%6d%6s%6d%6s%6s%6d%12.7f%8.3f\n' %(atom_num, 'OG', 1, 'GPH', 'OG', atom_num, delta-0.789, 15.9994))
                else:
                    myfile.write('%6d%6s%6d%6s%6s%6d%12.7f%8.3f\n' %(atom_num, 'HG', 1, 'GPH', 'HG', atom_num, delta+0.416, 1.008))
        return fspot, fgroup, ny

def write_bonds(myfile, fspot, fgroup, ny):
    '''Function to write basic information about atoms'''
    myfile.write('[ bonds ]\n')
    for i, item in enumerate(fspot):
        myfile.write('%5d%5d%5d\n' %(item, item-ny, 1))
        myfile.write('%5d%5d%5d\n' %(item, item-1, 1))
        myfile.write('%5d%5d%5d\n' %(item, item+1, 1))
        myfile.write('%5d%5d%5d\n' %(item, fgroup[2*i], 1))
        myfile.write('%5d%5d%5d\n' %(fgroup[2*i], fgroup[2*i+1], 1))

def write_angles(myfile, fspot, fgroup, ny):
    '''Function to write angle information about atoms'''
    myfile.write('[ angles ]\n')
    for i, item in enumerate(fspot):
        myfile.write('%5d%5d%5d%5d\n' %(item-ny, item, fgroup[2*i], 1))
        myfile.write('%5d%5d%5d%5d\n' %(item-1, item, fgroup[2*i], 1))
        myfile.write('%5d%5d%5d%5d\n' %(item+1, item, fgroup[2*i], 1))
        myfile.write('%5d%5d%5d%5d\n' %(item, fgroup[2*i], fgroup[2*i+1], 1))

def write_dihedrals(myfile, fspot, fgroup, ny):
    '''Function to write dihedral information about atoms'''
    myfile.write('[ dihedrals ]\n')
    for i, item in enumerate(fspot):
        myfile.write('%5d%5d%5d%5d%5d\n' %(item-ny, item, fgroup[2*i], fgroup[2*i+1], 1))
        myfile.write('%5d%5d%5d%5d%5d\n' %(item-1, item, fgroup[2*i], fgroup[2*i+1], 1))
        myfile.write('%5d%5d%5d%5d%5d\n' %(item+1, item, fgroup[2*i], fgroup[2*i+1], 1))

################### Main #####################

#### Inputs ####
global bondlen
global unit
global atom_num
global surc
bondlen = 0.142 # Specify the bond length for graphene
unit = [4, 4] # Repeat unit for functional groups
surc = input('Please specify the surface charge density (e/nm^2):\n')

# Lengths for the out circle
lx = 10
ly = 4
lz = 3.5

for elect in ['NEG', 'POS']:
    atom_num = 0
    myfile = open('%s.itp' %elect, 'w')
    myfile.write('[ moleculetype ]\n; molname   nrexcl\n')
    myfile.write('%s   3\n' %elect)
    myfile.write('[ atoms ]\n')
    myfile.write('; id type    res residu  at  cg  charge   mass\n')
    
    surc *= -1
    gen_circle([0, 0], lx, ly, lz, myfile, 3)
    gen_circle([0.341, 0.341], lx-0.682, ly, lz-0.682, myfile, 1)
    gen_circle([0.682, 0.682], lx-0.682*2, ly, lz-0.682*2, myfile, 1)
    #write_bonds(myfile, fspot, fgroup, ny)
    #write_angles(myfile, fspot, fgroup, ny)
    #write_dihedrals(myfile, fspot, fgroup, ny)
    myfile.close()
