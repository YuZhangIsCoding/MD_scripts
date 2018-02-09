#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: gen_slit_pore_functionalized.py
# This script is to generate functionlized slit pore systems. 
# The outer-most layer has dimensions of 28*28*16,
# and the smallest unit cell is 4*4, so the percentage of 
# the functional group (here just use hydroxyl group) is 6.25% or 
# several times of it.
# Date: 6-1-2015 Created

import numpy as np

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
    

def gen_circle(pstart, lx, ly, lz, myfile, gtype, symm_choice):
    '''Generate a circle of graphene sheets'''

    print('Input circle size:',lx,ly,lz)
    nx = int(np.floor(lx/bondlen/3*2))+1
    ny = int(np.ceil(ly/bondlen/3**0.5))*2
    nz = int(np.floor(lz/bondlen/3*2))
    
    # Button and top sheets
    lfspot = []
    vfspot = []
    while True:
        if gtype == 1:
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
            print('Functionalization type unknown')
        dx = lx-x_real
        print('\tdx',dx)
        if dx > 2*bondlen:
            nx += 1
        else:
            break
    # Left and right sheets
    while True:
        if gtype == 1:
            vframe = gen_frame(nz, ny, 1)
            lz_real = max(vframe[-1][0],vframe[-2][0])
        elif gtype == 2:
            vframe, vfframe, vfspot = gen_functionalized(nz, ny)
            lz_real = max(vframe[-1][0],vframe[-2][0])
        else:
            print('Functionalization type unknown')
        dz = lz-lz_real
        print('\tdz',dz)
        if dz > 2*bondlen:
            nz += 1
        else:
            break
    
    print('Real row and column numbers:', nx-2, ny, nz)
    for i, item in enumerate(lframe):
        if i in lfspot:
            myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'COH',1,pstart[0]+item[0]+bondlen+dx/2,item[1],pstart[1]+item[2]))
        else:
            myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'CG',1,pstart[0]+item[0]+bondlen+dx/2,item[1],pstart[1]+item[2]))
    for i, item in enumerate(lframe):
        if i in lfspot:
            myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'COH',1,pstart[0]+item[0]+bondlen+dx/2,item[1],lz+pstart[1]+item[2]*(-1)**symm_choice))
        else:
            myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'CG',1,pstart[0]+item[0]+bondlen+dx/2,item[1],lz+pstart[1]+item[2]*(-1)**symm_choice))
    for i, item in enumerate(vframe):
        if i in vfspot:
            myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'COH',1, pstart[0]+item[2], item[1], pstart[1]+dz/2+item[0]))
        else:
            myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'CG',1, pstart[0]+item[2], item[1], pstart[1]+dz/2+item[0]))
    for i, item in enumerate(vframe):
        if i in vfspot:
            myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'COH',1, lx+pstart[0]-item[2], item[1], pstart[1]+dz/2+item[0]))
        else:
            myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'CG',1, lx+pstart[0]-item[2], item[1], pstart[1]+dz/2+item[0]))
    if gtype == 2:
        for i, item in enumerate(lfframe):
            if i%2 == 0:
                myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'OG',1,pstart[0]+item[0]+bondlen+dx/2,item[1],pstart[1]+item[2]))
            else:
                myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'HG',1,pstart[0]+item[0]+bondlen+dx/2,item[1],pstart[1]+item[2]))
        for i, item in enumerate(lfframe):
            if i%2 == 0:
                myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'OG',1,pstart[0]+item[0]+bondlen+dx/2,item[1],lz+pstart[1]+item[2]*(-1)**symm_choice))
            else:
                myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'HG',1,pstart[0]+item[0]+bondlen+dx/2,item[1],lz+pstart[1]+item[2]*(-1)**symm_choice))
        for i, item in enumerate(vfframe):
            if i%2 == 0:
                myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'OG',1, pstart[0]+item[2], item[1], pstart[1]+dz/2+item[0]))
            else:
                myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'HG',1, pstart[0]+item[2], item[1], pstart[1]+dz/2+item[0]))
        for i, item in enumerate(vfframe):
            if i%2 == 0:
                myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'OG',1, lx+pstart[0]-item[2], item[1], pstart[1]+dz/2+item[0]))
            else:
                myfile.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(1,'GPH', 'HG',1, lx+pstart[0]-item[2], item[1], pstart[1]+dz/2+item[0]))
    return lframe[-1][1]+3**0.5/2*bondlen

#################### Main #####################

#### Inputs ####
global bondlen
global unit
bondlen = 0.142 # Specify the bond length for graphene
unit = [4, 4] # Repeat unit for functional groups

# Lengths for the out circle
symm_choice = input('Please select the symmetry of the slit pore:\n\
                    1. Symmetric\n\
                    2. Asymmetric\n')
lx = 12
ly = 2.9
lz = 3.5
myfile = open('ohslit_out.gro','w')
myfile.write('Undefined Name\n')
myfile.write('%5d\n' %(00000))

ly_real = gen_circle([0, 0], lx, ly, lz, myfile, 2, symm_choice)
gen_circle([0.341, 0.341], lx-0.682, ly, lz-0.682, myfile, 1, symm_choice)
gen_circle([0.682, 0.682], lx-0.682*2, ly, lz-0.682*2, myfile, 1, symm_choice)

myfile.write("%10.5f%10.5f%10.5f\n" %(lx ,ly_real ,lz))
myfile.close()

