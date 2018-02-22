#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: gen_functionalized.py
# Description: This is a python script that creates the configuration file 
# Date: 03-02-2016 Created

import numpy as np
import argparse

############ bonded parameters ##########
bondlen = 0.142
bcc = 0.151
bco = 0.141
boh = 0.096
acoh = 108.5
resname = 'GOH'

############ argparse ###########
def myinput(inote):
    temp = raw_input(inote)
    return temp.split()
parser = argparse.ArgumentParser(description = 'parameters for the functionalized carbon layer')
parser.add_argument('-d', '--dimention', dest = 'dim', nargs = '*', type = int, help = 'dimension of the carbon frame')
parser.add_argument('-n', '--nlayers', dest = 'n', type = float, default = 3, help = 'total number of graphene layers')
parser.add_argument('-r', '--random', dest = 'random', action = 'store_true', help = 'random configuration of functional groups')
parser.add_argument('-p', '--percent', dest = 'percent', type = float, help = 'percentage of functional groups')
parser.add_argument('-u', '--unitcell', dest = 'unit', nargs = '*', type = int, help = 'unitcell for regular distributed functionalized graphene')
parser.add_argument('-g', '--group', dest = 'group', choices = ['H', 'O', 'OH'], help = 'Select the functional group')
#args = parser.parse_args(['-d', '24', '40', '-r', '-g', 'OH', '-p', '5'])
args = parser.parse_args(['-d', '48', '80', '-u', '4', '4', '-g', 'OH', '-p', '5'])
if args.dim == None:
    temp = raw_input('Please input the dimension for the graphene layer:\n').split()
    dim = [int(i) for i in temp]
else:
    dim = args.dim
groupname = ['pristine', 'H', 'O', 'OH']
if args.group == None:
    group = int(raw_input('Please select the functional group:\n\
    0.  pristine\n\
    1.  H\n\
    2.  O\n\
    3.  OH\n'))
else:
    group = args.group
    group = groupname.index(group)

########### functions ##########
def gen_frame(dim, bondlen, grofile, df = []):
    '''gen_frame(dim) -> gframe
    dim:        a list containing the dimension for the frame
    bondlen:    bond length of C-C bond in graphene
    grofile:     gro file to write the coordinate file
    df:         the list of functional groups

    Generate general frame for a single layer of graphene.
    Return a list of coordinates according to the atom index.'''
    gframe = np.array([np.zeros(dim)]*3)
    if df != []:
        for coord in df:
            gframe[2][tuple(coord[:2])] += np.sqrt(bcc**2-bondlen**2)*coord[2]
    count = 0
    for row in range(dim[0]):
        for col in range(dim[1]):
            count += 1
            gframe[0][row, col] = ((1-(-1)**(row+col))/4.0+1.5*row)*bondlen
            gframe[1][row, col] = (np.sqrt(3)/2*col)*bondlen
            if [row, col] in [item[:2] for item in df]:
                grofile.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" %(1, resname, 'COH', count, gframe[0][row, col], gframe[1][row,col], gframe[2][row, col]))
                itpfile.write("%5d%5s%5d%5s%5s%5d%12.7f%12.4f\n" %(count, 'COH', 1, 'G'+groupname[group], 'COH', count, 0.373, 12.011))
            else:
                grofile.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" %(1, resname, 'CG', count, gframe[0][row, col], gframe[1][row,col], gframe[2][row, col]))
                itpfile.write("%5d%5s%5d%5s%5s%5d%12.7f%12.4f\n" %(count, 'CG', 1, 'G'+groupname[group], 'CG', count, 0, 12.011))
    return gframe, count

def append_OH(gframe, df):
    '''add hydroxyl groups'''
    fframe = [np.zeros((len(df), 3)) for row in range(2)]
    for i, item in enumerate(df):
        for j in range(3):
            fframe[0][i, j] = gframe[j][tuple(item[:2])]
            fframe[1][i, j] = gframe[j][tuple(item[:2])]
        fframe[0][i, 2] += item[2]*bco
        fframe[1][i, 0] += boh*np.cos(acoh/180*np.pi)
        fframe[1][i, 2] = fframe[0][i, 2] + item[2]*boh*np.sin(acoh/180*np.pi)
    return fframe

############ main ############
grofile = open('graphene_funct.gro', 'w')
grofile.write('functionalized graphene layer\n')
itpfile = open('graphene_funct.itp', 'w')
itpfile.write('[ moleculetype ]\n; molname   nrexcl\nG%s    3\n[ atoms ]\n'  %groupname[group])
itpfile.write('; id type    res residu  at  cg  charge   mass\n')

df = []
if args.random == True:
    nf = int(args.percent/100*dim[0]*dim[1])
    for i in range(nf):
        while True:
            temp = [np.random.random_integers(0, j-1) for j in dim]
            if temp not in [i[:2] for i in df]:
                df.append(temp+[(-1)**np.random.random_integers(0, 1)])
                break
else:
    temp = 0
    for row in range(dim[0]):
        for col in range(dim[1]):
            if row%args.unit[0] == 1 and col%args.unit[1] == 1:
                temp += 1
                #df.append([row, col, (-1)**np.random.random_integers(0, 1)])
                df.append([row, col, (-1)**(temp)])
grofile.write('%5d\n' %(dim[0]*dim[1]+2*len(df)))
gframe, count = gen_frame(dim, bondlen, grofile, df)
if group == 3:
    fframe = append_OH(gframe, df)
    atomname = ['OG', 'HG']
    atommass = [15.999, 1.008]
    atomcharge = [-0.7890, 0.416]
    for i, subgroup in enumerate(fframe):
        for item in subgroup:
            count += 1
            grofile.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" %(1, resname, atomname[i], count, item[0], item[1], item[2]))
            itpfile.write("%5d%5s%5d%5s%5s%5d%12.7f%12.4f\n" %(count,atomname[i], 1, 'G'+groupname[group], atomname[i], count, atomcharge[i], atommass[i]))
grofile.write("%10.5f%10.5f%10.5f\n" %(dim[0]*3/2*bondlen, dim[1]*3**0.5/2*bondlen, 5))
#add_bond()
itpfile.write('[ bonds ]\n')
indx = []
oidx = []
hidx = []
up = []
down = []
side = []
for i, item in enumerate(df):
    indx.append(item[0]*dim[1]+item[1]+1)
    oidx.append(dim[0]*dim[1]+i+1)
    hidx.append(oidx[-1]+len(df))
    if (item[0]+item[1])%2 == 0:
        if item[0] == 0:
            side.append((dim[0]-1)*dim[1]+item[1]+1)
        else:
            side.append(indx[-1]-dim[1])
    else:
        if item[0] == dim[0]-1:
            side.append(item[1]+1)
        else:
            side.append(indx[-1]+dim[1])
    if item[1]+1 == dim[1]:
        up.append(indx[-1]+1-dim[1])
    else:
        up.append(indx[-1]+1)
    if item[1] == 0:
        down.append(indx[-1]-1+dim[1])
    else:
        down.append(indx[-1]-1)
    itpfile.write('%5d%5d%5d\n' %(indx[-1], side[-1], 1))
    itpfile.write('%5d%5d%5d\n' %(indx[-1], up[-1], 1))
    itpfile.write('%5d%5d%5d\n' %(indx[-1], down[-1], 1))
    itpfile.write('%5d%5d%5d\n' %(indx[-1], oidx[-1], 1))
    itpfile.write('%5d%5d%5d\n' %(oidx[-1], hidx[-1], 1))
    #print 'indx, side, up, down, oidx, hidx', indx[-1], side[-1], up[-1], down[-1], oidx[-1], hidx[-1]
#angles
itpfile.write('[ angles ]\n')
for i in range(len(indx)):
    itpfile.write('%5d%5d%5d%5d\n' %(side[i], indx[i], oidx[i], 1))
    itpfile.write('%5d%5d%5d%5d\n' %(up[i], indx[i], oidx[i], 1))
    itpfile.write('%5d%5d%5d%5d\n' %(down[i], indx[i], oidx[i], 1))
    itpfile.write('%5d%5d%5d%5d\n' %(indx[i], oidx[i], hidx[i], 1))
# dihedrals
itpfile.write('[ dihedrals ]\n')
dih = []
for i in range(len(indx)):
    for j in [side[i], up[i], down[i]]:
        if j in indx and min(indx[i], j) not in dih:
            print indx[i], j
            dih.append(min(indx[i], j))
grofile.close()
itpfile.close()
