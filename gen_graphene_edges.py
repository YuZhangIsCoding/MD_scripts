#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: gen_graphene_edges.py
# Description:  A python script that generates the coordinate file for graphene
#               edge planes, including armchair and zigzag. This script will 
#               generate 2 electrodes, cathode and anode respectively. These
#               two electrode are symmetric in the z direction. Along with a
#               gro file, the script will also produce two itp files, where 
#               the surface charge density is opposite from each other.

import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'User specified substrate size')
parser.add_argument('-s', '--size', dest = 'size', nargs = '*', help = 'define the size of the electrode')
parser.add_argument('-n', '--name', dest = 'name', help = 'output file name')
parser.add_argument('-e', '--charge', dest = 'charge', default = 0, type = float, help = 'surface charge density')
parser.add_argument('-c', '--choice', dest = 'choice', default = 'arm', choices = ['armchair', 'zigzag'], help = 'define the type of edge')
parser.add_argument('-v', action = 'store_true', help = 'show details')
args = parser.parse_args()

########## functions ##############
def load_armchair():
    dim = [int(round(size[0]/bondlen*2/6)*2), int(round(size[1]/gap)), int(size[2]/bondlen*2/3**0.5/2)*2+1]
    box = [bondlen*3/2*dim[0], gap*dim[1], 5]
    n_t = np.prod(dim)+2*dim[0]*dim[1]
    delta = args.charge*box[0]*box[1]/dim[0]/dim[1]*6.2415
    return dim, box, n_t, delta

def load_zigzag():
    dim = [int(round(size[0]/bondlen/3**0.5)*2), int(round(size[1]/gap)), int(round(size[2]/bondlen*2/6)*2)]
    box = [bondlen*3**0.5/2*dim[0], gap*dim[1], 5]
    n_t = np.prod(dim)+dim[0]*dim[1]
    delta = args.charge*box[0]*box[1]/dim[0]/dim[1]*2*6.2415
    return dim, box, n_t, delta

def write_armchair(grofile, posfile, negfile, dim, delta):
    global he, bondlen, hcbond, gap, dx, dz
    count = 0
    for n in range(dim[1]):
        for l in range(dim[0]):
            for h in range(dim[2]):
                count += 1
                x = bondlen*(0.5*((l+h+1)%2)+1.5*l+n%2)
                y = gap*n
                z = bondlen*(h*0.5*3**(0.5))
                grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'CG', count, x, y, z))
                if h == 0 or h == dim[2]-1:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, -he, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, -he, 12.011))
                else:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, 0, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, 0, 12.011))
            count += 1
            grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'HG', count, x-(-1)**l*dx, y, -dz))
            posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GPO', 'HG', count, he, 1.008))
            negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GNE', 'HG', count, he-delta, 1.008))
            count += 1
            grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'HG', count, x-(-1)**l*dx, y, z+dz))
            posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GPO', 'HG', count, he+delta, 1.008))
            negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GNE', 'HG', count, he, 1.008))

def write_zigzag(grofile, posfile, negfile, dim, delta):
    global he, bondlen, hcbond, gap, dx, dz
    count = 0
    for n in range(dim[1]):
        for h in range(dim[2]):
            for l in range(dim[0]):
                count += 1
                x = bondlen*(l*0.5*3**0.5)
                y = gap*n
                z = bondlen*(0.5*((l+h)%2)+1.5*h+n%2)
                grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'CG', count, x, y, z))
                if h == 0 and l%2 == 0:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, -he, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, -he, 12.011))
                    count += 1
                    grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'HG', count, x, y, n%2*bondlen-hcbond))
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GPO', 'HG', count, he, 1.008))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GNE', 'HG', count, he, 1.008))
                elif h == dim[2]-1 and l%2 == 0:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, -he, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, -he, 12.011))
                    count += 1
                    grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'HG', count, x, y, z+hcbond))
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GPO', 'HG', count, he+delta, 1.008))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GNE', 'HG', count, he-delta, 1.008))
                else:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, 0, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, 0, 12.011))

############## main ##################
# some constants
he = 0.1268
bondlen = 0.142
hcbond = 0.108
gap = 0.341
dx = hcbond/2
dz = hcbond/2*3**(0.5)

if args.size:
    size = [float(i) for i in args.size]
else:
    size = [float(i) for i in input('Please enter the size of the electroe:\n').split()]

if args.choice == 'armchair':
    dim, box, n_t, delta = load_armchair()
else:
    dim, box, n_t, delta = load_zigzag()

if args.name:
    if args.name.endswith('.gro'):
        groname = args.name
    else:
        groname = args.name+'.gro'
else:
    groname = args.choice+'.gro'

if args.v:
    print('generating %s edges to file %s ...' %(args.choice, groname))
    print('generating force field files to GPO.itp and GNE.itp ...')
    print('dim of the electrode:', dim)
    if args.choice == 'armchair':
        print('Real size of the electrode: %8.3f, %8.3f, %8.3f' %(bondlen*3/2*dim[0], gap*dim[1], bondlen*3**0.5/2*dim[2]))
    else:
        print('Real size of the electrode: %8.3f, %8.3f, %8.3f' %(bondlen*3**0.5/2*dim[0], gap*dim[1], bondlen*3/2*dim[2]))
    print('partial charge on surface H atoms:', delta)

grofile = open(groname, 'w')
posfile = open('GPO.itp', 'w')
negfile = open('GNE.itp', 'w')

grofile.write('Grapheen edge with H\n')
grofile.write('%5d\n' %n_t)
posfile.write('[ moleculetype ]\n; molname       nrexcl\nGPO             1\n[ atoms ]\n')
negfile.write('[ moleculetype ]\n; molname       nrexcl\nGNE             1\n[ atoms ]\n')

if args.choice == 'armchair':
    write_armchair(grofile, posfile, negfile, dim, delta)
else:
    write_zigzag(grofile, posfile, negfile, dim, delta)

grofile.write("%10.5f%10.5f%10.5f\n" %(box[0], box[1], box[2]))
grofile.close()
posfile.close()
negfile.close()
