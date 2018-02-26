#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: gen_graphene_edges.py
# Description:  A python script that generates the coordinate file for graphene
#               edge planes, including armchair and zigzag. This script will 
#               generate 2 electrodes, cathode and anode respectively. These
#               two electrode are symmetric in the z direction. Along with a
#               gro file, the script will also produce two itp files, where 
#               the surface charge density is opposite from each other.
#               Zigzag edge:            Armchair edge:
#                 H  H  H               H    HH    H
#                 |  |  |                \__/  \__/ 
#                / \/ \/ \               /  \__/  \

import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'User specified substrate size')
parser.add_argument('-s', '--size', dest = 'size', nargs = '*', type = float, help = 'define the size of the electrode')
parser.add_argument('-n', '--name', dest = 'name', help = 'output file name')
parser.add_argument('-e', '--charge', dest = 'charge', default = 0, type = float, help = 'surface charge density')
parser.add_argument('-c', '--choice', dest = 'choice', default = 'arm', choices = ['armchair', 'zigzag'], \
        help = 'define the type of edge')
parser.add_argument('-v', action = 'store_true', help = 'show details')
args = parser.parse_args()

########## functions ##############
class consts(object):
    '''Class that stores constants'''
    def __init__(self, args = []):
        self.he = 0.1268
        self.bondlen = 0.142
        self.hcbond = 0.108
        self.gap = 0.341
        self.dx = self.hcbond/2
        self.dz = self.hcbond/2*3**(0.5)
        if args:
            self.dim, self.box, self.n_t, self.delta, self.groname = self.load_consts(args)
            if args.v:
                print('generating %s edges to file %s ...' %(args.choice, self.groname))
                print('generating force field files to GPO.itp and GNE.itp ...')
                print('dim of the electrode:', self.dim)
                if args.choice == 'armchair':
                    print('Real size of the electrode: %.3f, %.3f, %.3f' \
                            %(self.bondlen*3/2*self.dim[0], self.gap*self.dim[1], self.bondlen*3**0.5/2*self.dim[2]))
                else:
                    print('Real size of the electrode: %.3f, %.3f, %.3f' 
                            %(self.bondlen*3**0.5/2*self.dim[0], self.gap*self.dim[1], self.bondlen*3/2*self.dim[2]))
                print('partial charge on surface H atoms:', self.delta)
    
    def load_consts(self, args):
        if args.size:
            size = args.size
        else:
            size = [float(i) for i in input('Please enter the size of the electroe:\n').split()]
    
        if args.name:
            if args.name.endswith('.gro'):
                groname = args.name
            else:
                groname = args.name+'.gro'
        else:
            groname = args.choice+'.gro'

        if args.choice == 'armchair':
            dim = [int(round(size[0]/self.bondlen*2/6)*2), int(round(size[1]/self.gap)), int(size[2]/self.bondlen*2/3**0.5/2)*2+1]
            box = [self.bondlen*3/2*dim[0], self.gap*dim[1], 5]
            n_t = np.prod(dim)+2*dim[0]*dim[1]
            delta = args.charge*box[0]*box[1]/dim[0]/dim[1]*6.2415
        else:
            dim = [int(round(size[0]/self.bondlen/3**0.5)*2), int(round(size[1]/self.gap)), int(round(size[2]/self.bondlen*2/6)*2)]
            box = [self.bondlen*3**0.5/2*dim[0], self.gap*dim[1], 5]
            n_t = np.prod(dim)+dim[0]*dim[1]
            delta = args.charge*box[0]*box[1]/dim[0]/dim[1]*2*6.2415
        return dim, box, n_t, delta, groname

def write_armchair(grofile, posfile, negfile, c):
    count = 0
    for n in range(c.dim[1]):
        for l in range(c.dim[0]):
            for h in range(c.dim[2]):
                count += 1
                x = c.bondlen*(0.5*((l+h+1)%2)+1.5*l+n%2)
                y = c.gap*n
                z = c.bondlen*(h*0.5*3**(0.5))
                grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'CG', count, x, y, z))
                if h == 0 or h == c.dim[2]-1:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, -c.he, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, -c.he, 12.011))
                else:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, 0, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, 0, 12.011))
            count += 1
            grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'HG', count, x-(-1)**l*c.dx, y, -c.dz))
            posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GPO', 'HG', count, c.he, 1.008))
            negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GNE', 'HG', count, c.he-c.delta, 1.008))
            count += 1
            grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'HG', count, x-(-1)**l*c.dx, y, z+c.dz))
            posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GPO', 'HG', count, c.he+c.delta, 1.008))
            negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GNE', 'HG', count, c.he, 1.008))

def write_zigzag(grofile, posfile, negfile, c):
    count = 0
    for n in range(c.dim[1]):
        for h in range(c.dim[2]):
            for l in range(c.dim[0]):
                count += 1
                x = c.bondlen*(l*0.5*3**0.5)
                y = c.gap*n
                z = c.bondlen*(0.5*((l+h)%2)+1.5*h+n%2)
                grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'CG', count, x, y, z))
                if h == 0 and l%2 == 0:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, -c.he, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, -c.he, 12.011))
                    count += 1
                    grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'HG', count, x, y, n%2*c.bondlen-c.hcbond))
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GPO', 'HG', count, c.he, 1.008))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GNE', 'HG', count, c.he, 1.008))
                elif h == c.dim[2]-1 and l%2 == 0:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, -c.he, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, -c.he, 12.011))
                    count += 1
                    grofile.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n' %(1, 'GPH', 'HG', count, x, y, z+c.hcbond))
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GPO', 'HG', count, c.he+c.delta, 1.008))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'HG', 1, 'GNE', 'HG', count, c.he-c.delta, 1.008))
                else:
                    posfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GPO', 'CG', count, 0, 12.011))
                    negfile.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(count, 'CG', 1, 'GNE', 'CG', count, 0, 12.011))

############## main ##################
c = consts(args)
grofile = open(c.groname, 'w')
posfile = open('GPO.itp', 'w')
negfile = open('GNE.itp', 'w')

grofile.write('Grapheen edge with H\n')
grofile.write('%5d\n' %c.n_t)
posfile.write('[ moleculetype ]\n; molname       nrexcl\nGPO             1\n[ atoms ]\n')
negfile.write('[ moleculetype ]\n; molname       nrexcl\nGNE             1\n[ atoms ]\n')

if args.choice == 'armchair':
    write_armchair(grofile, posfile, negfile, c)
else:
    write_zigzag(grofile, posfile, negfile, c)

grofile.write("%10.5f%10.5f%10.5f\n" %(c.box[0], c.box[1], c.box[2]))
grofile.close()
posfile.close()
negfile.close()
