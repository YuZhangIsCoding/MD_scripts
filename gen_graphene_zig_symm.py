#!/Users/yuzhang/anaconda/bin/python
# This is a python script that generates the coordinate file for graphene sheets

import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'User specified substrate size')
parser.add_argument('-s', '--size', dest = 'size', nargs = '*', help = 'define the size of the electrode')
parser.add_argument('-n', '--name', dest = 'name', default = 'graphene_zig', help = 'output file name')
parser.add_argument('-e', '--charge', dest = 'charge', default = 0, type = float, help = 'surface charge density')
#parser.add_argument('-c', '--choice', dest = 'choice', choices = ['arm', 'zig'], help = 'define the type of edge')
args = parser.parse_args()
if args.size == None:
    size = [float(i) for i in raw_input('Please enter the size of the electroe:\n').split()]
else:
    size = [float(i) for i in args.size]
he = 0.1268
bondlen = 0.142
hcbond = 0.108
gap = 0.341
dim = [int(round(size[0]/bondlen/3**0.5)*2), int(round(size[1]/gap)), int(round(size[2]/bondlen*2/6)*2)]
box = [bondlen*3**0.5/2*dim[0], gap*dim[1], 5]
print 'dim of the electrode:', dim
print 'Real size of the electrode:', bondlen*3**0.5/2*dim[0], gap*dim[1], bondlen*3/2*dim[2]
n_t = np.prod(dim)+dim[0]*dim[1]
grofile = open(args.name+'.gro', 'w')
grofile.write('Grapheen edge with H\n')
grofile.write('%5d\n' %n_t)
count = 0
dx = hcbond/2
dz = hcbond/2*3**(0.5)
delta = args.charge*box[0]*box[1]/dim[0]/dim[1]*2*6.2415
print 'partial charge on surface H atoms:', delta

posfile = open('GPO.itp', 'w')
negfile = open('GNE.itp', 'w')
posfile.write('[ moleculetype ]\n; molname       nrexcl\nGPO             1\n[ atoms ]\n')
negfile.write('[ moleculetype ]\n; molname       nrexcl\nGNE             1\n[ atoms ]\n')
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
grofile.write("%10.5f%10.5f%10.5f\n" %(box[0], box[1], box[2]))
grofile.close()
posfile.close()
negfile.close()
