#!/Users/yuzhang/anaconda/envs/py3/bin/python
# This is a python script that generates the coordinate file for graphene sheets with holes
# The idea is to first generate contact graphene sheets, and remove atoms found within the
# hole radius. Build small unit first, and expand the pattern to a larger plane.
# Note that we don't consider periodic condition here.

import argparse
import numpy as np

parser = argparse.ArgumentParser(description = 'Specify the size of the electrode')
parser.add_argument('-f', '--filename', dest = 'filename', help = 'Specify the filename of unitcell of graphene')
parser.add_argument('-r', '--radius', dest = 'radius', type = float, help = 'cutoff radius of the hole')
parser.add_argument('-c', '--center', dest = 'center', nargs = '*', type = float, help = 'center position of the hole')
parser.add_argument('-n', '--nlayer', type = int, default = 3, help = 'Specify the number of layers, default is 3')
parser.add_argument('-rep', '--repeat', dest = 'rep', nargs = 2, type = int, help = 'repeat rate in x and y dimension')
args = parser.parse_args('-rep 5 4 -c 0.5 0.5 -r 0.25 -f graphene_sheets.gro'.split())

if args.filename:
    filename = args.filename
else:
    filename = input('Please input the filename of unitcell:\n')

if args.center:
    center = np.array(args.center)
else:
    center = np.array([float(i) for i in input('Please input center coordinates:\n')])

if args.radius:
    radius = args.radius
else:
    radius = float(input('Please input the radius of the hole:\n'))

def load_gro(filename):
    myfile = open(filename, 'r')
    coord = []
    for i, line in enumerate(myfile):
        temp = line.split()
        if i > 1 and len(temp) >3:
            coord.append([line[:20], np.array([float(item) for item in temp[-3:]])])
    coord.append(np.array([float(item) for item in temp]))
    myfile.close()
    return coord

def check_range(coord, center, radius):
    n = len(center)
    return sum((coord[:n]-center)**2) < radius**2

coord = load_gro(filename)
i = 0
while i < len(coord)-1:
    if check_range(coord[i][1], center, radius):
        coord.pop(i)
    else:
        i += 1

myfile = open('graphene_holes.gro', 'w')
myfile.write('This is a graphene sheet with holes\n')
myfile.write('%5d\n' %((len(coord)-1)*args.rep[0]*args.rep[1]*args.nlayer))
for n in xrange(args.nlayer):
    for i in range(args.rep[0]):
        for j in range(args.rep[1]):
            base = [i*coord[-1][0], j*coord[-1][1]]
            for item in coord[:-1]:
                myfile.write('%s%8.3f%8.3f%8.3f\n' %(item[0], item[1][0]+base[0]+n%2*0.142, item[1][1]+base[1], item[1][2]+0.341*n))
myfile.write('%10.5f%10.5f%10.5f\n' %(coord[-1][0]*args.rep[0], coord[-1][1]*args.rep[1], coord[-1][2]))
myfile.close()
