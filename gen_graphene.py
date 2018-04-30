#!/Users/yuzhang/anaconda3/bin/python
# This is a python script that generates the coordinate file for graphene sheets.
# For MD simulations, graphene is usually modeled as an infinite plane. Due to 
# the periodic boundary conditions, the number of points in each edge must be
# even. If you want to build a free standing graphene without the boundary
# limits, please use --real flag.

import argparse

parser = argparse.ArgumentParser(description = 'Specify the size of the electrode')
parser.add_argument('-s', '--size', nargs = '*', type = float, help = 'Specify the approximate size of one graphene sheet')
parser.add_argument('-n', '--nlayer', type = int, default = 3, help = 'Specify the number of layers, default is 3')
parser.add_argument('--bd', type = float, default = 0.142, help = 'Specify the bondlength of C-C bond')
parser.add_argument('--real', action = 'store_false', help = 'won\'t use periodic boundary conditions')
parser.add_argument('-v', action = 'store_true', help = 'show details')
args = parser.parse_args()
bondlen = args.bd

############ functions #################
def gf_iter(l, h, dx = 0.5, dy = 0.5*3**0.5):
    '''generator that returns x and y coordinates'''
    for i in range(l):
        if i%2 == 0:
            x = 3*i*dx+dx
            carry = 1
        else:
            x = 3*i*dx
            carry = -1
        for j in range(h):
            carry = -carry
            x += dx*carry
            yield x, dy*j

########### main ###################3
if args.size:
    size = [float(i) for i in args.size]
else:
    size = [float(item) for item in input('Please input the x and y size of the graphene:').split()]

if args.real:
    l = int(round(size[0]/bondlen*2/6)*2)
    h = int(size[1]/bondlen/3**0.5)*2
else:
    l = int(round(size[0]/bondlen*2/3))
    h = int(round(2*size[1]/bondlen/3**0.5))

if args.v:
    print('Dimension size of the sheet:', l, h, args.nlayer)
    print('Real size of the graphene sheet: %.3f, %.3f' %(bondlen*l*3/2, bondlen*h*3**0.5/2))
    print('Writing to file graphene_sheets.gro ...')

grofile = open('graphene_sheets.gro','w')
grofile.write('%d-layer %dx%d graphene\n' %(args.nlayer, l, h))
grofile.write('%5d\n' %(l*h*args.nlayer))

count = 0
for n in range(args.nlayer):
    for (x, y) in gf_iter(l, h):
        count += 1
        grofile.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" %(count, 'GPH', 'CG', count, (x+n%2)*bondlen, y*bondlen, n*0.341))

grofile.write("%10.5f%10.5f%10.5f\n" %(l*3/2*bondlen, h/2*3**0.5*bondlen, args.nlayer*0.341))
grofile.close()
