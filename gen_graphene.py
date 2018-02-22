#!/Users/yuzhang/anaconda/envs/py3/bin/python
# This is a python script that generates the coordinate file for graphene sheets

import argparse

parser = argparse.ArgumentParser(description = 'Specify the size of the electrode')
parser.add_argument('-s', '--size', nargs = '*', help = 'Specify the approximate size of one graphene sheet')
parser.add_argument('-n', '--nlayer', type = int, default = 3, help = 'Specify the number of layers, default is 3')
parser.add_argument('--bd', type = float, default = 0.142, help = 'Specify the bondlength of C-C bond')
args = parser.parse_args()
bondlen = args.bd

size = [float(i) for i in args.size]
l = int(round(size[0]/bondlen*2/6)*2)
h = int(size[1]/bondlen/3**0.5)*2
n_layer = args.nlayer

print('Dimension size of the sheet:', l, h, n_layer)
print('Real size of the graphene sheet:', bondlen*l*3/2, bondlen*h*3**0.5/2)

wholepair = [[0]*4 for row in range(l*h)]

wholepair[0][0]=0
wholepair[0][1]=0

# general frame

for i in range(0,l):
    if (i%2 == 0):
        wholepair[h*i][0]=1.5*i
        wholepair[h*i+1][0]=wholepair[h*i][0]+0.5
    else:
        wholepair[h*i][0]=wholepair[h*(i-1)][0]+2
        wholepair[h*i+1][0]=wholepair[h*i][0]-0.5
    wholepair[h*i][1]=0
    wholepair[h*i+1][1]=0.5*3**0.5
    for j in range(2,h):
        wholepair[h*i+j][0]=wholepair[h*i+j-2][0]
        wholepair[h*i+j][1]=wholepair[h*i+j-2][1]+3**0.5

# pure graphene
grofile = open('graphene_sheets.gro','w')
grofile.write('%d-layer %dx%d graphene\n' %(n_layer,l,h))
atoms=l*h*n_layer
grofile.write('%5d\n' %(atoms))

for n in range(n_layer):
    for i in range(0,l*h):
        wholepair[i][3]='CG'
        grofile.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n" %(n*l*h+i+1,'GPH',wholepair[i][3], n*l*h+i+1,(wholepair[i][0]+n%2)*bondlen, wholepair[i][1]*bondlen, n*0.341))

x=l*3/2*bondlen
y=h/2*3**0.5*bondlen
z=(x+y)/2
grofile.write("%10.5f%10.5f%10.5f\n" %(x, y, z))
