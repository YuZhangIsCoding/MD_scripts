#!/Users/yuzhang/anaconda/bin/python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', dest = 'surfc', default = 0, type = float, help = 'Surface charge density')
parser.add_argument('-i', '--input', dest = 'filename', default = 'graphene_funct.itp', help = 'itp file for the neutral graphene sheet')
#parser.add_argument('-u', '--unit', dest = 'unit', default = [24, 36, 24], narg = '*', help = 'deminsion of the oxidized graphene')
parser.add_argument('-b', '--bond', dest = 'bond', default = 0.142, help = 'bondlength of C-C')
args = parser.parse_args()
unit = [24, 36, 24]

delta = args.surfc*unit[0]*unit[1]*3/4*3**0.5/(unit[0]*unit[1]+2*unit[2])*6.242*args.bond**2
print delta

file_in = open(args.filename)

pos = open('POS.itp', 'w')
neg = open('NEG.itp', 'w')
flag = 0
for count, line in enumerate(file_in):
    if count == 2:
        pos.write('POS  3\n')
        neg.write('NEG  3\n')
        continue
    elif 'atoms' in line:
        flag = 1
        pos.write('[ atoms ]\n')
        neg.write('[ atoms ]\n')
        continue
    elif 'bonds' in line:
        flag = 0
    if flag == 1 and line[0] != ';':
        temp = line.split()
        temp[-2] = float(temp[-2])+delta
        pos.write("%5s%5s%5s%5s%5s%5s%12.7f%12s\n" %tuple(temp))
        temp[-2] -= 2*delta
        neg.write("%5s%5s%5s%5s%5s%5s%12.7f%12s\n" %tuple(temp))
    else:
        pos.write(line)
        neg.write(line)
neg.close()
pos.close()
file_in.close()
