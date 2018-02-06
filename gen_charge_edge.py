#!/Users/yuzhang/anaconda/bin/python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', dest = 'surfc', default = 0, type = float, help = 'Surface charge density')
parser.add_argument('-i', '--input', dest = 'filename', default = 'graphene_arm.gro', help = 'gro file for the neutral graphene sheet')
#parser.add_argument('-u', '--unit', dest = 'unit', default = [24, 36, 24], narg = '*', help = 'deminsion of the oxidized graphene')
parser.add_argument('-b', '--bond', dest = 'bond', default = 0.142, help = 'bondlength of C-C')
args = parser.parse_args()

he = 0.1268

file_in = open(args.filename)
count = 0
#cidx = []
hidx = []
for idx, line in enumerate(file_in):
    if 'HG' in line:
        hidx.append(idx-1)
#    elif 'CG' in line:
#        cidx.append(idx-1)
box = [float(i) for i in line.split()]
delta = args.surfc*box[0]*box[1]/len(hidx)*2*6.2415
file_in.close()
print delta

pos = open('GPO.itp', 'w')
neg = open('GNE.itp', 'w')
flag = 0
pos.write('[ moleculetype ]\n; molname       nrexcl\nGPO             1\n[ atoms ]\n')
neg.write('[ moleculetype ]\n; molname       nrexcl\nGNE             1\n[ atoms ]\n')

#for i in range(idx-2):
#    if i in hidx:
#        pos.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(i+1, 'HG', 1, 'GPO', 'HG', i+1, he+delta, 1.008))
#    else:
#        pos.write('%5d%5s%5d%5s%5s%5d%12.7f%8.3f\n' %(i+1, 'CG', 1, 'GPO', 'CG', i+1, he+delta, 1.008))

file_in = open(args.filename)
for count, line in enumerate(file_in):
    if count == 2:
        pos.write('GPO  3\n')
        neg.write('GNE  3\n')
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
