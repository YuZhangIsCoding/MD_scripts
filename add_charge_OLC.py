#!/Users/yuzhang/anaconda/envs/py3/bin/python
import calc_common as comm
import numpy as np


print('Reading file %s ...' %comm.args.filename)

r = 1.018
num = 500
def charge_density():
    print('Adding surface charge density around %f C/m^2' %comm.args.sc)
    extra_charge = comm.args.sc*6.242*4*np.pi*r**2
    if extra_charge > 0:
        if extra_charge < 1:
            extra_charge = 1
        else:
            extra_charge = int(extra_charge)
    elif extra_charge < 0:
        if extra_charge > -1:
            extra_charge = -1
        else:
            extra_charge = int(extra_charge)
    return extra_charge

def charge_e():
    return comm.args.sc

extra_charge = charge_e()
delta = float(extra_charge)/num
print('Each surface atom has partial charge of %.8f e' %delta)
real_sc = float(extra_charge)/6.242/4/np.pi/r**2
print('Real surface charge density: %.8f C/m^2' %real_sc)
print('Net charge to remove from electrolytes: %d e' %extra_charge)

myfile = open(comm.args.filename, 'r')
outfile = open('OCO.itp', 'w')
print('Writing to OCO.itp ...')
outfile.write('; surface charge density: %.8f C/m^2\n' %real_sc)
flag = 0
for i, line in enumerate(myfile):
    if line[0] == ';':
        continue
    elif '[' in line and ']' in line and 'atoms' in line:
        flag = 1
        outfile.write(line)
    elif '[' in line and ']' in line:
        flag = -flag
        outfile.write(line)
    elif flag == 1:
        temp = line.split()
        outfile.write('%5s%5s%5s%5s%5s%5s%12.8f%12s\n' %(tuple(temp[:6])+(delta, temp[-1])))
    else:
        outfile.write(line)

outfile.close()
myfile.close()
