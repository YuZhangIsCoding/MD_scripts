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
real_sc = float(extra_charge)/6.242/4/np.pi/r**2
print('Real surface charge density: %.8f C/m^2' %real_sc)
print('Net charge to remove from electrolytes: %d e' %extra_charge)

def get_delta(ind):
    delta = float(extra_charge)/num*(-1)**ind
#    print('Each surface atom has partial charge of %.8f e' %delta)
    return delta
typenames = ['OCP', 'OCN']
file_list = [open(item+'.itp', 'w') for item in typenames]

def write_line(file_list, line, write_charge = False, write_title = False):
    for ind, file_iter in enumerate(file_list):
        temp = line.replace('OCO', typenames[ind])
        if write_charge:
            temp = temp.split()
            file_iter.write('%5s%5s%5s%5s%5s%5s%12.8f%12s\n' %(tuple(temp[:6])+(get_delta(ind), temp[-1])))
        elif write_title:
            if ind == 0:
                file_iter.write(line)
            else:
                file_iter.write(line[:26]+'-'+line[26:])
        else:
            file_iter.write(temp)

myfile = open(comm.args.filename, 'r')
write_line(file_list, '; surface charge density: %.8f C/m^2\n' %real_sc, write_title = True)

flag = 0
for i, line in enumerate(myfile):
    if line[0] == ';':
        continue
    elif '[' in line and ']' in line and 'atoms' in line:
        flag = 1
        write_line(file_list, line)
    elif '[' in line and ']' in line:
        flag = -flag
        write_line(file_list, line)
    elif flag == 1:
        write_line(file_list, line, True)
    else:
        write_line(file_list, line)

myfile.close()
for item in file_list:
    item.close()
