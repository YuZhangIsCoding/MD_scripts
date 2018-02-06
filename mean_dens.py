#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: mean_dens.py
# Description:  This is a python script to get the mean number density profile at central (specified) region
#               If number density (.txt) file is input, the averaged number density will be output
#               If mass density (.xvg) file is input, please also indicate the molar mass for each group
# Date:     02-15-2016  Created 
import numpy as np
import sys, os
import calc_common as comm

file_name = comm.args.filename
if not os.path.isfile(file_name):
    sys.exit('Exit: no file named %s found!' %file_name)

if comm.args.bound != None:
    bound = [float(i) for i in comm.args.bound]
else:
    bound = [float(i) for i in input('Please input the boundary of interest').split()]
if len(bound)%2 != 0:
    sys.exit('Exit: boundary not recognized')

if file_name.endswith('.txt'):
    data_in = np.loadtxt(file_name)
elif file_name.endswith('xvg'):
    data_in = comm.load_xvg(file_name)

data = [[] for i in range(int(len(bound)/2))]
for temp in data_in:
    for i in range(int(len(bound)/2)):
        if float(temp[0]) > bound[2*i] and float(temp[0]) < bound[2*i+1]:
            data[i].append([float(i) for i in temp])
data = np.array(data)
for i in range(len(data)):
    data[i] = np.array(data[i])
if file_name.endswith('.xvg'):
    if comm.args.mmass != None:
        mmass = comm.args.mmass
    else:
        mmass = [float(i) for i in input('Please input the molar mass for %d group(s):' %(len(data[0][0])-1)).split()]
    if len(mmass) != len(data[0][0])-1:
        sys.exit('Exit: please input the exact number of molmasses for the groups selected!')
    for i in range(int(len(bound)/2)):
        print('Mean number density for region %d:' %(i+1), [np.mean(data[i][:, j])/mmass[j-1] for j in range(1,len(data[i][0]))])
else:
    for i in range(int(len(bound)/2)):
        print('Mean number density for region %d:' %(i+1), [np.mean(data[i][:, j]) for j in range(1,len(data[i][0]))])
