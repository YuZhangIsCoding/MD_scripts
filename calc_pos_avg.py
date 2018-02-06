#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_pos_avg.py
# Description:  This is a python script that calculate the average position of the charge plane
#               for the cases when not all surface charges are located at the same plane.
# Date: 05-30-2017 Created

import numpy as np
import calc_common as comm
myfile = open(comm.args.filename, 'r')
mylist = []
for count, line in enumerate(myfile):
    if count > 1:
        temp = line.split()
        if temp[0] == '1GPO':
            if 'CG' not in temp[1]:
                mylist.append(float(temp[5]))
        else:
            break
myfile.close()
print np.mean(mylist)
