#!/Users/yuzhang/anaconda/bin/python
# Filename: gen_samp_dist.py
# Discription: this is a python script to generate different configurations for umbrella samling.
# Date: 09-21-2015 created

import os
import numpy as np

sep = [0.05, 0.1]
file_in = open('pullx.xvg', 'r')
file_out = open('sample_dist.txt','w')
p_start = None
p_sep = [1, 5]
zone = 0
count = 0
for line in file_in:
    if line[0] != '#' and line[0] != '@':
        temp = [float(i) for i in line.split()]
        if p_start == None:
            p_start = temp[1]+np.mean(temp[2:])
            p_jdg = p_start
        if temp[0]*5%1 == 0:
            p_now = temp[1]+np.mean(temp[2:])
            if p_now-p_start >= p_sep[0]:
                zone = 1
            if p_now-p_start > p_sep[1]:
                break
            if p_now >= p_jdg:
                count += 1
                if p_now < p_jdg+sep[zone]:
                    print temp[0], p_now
                    os.system("echo 0| gmx_mpi trjconv -b %f -e %f -o whatever/conf%d.gro" %(temp[0], temp[0], count))
                    file_out.write('%5d%12s%10.5f%10.5f\n' %(count, temp[0], p_now, p_jdg))
                p_jdg += sep[zone]


file_in.close()
file_out.close()
