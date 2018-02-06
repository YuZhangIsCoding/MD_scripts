import numpy as np
import pdb

numd = np.loadtxt('NumberDensity.txt',dtype=[('distance',float),('emim',float),('bf4',float),('tf2n',float)])
dz = numd[1][0]-numd[0][0]
mark_b = 0
sum_z = 0
sum_n = 0
for item in numd:
    sum_z += item[0]*(item[2]+item[3])*dz
    sum_n += (item[2]+item[3])*dz
    if mark_b == 0 and item[2] > 1:
        mark_b = 1
    if mark_b == 1 and max(item[2], item[3])<0.01:
        break
print sum_z/sum_n-0.682


mark_b = 0
sum_z = 0
sum_n = 0
for i in range(len(numd)):
    item = numd[len(numd)-i-1]
    item[0] = numd[i][0]
    sum_z += item[0]*item[1]*dz
    sum_n += item[1]*dz
    if mark_b == 0 and item[1] > 1:
        mark_b = 1
    if mark_b == 1 and item[1] < 0.1:
        break
print sum_z/sum_n-0.682
