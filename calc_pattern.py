#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_pattern.py
# Description:  This is a python script to get the heat map of the molecules in center layer.
#               This script uses 2d histogram to obtain a heat map of all atoms within a slice.
# Dates:    06-06-2017 Created

import mdtraj as md
import numpy as np
import calc_common as comm
import matplotlib.pyplot as plt

########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
#res_targ, res_name = comm.select_groups(traj0)
#massind, chargeind = comm.load_atomic(res_targ, topology, para)
direction = comm.direction
bin_size = comm.bin_size
bound = comm.load_bound(traj0, direction)

xedges = np.arange(0, traj0.unitcell_lengths[0, 0]+bin_size, bin_size)
yedges = np.arange(0, traj0.unitcell_lengths[0, 1]+bin_size, bin_size)
H = np.zeros((len(xedges)-1, len(yedges)-1))
chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xy = []
        for item in frame.xyz[0]:
            if item[2] <= bound[1] and item[2] >= bound[0]:
                xy.append(item[:2])
        xy = np.array(xy)
        temp, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins = (xedges, yedges))
        H += temp
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
    if comm.args.test == True and chunk_index == 0:
        break
frames_read = chunk_index*chunk_size+sub_ind+1

###### output ######
np.savetxt('Pattern_all.txt', H.T[::-1], fmt = ['%8d' for i in H])
print frames_read
fig = plt.figure(figsize = (16, 12), dpi = 300)
fsize = 28
plt.imshow(H.T[::-1], cmap = 'jet', interpolation='nearest', alpha = 1, extent = [0, 5.12, 0, 4.92])
plt.xlabel('x (nm)', fontsize = fsize)
plt.ylabel('y (nm)', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10, direction='out')
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()
plt.savefig('Pattern_all.pdf')
