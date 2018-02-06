#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_density_planar_iter.py
# Description: This is a python script to get the density profiles of a channel simulation
# Dates:    10-08-2015 Created
#           02-05-2015 Modified to split the trajectory by chunks to save memory

import sys, os, pdb, time
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

lwidth = 4
fsize = 28
########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
direction = comm.direction
bin_size = comm.bin_size
bound = comm.load_bound(traj0, direction)

n_all = [[[] for row in range(len(res_targ))] for row in range(4)]

chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        for j, item in enumerate(xyz_com):
            for i in range(4):
                n_all[i][j].append(0)
            for xyz in item[:, direction]:
                if xyz >= 1 and xyz <= 13:
                    n_all[0][j][-1] += 1
                elif xyz > 13 and xyz < 20:
                    n_all[1][j][-1] += 1
                elif xyz >= 20 and xyz<= 32:
                    n_all[2][j][-1] += 1
                else:
                    n_all[3][j][-1] += 1
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if chunk_index == 0:
#        break

########## Outputs ##########
for i in range(4):
    n_all[i] = np.array(n_all[i]).T
    if i == 0:
        data = np.column_stack((np.arange(len(n_all[0]))*comm.args.dt/1000.0, n_all[i]))
    else:
        data = np.column_stack((data, n_all[i]))
np.savetxt('N_all.txt', data, fmt = ['%12.4f' for i in range(len(data[0]))])

fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (16, 16), dpi = 1000)
region = ['slit 1', 'bulk 1', 'slit 2', 'bulk 2']
for i in range(2):
    for j in range(2):
        c = 2*i+j
        axes[i, j].plot(data[:, 0], data[:, 2*c+1], 'r-', linewidth = lwidth, label = 'OMI - '+region[c])
        axes[i, j].plot(data[:, 0], data[:, 2*c+2], 'b-', linewidth = lwidth, label = 'Tf2N - '+region[c])
axes[1, 0].set_xlabel('Time (ns)', fontsize = fsize)
axes[1, 1].set_xlabel('Time (ns)', fontsize = fsize)
axes[0, 0].set_ylabel('Number in the slit', fontsize = fsize)
axes[1, 0].set_ylabel('Number in the slit', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
    axes.legend(loc = 'best', fontsize = 16)
    axes.grid()
plt.tight_layout()
plt.savefig('N_all.pdf')
