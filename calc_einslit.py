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

e_all = [[] for row in range(4)]

chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
#        for i in range(4):
#            e_all[i].append(0)
        e_all[0].append(32.73)
        e_all[1].append(0)
        e_all[2].append(-32.73)
        e_all[3].append(0)
        for xyz_ind, xyz in enumerate(frame.xyz[0, :, direction]):
            if xyz >= 0 and xyz <= 12:
                e_all[0][-1] += chargeind[xyz_ind]
            elif xyz > 12 and xyz < 19:
                e_all[1][-1] += chargeind[xyz_ind]
            elif xyz >= 19 and xyz<= 31:
                e_all[2][-1] += chargeind[xyz_ind]
            else:
                e_all[3][-1] += chargeind[xyz_ind]
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if chunk_index == 0:
#        break

########## Outputs ##########
for i in range(4):
    e_all[i] = np.array(e_all[i])
    if i == 0:
        data = np.column_stack((np.arange(len(e_all[0]))*comm.args.dt/1000.0, e_all[i]))
    else:
        data = np.column_stack((data, e_all[i]))
np.savetxt('E_all.txt', data, fmt = ['%12.4f' for i in range(len(data[0]))])

fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (16, 16), dpi = 1000)
region = ['slit 1', 'bulk 1', 'slit 2', 'bulk 2']
for i in range(2):
    for j in range(2):
        c = 2*i+j
        axes[i, j].plot(data[:, 0], data[:, c+1], 'k-', linewidth = lwidth, label = 'Charge - '+region[c])
axes[1, 0].set_xlabel('Time (ns)', fontsize = fsize)
axes[1, 1].set_xlabel('Time (ns)', fontsize = fsize)
axes[0, 0].set_ylabel('Charge in the slit (e)', fontsize = fsize)
axes[1, 0].set_ylabel('Charge in the slit (e)', fontsize = fsize)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
    axes.legend(loc = 'best', fontsize = 16)
    axes.grid()
plt.tight_layout()
plt.savefig('E_all.pdf')
