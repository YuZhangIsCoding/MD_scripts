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

last_frame = md.load_frame('last_frame.xtc', 0, top = 'topol.pdb')
xyz_com = comm.calc_xyz_com(res_targ, last_frame, topology, para)

res_sel = [[], []]
x_t = [[], []]
if comm.args.rev == False:
    rev = 0
else:
    rev = 1

bound_sel = [0, 12, 19, 31] # boundaries for slit 1 and slit 2
for i in range(2):
    jdg = (i+rev)%2
    for j, xyz in enumerate(xyz_com[i][:, direction]):
        if xyz >= bound_sel[2*jdg] and xyz <= bound_sel[2*jdg+1]:
            res_sel[jdg].append(res_targ[i][j])
            x_t[jdg].append([])

chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_com = comm.calc_xyz_com(res_sel, frame, topology, para)
        for i in range(2):
            for j in range(len(res_sel[i])):
                x_t[i][j].append(xyz_com[i][j, direction])
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if chunk_index == 0:
#        break
x_t = np.array(x_t)
########## Outputs ##########
data1 = np.column_stack((np.arange(len(x_t[0][0]))*comm.args.dt/1000.0, x_t[0].T))
data2 = np.column_stack((np.arange(len(x_t[1][0]))*comm.args.dt/1000.0, x_t[1].T))
np.savetxt('X_t_slit1.txt', data1, fmt = ['%12.4f' for i in range(len(data1[0]))])
np.savetxt('X_t_slit2.txt', data2, fmt = ['%12.4f' for i in range(len(data2[0]))])
data = [data1, data2]

fig, axes = plt.subplots(2, figsize = (16, 16), dpi = 1000, sharex = True)
for i in range(2):
    for j in range(len(res_sel[i])):
        axes[i].plot(data[i][:, 0], data[i][:, j+1], linewidth = lwidth, label = 'Mol '+str(j+1))
    for j in bound_sel:
        axes[i].plot([0, data[i][-1, 0]], [j, j], '--', color = 'gray', linewidth = lwidth)
plt.xlabel('Time (ns)', fontsize = fsize)
for axes in fig.axes:
    axes.set_ylim([0, 38])
    axes.set_ylabel('X (nm)', fontsize = fsize)
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    axes.legend(loc = 'best', fontsize = 16, ncol = 6)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.grid()
plt.tight_layout()
plt.savefig('X_t_slit_all.pdf')
