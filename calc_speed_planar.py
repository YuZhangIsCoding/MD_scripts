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


########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
direction = comm.direction
bin_size = comm.bin_size+1
bound = comm.load_bound(traj0, direction)
distances = np.linspace(bound[0], bound[1], bin_size)

print 'Using dt (ps)', comm.args.dt

chunk_size = 100
xyz_pre = comm.calc_xyz_com(res_targ, traj0, topology, para)
box = traj0.unitcell_lengths[0, 0]
x_st = 1 # interface that flux cross
speed = [[] for _ in res_targ]
#per = [[0 for i in range(len(j))] for j in res_targ]
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        for res_ind, xyz_res in enumerate(xyz_com):
            diff = xyz_res[:, 0]-xyz_pre[res_ind][:, 0]
            for i, item in enumerate(diff):
                if item > box/2:
                    diff[i] -= box
                elif item < -box/2:
                    diff[i] += box
            temp_1, _ = np.histogram(xyz_res[:, 2], weights = diff, bins = distances)
            temp_2, _ = np.histogram(xyz_res[:, 2], bins = distances)
            for i, item in enumerate(temp_2):
                if item != 0:
                    temp_1[i] /= item
            speed[res_ind].append(temp_1)
        xyz_pre = xyz_com
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
    if comm.args.test == True:
        if chunk_index == 0:
            break

fig, axes = plt.subplots(len(res_targ), sharex = True)
for res in range(len(speed)):
    speed[res] = np.array(speed[res])
    for j in range(len(speed[res][0])):
        axes[res].plot(np.arange(len(speed[res]))*comm.args.dt/1000, speed[res][:, j]/comm.args.dt)
    axes[res].set_ylabel('m/s')
plt.xlabel('ns')
plt.savefig('speed_t.pdf')

speed_avg = []
fig, axes = plt.subplots(len(res_targ), sharex = True)
for res in range(len(speed)):
    speed_avg.append(np.mean(speed[res], axis = 0)/comm.args.dt*1000)
    axes[res].plot(distances[:-1], speed_avg[-1])
    axes[res].set_ylabel('m/s')
plt.xlabel('nm')
plt.savefig('speed.pdf')
speed_avg = np.array(speed_avg).T
np.savetxt('speed_avg.txt', np.column_stack((distances[:-1], speed_avg)), fmt = ['%14.6f' for _ in range(len(speed_avg[0])+1)])
