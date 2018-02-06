#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_template.py
# Description: This is a python script to calculate the angle distribution of a molecule near a planar surface
# Date: 02-22-2016 Created

import sys, os, pdb, time, argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

time_start = time.time()

########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
direction = comm.direction
bound = comm.load_bound(traj0, direction)
bin_size = comm.bin_size
angles = np.arange(-1, 1, bin_size)
vec_ind = comm.select_vec(res_targ, topology)
#n_angle = [np.zeros((len(angles), len(res_targ))) for row in range(len(bound)/2)]
pre_time = time.time()
distances = np.arange(bound[0], bound[1], bin_size)
n_angle = np.zeros((len(distances), len(res_targ))).T
n_count = np.zeros((len(distances), len(res_targ))).T
if comm.args.suffix == None:
    outname = 'PAngle'
else:
    outname = 'pAngle_'+comm.args.suffix
chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_dir = comm.calc_xyz_direction(res_targ, frame, topology, para, direction)
        angle = comm.calc_pangle(res_targ, vec_ind, frame, topology, para, direction, bound)
        n_temp, _ = np.histogram(xyz_dir[0], weights = angle[0][0], bins = np.append(distances, max(distances)+bin_size))
        n_angle += n_temp
        n_temp, _ = np.histogram(xyz_dir[0], bins = np.append(distances, max(distances)+bin_size))
        n_count += n_temp
        #for i in range(len(bound)/2):
        #    for j in range(len(res_targ)):
        #        n_temp , _ = np.histogram(angle[i][j], bins = np.append(angles, max(angles)+bin_size))
        #        n_angle[i][:, j] += n_temp
        if sub_ind %10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if chunk_index == 0:
#        break
frames_read = chunk_index*chunk_size+sub_ind+1
print frames_read
for i in range(len(n_count[0])):
    if n_count[0][i] != 0:
        n_angle[0][i] /= n_count[0][i]

for i, item in enumerate(n_angle):
    np.savetxt('%s_%s.txt' %(outname, i), np.column_stack((distances, n_angle[i])), fmt = ['%12.4f' for row in range(len(res_targ)+1)])
    fig = plt.figure()
    plt.plot(distances, n_angle[i])
    plt.xlabel('Distance (nm)')
    plt.ylabel('P2')
    plt.savefig('%s_%s.pdf' %(outname, i))
