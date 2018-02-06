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
n_angle = [np.zeros((len(angles), len(res_targ))) for row in range(len(bound)/2)]
pre_time = time.time()
if comm.args.suffix == None:
    outname = 'cosangle'
else:
    outname = 'cosangle_'+comm.args.suffix
chunk_size = 100

for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        cosangle = comm.calc_cosangle(res_targ, vec_ind, frame, topology, para, direction, bound)
        for i in range(len(bound)/2):
            for j in range(len(res_targ)):
                n_temp , _ = np.histogram(cosangle[i][j], bins = np.append(angles, max(angles)+bin_size))
                n_angle[i][:, j] += n_temp
        if sub_ind %10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
    if comm.args.test == True:
        if chunk_index == 0:
            break
frames_read = chunk_index*chunk_size+sub_ind+1
print frames_read
for i, item in enumerate(n_angle):
    n_angle[i] /= sum(n_angle[i])*bin_size
    np.savetxt('%s_%s.txt' %(outname, i), np.column_stack((angles, n_angle[i])), fmt = ['%12.4f' for row in range(len(res_targ)+1)])
    fig = plt.figure()
    plt.plot(angles, n_angle[i])
    plt.xlabel('cos(${\\theta}$)')
    plt.ylabel('Probability density')
    plt.savefig('%s_%s.pdf' %(outname, i))
