#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_template.py
# Description:  This is a python script to calculate the angle distribution of a molecule near a planar surface
#               Usuage example: calc_angle_planar.py -i test.xtc -bs 2 -g 2 3 --aname C3 C5 C6 C1T C5T -b 0 2.42 2.42 3.687 2.687 10 --suffix reg
# Date: 02-22-2016 Created

import time
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
angles = np.arange(0, 180, bin_size)
vec_ind = comm.select_vec(res_targ, topology)
n_angle = [np.zeros((len(angles), len(res_targ))) for row in range(len(bound)/2)]
pre_time = time.time()
if comm.args.suffix == None:
    outname = 'Angle'
else:
    outname = 'Angle_'+comm.args.suffix
chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        angle = comm.calc_angle(res_targ, vec_ind, frame, topology, para, direction, bound)
        for i in range(len(bound)/2):
            for j in range(len(res_targ)):
                n_temp , _ = np.histogram(angle[i][j], bins = np.append(angles, max(angles)+bin_size))
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
    plt.xlabel('$\mathregular{\\theta (^{\circ})}$')
    plt.ylabel('Probability density')
    plt.savefig('%s_%s.pdf' %(outname, i))
