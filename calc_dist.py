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
direction = comm.direction
bound = comm.load_bound(traj0, direction)
hlist = []
for atom in topology.atoms:
    if atom.name[0] == 'H':
        hlist.append(atom.index)
pair_list = []
for i in range(len(hlist)):
    for j in range(i+1, len(hlist)):
        pair_list.append([hlist[i], hlist[j]])
pair_list = np.array(pair_list)
chunk_size = 100
bin_size = comm.bin_size
distances = np.arange(bound[0], bound[1], bin_size)
pr = np.zeros(len(distances))
occ = np.zeros(len(pair_list))
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    dist = md.compute_distances(traj, pair_list)
    for idx, i in enumerate(np.amax(dist, axis = 1)):
        if i >= 1.87-bin_size/2 and i < 1.87+bin_size/2:
            print 'chunk_index and sequence:', chunk_index, idx
    temp, _ = np.histogram(np.argmax(dist, axis = 1), bins = range(len(dist[0])+1))
    occ += temp
    temp, _ = np.histogram(np.amax(dist, axis = 1), bins = np.append(distances, max(distances)+bin_size))
    pr += temp
    print chunk_index
#    if chunk_index == 10:
#        break

frames_read = (chunk_index+1)*chunk_size
print frames_read

########## Outputs ##########
fig, (ax1, ax2) = plt.subplots(2, figsize = (12, 12), dpi = 1000)
pr /= sum(pr)*bin_size
np.savetxt('dist.txt', np.column_stack((distances, pr)), fmt=['%12.4f', '%12.4f'])
ax1.plot(distances, pr)
ax1.set_xlim([bound[0], bound[1]])
ax2.plot(range(len(dist[0])), occ)
plt.savefig('1st.pdf')
print pair_list[np.argmax(occ)]
