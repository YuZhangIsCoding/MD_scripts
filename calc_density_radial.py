#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: calc_density_planar_iter.py
# Description:  This is a python script to get the radial density profiles in 
#               Onion like carbon system.
# Dates:    11-27-2017 Created

import mdtraj as md
import numpy as np
import calc_common as comm

########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
bin_size = comm.bin_size

box = traj0.unitcell_lengths[0]
distances = np.arange(0, np.min(box/2), bin_size)

num_density = np.zeros((len(distances), len(res_targ)))
mass_density = np.zeros(len(distances))
charge_density = np.zeros(len(distances))

chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'begin.gro')):
    for sub_ind, frame in enumerate(traj):
        dist2 = np.sum((frame.xyz[0]-box/2)**2, axis = 1)
        dist = np.sqrt(dist2)
        mass_weights = massind/(4*np.pi*dist2)
        charge_weights = chargeind/(4*np.pi*dist2)
        temp, _ = np.histogram(dist, weights = mass_weights, bins = np.append(distances, max(distances)+bin_size))
        mass_density += temp
        temp, _ = np.histogram(dist, weights = charge_weights, bins = np.append(distances, max(distances)+bin_size))
        charge_density += temp
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        for j, item in enumerate(xyz_com):
            dist2 = np.sum((item-box/2)**2, axis = 1)
            dist = np.sqrt(dist2)
            weights = 1/(4*np.pi*dist2)
            temp, _ = np.histogram(dist, weights = weights, bins = np.append(distances, max(distances)+bin_size))
            num_density[:, j] += temp
        if sub_ind%10 == 9:
            print('Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1)
    if comm.args.test == True and chunk_index == 0:
        break

frames_read = chunk_index*chunk_size+sub_ind+1
print(frames_read)
mass_density /= bin_size*frames_read/10*6.022
charge_density /= bin_size*frames_read
num_density /= bin_size*frames_read

########## Outputs ##########
comm.out_nd(distances, num_density, res_name)
comm.out_density(distances, mass_density, charge_density)
