#!/Users/yuzhang/anaconda3/bin/python
# Filename: calc_density_planar_iter.py
# Description: This is a python script to get the density profiles of a channel simulation
# Dates:    10-08-2015 Created
#           02-05-2015 Modified to split the trajectory by chunks to save memory

import mdtraj as md
import numpy as np
import calc_common as comm

########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
direction = comm.direction
bin_size = comm.bin_size
bound = comm.load_bound(traj0, direction)
distances = np.arange(bound[0], bound[1], bin_size)

num_density = np.zeros((len(distances), len(res_targ)))
mass_density = np.zeros(len(distances))
charge_density = np.zeros(len(distances))

chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'begin.gro')):
    for sub_ind, frame in enumerate(traj):
        temp, _ = np.histogram(frame.xyz[0, :, direction], weights = massind, bins = np.append(distances, max(distances)+bin_size))
        mass_density += temp
        temp, _ = np.histogram(frame.xyz[0, :, direction], weights = chargeind, bins = np.append(distances, max(distances)+bin_size))
        charge_density += temp
        xyz_dir = comm.calc_xyz_direction(res_targ, frame, topology, para, direction)
        for j, item in enumerate(xyz_dir):
            temp, _ = np.histogram(item, bins = np.append(distances, max(distances)+bin_size))
            num_density[:, j] += temp
        if sub_ind%10 == 9:
            print('Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1)
    if comm.args.test == True and chunk_index == 0:
        break

frames_read = chunk_index*chunk_size+sub_ind+1
print(frames_read)
mass_density /= bin_size*np.prod(traj.unitcell_lengths[0])/traj.unitcell_lengths[0][direction]*frames_read/10*6.022
charge_density /= bin_size*np.prod(traj.unitcell_lengths[0])/traj.unitcell_lengths[0][direction]*frames_read
num_density /= bin_size*np.prod(traj.unitcell_lengths[0])/traj.unitcell_lengths[0][direction]*frames_read

########## Outputs ##########
comm.out_nd(distances, num_density, res_name)
comm.out_density(distances, mass_density, charge_density)
