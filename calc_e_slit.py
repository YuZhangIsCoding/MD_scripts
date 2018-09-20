#!/Users/yuzhang/anaconda3/bin/python
# Filename: calc_density_planar_iter.py
# Description: This is a python script to get the charge density profiles inside a slit pore
# Dates:    08-24-2015 Created

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

# hard code
dir_bound = 2
dir_result = 0
distances = np.arange(0, 10, bin_size)

mass_density = np.zeros(len(distances))
charge_density = np.zeros(len(distances))

chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'begin.gro')):
    for sub_ind, frame in enumerate(traj):
        indicator = (traj0.xyz[0, :, dir_bound] > bound[0]) & (traj0.xyz[0, :, dir_bound]< bound[1]).astype(int)
        temp, _ = np.histogram(frame.xyz[0, :, dir_result], weights = massind*indicator,
                               bins = np.append(distances, max(distances)+bin_size))
        mass_density += temp
        temp, _ = np.histogram(frame.xyz[0, :, dir_result], weights = chargeind*indicator,
                               bins = np.append(distances, max(distances)+bin_size))
        charge_density += temp
        if sub_ind%10 == 9:
            print('Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1)
    if comm.args.test == True and chunk_index == 0:
        break

frames_read = chunk_index*chunk_size+sub_ind+1
print(frames_read)
mass_density /= bin_size*np.prod(traj.unitcell_lengths[0])/traj.unitcell_lengths[0][dir_result]*frames_read/10*6.022
charge_density /= bin_size*np.prod(traj.unitcell_lengths[0])/traj.unitcell_lengths[0][dir_result]*frames_read
comm.out_density(distances, mass_density, charge_density)
