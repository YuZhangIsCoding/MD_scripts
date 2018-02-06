#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: calc_density_planar_iter.py
# Description:  This is a python script to get the radial density profiles in 
#               Onion like carbon system with both cathode and anode. The
#               electrodes were arranged in the z direction. The lower one is
#               anode and upper one is cathode.
# Dates:    01-18-2018 Created

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

num_density = [np.zeros((len(distances), len(res_targ))) for _ in range(2)]
mass_density = [np.zeros(len(distances)) for _ in range(2)]
charge_density = [np.zeros(len(distances)) for _ in range(2)]

chunk_size = 10
center = [[4.386006683, 4.392999283, 4.37500211], [4.386006683, 4.392999283, 13.35800205]]
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'begin.gro')):
    for sub_ind, frame in enumerate(traj):
        for block in range(2):
            if block == 0:
                bool_ind = frame.xyz[0][:,2] < box[2]/2
            else:
                bool_ind = frame.xyz[0][:,2] >= box[2]/2
            dist2 = np.sum((frame.xyz[0][bool_ind]-center[block])**2, axis = 1)
            dist = np.sqrt(dist2)
            temp, _ = np.histogram(dist, weights = massind[bool_ind]/(4*np.pi*dist2),\
                    bins = np.append(distances, max(distances)+bin_size))
            mass_density[block] += temp

            temp, _ = np.histogram(dist, weights = chargeind[bool_ind]/(4*np.pi*dist2),\
                    bins = np.append(distances, max(distances)+bin_size))
            charge_density[block] += temp
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        for j, item in enumerate(xyz_com):
            for block in range(2):
                if block == 0:
                    dist2 = np.sum((item[item[:, 2] < box[2]/2]-center[block])**2, axis = 1)
                else:
                    dists = np.sum((item[item[:, 2] >= box[2]/2]-center[block])**2, axis = 1)
                dist = np.sqrt(dist2)
                temp, _ = np.histogram(dist, weights = 1/(4*np.pi*dist2), bins = np.append(distances, max(distances)+bin_size))
                num_density[block][:, j] += temp
        if sub_ind%10 == 9:
            print('Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1)
    if comm.args.test == True and chunk_index == 0:
        break

frames_read = chunk_index*chunk_size+sub_ind+1
print(frames_read)
mass_density = [item/(bin_size*frames_read/10*6.022) for item in mass_density]
charge_density = [item/(bin_size*frames_read) for item in charge_density]
num_density = [item/(bin_size*frames_read) for item in num_density]

########## Outputs ##########
outnames = ['Anode', 'Cathode']
for block in range(2):
    comm.out_nd(distances, num_density[block], res_name, out_name = 'NumDen_'+outnames[block])
    comm.out_density(distances, mass_density[block], charge_density[block], out_name = 'Den_'+outnames[block])
