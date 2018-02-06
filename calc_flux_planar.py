#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_density_planar_iter.py
# Description: This is a python script to get the density profiles of a channel simulation
# Dates:    10-08-2015 Created
#           02-05-2015 Modified to split the trajectory by chunks to save memory

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

chunk_size = 100
xyz_pre = comm.calc_xyz_com(res_targ, traj0, topology, para)
x_st = 1 # interface that flux cross
flux = [[] for _ in res_targ]
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        for res_ind, xyz_res in enumerate(xyz_com):
            flux_in = []
            flux_out = []
            for mol_ind, xyz in enumerate(xyz_res):
                if xyz[0] >= x_st and xyz_pre[res_ind][mol_ind, 0] < x_st:
                    flux_in.append(xyz[2])
                elif xyz[0] <= x_st and xyz_pre[res_ind][mol_ind, 0] > x_st:
                    flux_out.append(xyz[2])
            temp_in, _ = np.histogram(flux_in, bins = distances)
            temp_out, _ = np.histogram(flux_out, bins = distances)
            flux[res_ind].append(temp_in-temp_out)
        xyz_pre = xyz_com
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if comm.args.check == True:
#        if chunk_index == 0:
#            break
fig, axes = plt.subplots(len(res_targ), sharex = True)
for res in range(len(flux)):
    flux[res] = np.array(flux[res])
    for j in range(len(flux[res][0])):
        axes[res].plot(range(len(flux[res])), flux[res][:, j])
plt.savefig('flux_t.pdf')

fig, axes = plt.subplots(len(res_targ), sharex = True)
for res in range(len(flux)):
    axes[res].plot(range(len(flux[res][0])), np.mean(flux[res], axis = 0))
plt.savefig('flux.pdf')
