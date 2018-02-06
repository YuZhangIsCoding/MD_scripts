#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_energy_nb.py
# Description: This is a python script to calculate the nonbond interaction between ions and the electrode surface within a cutoff distance
# Date: 02-24-2016 Created

import pdb, time
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

############ Main #############
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name, wallindex, wallname = comm.select_groups(traj0)
direction = comm.direction
bound = comm.load_bound(traj0, direction)
bin_size  = comm.bin_size
cutoff = comm.args.cutoff

distances = np.arange(0, cutoff+1.5, bin_size)

#ndist_bot = np.zeros((len(distances),3), dtype = np.int)
#ndist_top = np.zeros((len(distances),3), dtype = np.int)
#npotential = np.zeros((len(distances),3))

chunk_size = 200
n_dist = np.zeros((len(distances), len(res_targ)))
vdw_sum = np.zeros((len(distances), len(res_targ)))
es_sum = np.zeros((len(distances), len(res_targ)))
pre_time = time.time()
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_com, vdw, es = comm.calc_nb(res_targ, wallindex, frame, topology, para, direction, cutoff)
        for i, item in enumerate(xyz_com):
            temp, _ = np.histogram(item, bins = np.append(distances, max(distances)+bin_size))
            n_dist[:, i] += temp
            temp, _ = np.histogram(item, weights = vdw[i], bins = np.append(distances, max(distances)+bin_size))
            vdw_sum[:, i] += temp
            temp, _ = np.histogram(item, weights = es[i], bins = np.append(distances, max(distances)+bin_size))
            es_sum[:, i] += temp
        print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
        print 'time',time.time()-pre_time
    if chunk_index == 0:
        break
for i in range(len(res_targ)):
    for j, item in enumerate(n_dist[:, i]):
        if item != 0:
            vdw_sum[j, i] /= item
            es_sum[j, i] /= item
########## Output ##########
np.savetxt('energy_vdw_%s.txt' %time.strftime("%d%m%y") ,np.column_stack((distances, vdw_sum)), fmt = ['%12.4f' for row in range(len(res_targ)+1)])
np.savetxt('energy_es_%s.txt' %time.strftime("%d%m%y"),np.column_stack((distances, es_sum)), fmt = ['%12.4f' for row in range(len(res_targ)+1)])
fig, (ax1, ax2) = plt.subplots(2, sharex= True)
for i in range(len(res_targ)):
    ax1.plot(distances, vdw_sum[:, i], label = res_name[i]+' vdw')
    ax2.plot(distances, es_sum[:, i], label = res_name[i]+' electrostatic')
ax1.legend()
ax2.legend()
plt.savefig('Energy_split.pdf')
