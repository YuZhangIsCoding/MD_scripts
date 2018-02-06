#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_energy_nb.py
# Description: This is a python script to calculate the nonbond interaction between ions and the electrode surface within a cutoff distance
# Date: 02-24-2016 Created

import pdb, time
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

bound_slit = [4, 8, 23, 27]
#bound_slit = [4, 5, 23, 24]

############ Main #############
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name, wallindex, wallname = comm.select_groups(traj0)
direction = comm.direction
#bound = comm.load_bound(traj0, direction)
bin_size  = comm.bin_size
cutoff = comm.args.cutoff

#res_targ = [i[:100] for i in res_targ]

distances = []
n_dist = []
ninslit = []
vdw_sum = []
es_sum = []
vdw_slit = []
es_slit = []

p_start = 1.75
dinslit = np.arange(p_start, p_start+2.6, bin_size)
for i in range(2):
    distances.append(np.arange(bound_slit[2*i], bound_slit[2*i+1], bin_size))
    n_dist.append(np.zeros((len(distances[i]), len(res_targ))))
    ninslit.append(np.zeros((len(dinslit), len(res_targ))))
    vdw_sum.append(np.zeros((len(distances[i]), len(res_targ))))
    es_sum.append(np.zeros((len(distances[i]), len(res_targ))))
    vdw_slit.append(np.zeros((len(dinslit), len(res_targ))))
    es_slit.append(np.zeros((len(dinslit), len(res_targ))))
#wallindex = [range(200, 1000), range(12000, 12800)]
chunk_size = 50
pre_time = time.time()
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_com, vdw, es = comm.calc_nb_slit(res_targ, wallindex, frame, topology, para, direction, cutoff, bound_slit)
        for i in range(2):
            for j in range(len(res_targ)):
                if len(xyz_com[i][j]) == 0:
                    continue
                xyz_com[i][j] = np.array(xyz_com[i][j])
                temp, _ = np.histogram(xyz_com[i][j][:, direction], bins = np.append(distances[i], max(distances[i])+bin_size))
                n_dist[i][:, j] += temp
                temp, _ = np.histogram(xyz_com[i][j][:, 2], bins = np.append(dinslit, max(dinslit)+bin_size))
                ninslit[i][:, j] += temp
                temp, _ = np.histogram(xyz_com[i][j][:, direction], weights = vdw[i][j], bins = np.append(distances[i], max(distances[i])+bin_size))
                vdw_sum[i][:, j] += temp
                temp, _ = np.histogram(xyz_com[i][j][:, direction], weights = es[i][j], bins = np.append(distances[i], max(distances[i])+bin_size))
                es_sum[i][:, j] += temp
                temp, _ = np.histogram(xyz_com[i][j][:, 2], weights = vdw[i][j], bins = np.append(dinslit, max(dinslit)+bin_size))
                vdw_slit[i][:, j] += temp
                temp, _ = np.histogram(xyz_com[i][j][:, 2], weights = es[i][j], bins = np.append(dinslit, max(dinslit)+bin_size))
                es_slit[i][:, j] += temp
        print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
        print 'time',time.time()-pre_time
    if chunk_index == 0:
        break

dinslit -= p_start
for i in range(2):
    for j in range(len(res_targ)):
        for k, item in enumerate(n_dist[i][:, j]):
            if item != 0:
                vdw_sum[i][k, j] /= item
                es_sum[i][k, j] /= item
        for k, item in enumerate(ninslit[i][:, j]):
            if item != 0:
                vdw_slit[i][k, j] /= item
                es_slit[i][k, j] /= item
########## Output ##########
for i in range(2):
    np.savetxt('energy_vdw_slit%d_x_%s.txt' %(i, time.strftime("%d%m%y")), np.column_stack((distances[i], vdw_sum[i])), fmt = ['%12.4f' for row in range(len(res_targ)+1)])
    np.savetxt('energy_es_slit%d_x_%s.txt' %(i, time.strftime("%d%m%y")), np.column_stack((distances[i], es_sum[i])), fmt = ['%12.4f' for row in range(len(res_targ)+1)])
    np.savetxt('energy_vdw_slit%d_z_%s.txt' %(i, time.strftime("%d%m%y")), np.column_stack((dinslit, vdw_slit[i])), fmt = ['%12.4f' for row in range(len(res_targ)+1)])
    np.savetxt('energy_es_slit%d_z_%s.txt' %(i, time.strftime("%d%m%y")), np.column_stack((dinslit, es_slit[i])), fmt = ['%12.4f' for row in range(len(res_targ)+1)])

fig, (ax1, ax2) = plt.subplots(2, sharex= True)
for i in range(len(res_targ)):
    ax1.plot(distances[0], vdw_sum[0][:, i], label = res_name[i]+' vdw - slit 1')
    ax1.plot(distances[1], vdw_sum[1][:, i], label = res_name[i]+' vdw - slit 2')
    ax2.plot(distances[0], es_sum[0][:, i], label = res_name[i]+' es - slit 1')
    ax2.plot(distances[1], es_sum[1][:, i], label = res_name[i]+' es - slit 2')
ax1.legend(loc = 'best')
ax2.legend(loc = 'best')
plt.savefig('Energy_split_x.pdf')

fig, (ax1, ax2) = plt.subplots(2, sharex= True)
for i in range(len(res_targ)):
    ax1.plot(dinslit, vdw_slit[0][:, i], label = res_name[i]+' vdw - slit 1')
    ax1.plot(dinslit, vdw_slit[1][:, i], label = res_name[i]+' vdw - slit 2')
    ax2.plot(dinslit, es_slit[0][:, i], label = res_name[i]+' es - slit 1')
    ax2.plot(dinslit, es_slit[1][:, i], label = res_name[i]+' es - slit 2')
ax1.legend(loc = 'best')
ax2.legend(loc = 'best')
plt.savefig('Energy_split_z.pdf')
