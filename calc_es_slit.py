#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_energy_nb.py
# Description: This is a python script to calculate the nonbond interaction between ions and the electrode surface within a cutoff distance
# Date: 02-24-2016 Created

import pdb, time, sys
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

bound_slit = [4, 8, 23, 27]

#### functions ######
step = 0.1
if comm.args.dim == None:
    print 'No boxsize specified, gonna use the default value: (38, 3.19, 6.1)'
    dim = [38.0, 3.19, 6.1]
else:
        dim = comm.args.dim
box = dim
dim = [int(np.ceil(i*np.round(1/step))) for i in dim]

grid = np.loadtxt('grid_all.txt')

def est_pot(xyz, grid, step, dim, box):
    '''Estimate the potential at a specific coordinate according to the grids nearby'''
    grid0 = (xyz/step).astype(int)
    grid1 = grid0+1
    coord0 = grid0*step
    coord1 = grid1*step
    for i, item in enumerate(grid0):
        if item == dim[i]:
            grid1[i] = 0
            coord1[i] = box[i]
    diff0 = (xyz-coord0)/(coord1-coord0)
    diff = [1-diff0, diff0]
    temp = [grid0, grid1]
    grid_all = []
    pot = 0
    for i in range(2):
        for j in range(2):
            for k in range(2):
                label = temp[i][0]*dim[1]*dim[2]+temp[j][1]*dim[2]+temp[k][2]
                frac = np.prod([diff[i][0], diff[j][1], diff[k][2]])
                pot += frac*grid[label]
    return pot


############ Main #############
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
direction = comm.direction
#bound = comm.load_bound(traj0, direction)
bin_size  = comm.bin_size
cutoff = comm.args.cutoff

distances = []
n_dist = []
ninslit = []
es_sum = []
es_slit = []

p_start = 1.75
dinslit = np.arange(p_start, p_start+2.6, bin_size)
for i in range(2):
    distances.append(np.arange(bound_slit[2*i], bound_slit[2*i+1], bin_size))
    n_dist.append(np.zeros((len(distances[i]), len(res_targ))))
    ninslit.append(np.zeros((len(dinslit), len(res_targ))))
    es_sum.append(np.zeros((len(distances[i]), len(res_targ))))
    es_slit.append(np.zeros((len(dinslit), len(res_targ))))
#wallindex = [range(200, 1000), range(12000, 12800)]
chunk_size = 100
pre_time = time.time()
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_com = [[[] for row in range(len(res_targ))] for slit in range(len(bound_slit)/2)]
        pot_com = [[[] for row in range(len(res_targ))] for slit in range(len(bound_slit)/2)]
        para_me_res = [[] for row in range(len(res_targ))]
        for res_type, res_list in enumerate(res_targ):
            for atom in topology.residue(res_list[0]).atoms:
                para_me_res[res_type].append(para[atom.name][:2])
            para_me_res[res_type] = np.array(para_me_res[res_type])
            for res_index in res_list:
                xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
                for slit_ind in range(len(bound_slit)/2):
                    if max(xyz_org[:, 0]) > bound_slit[2*slit_ind] and min(xyz_org[:, 0] < bound_slit[2*slit_ind+1]):
                        for i, item in enumerate(xyz_org.T):
                            if np.ptp(item) > box[i]/2:
                                item[item < box[i]/2] += box[i]
                                xyz_org[:,i] = np.array(item).T
                        temp = np.sum(np.multiply(xyz_org.T, para_me_res[res_type][:,0]).T, axis=0)/sum(para_me_res[res_type][:,0])
                        if temp[0] >= bound_slit[2*slit_ind] and temp[0] <= bound_slit[2*slit_ind+1]:
                            for item in range(3):
                                if temp[item] > box[item]:
                                    temp[item] -= box[item]
                            xyz_com[slit_ind][res_type].append(temp)
                            pot = 0
                            for atom_ind, xyz in enumerate(xyz_org):
                                pot += est_pot(xyz, grid, step, dim, box)*para_me_res[res_type][atom_ind, 1]
                            pot_com[slit_ind][res_type].append(pot*1.6*6.022*10)
        for i in range(len(bound_slit)/2):
            for j in range(len(res_targ)):
                if len(xyz_com[i][j]) == 0:
                    continue
                xyz_com[i][j] = np.array(xyz_com[i][j])
                temp, _ = np.histogram(xyz_com[i][j][:, 0], bins = np.append(distances[i], max(distances[i])+bin_size))
                n_dist[i][:, j] += temp
                temp, _ = np.histogram(xyz_com[i][j][:, 2], bins = np.append(dinslit, max(dinslit)+bin_size))
                ninslit[i][:, j] += temp
                temp, _ = np.histogram(xyz_com[i][j][:, direction], weights = pot_com[i][j], bins = np.append(distances[i], max(distances[i])+bin_size))
                es_sum[i][:, j] += temp
                temp, _ = np.histogram(xyz_com[i][j][:, 2], weights = pot_com[i][j], bins = np.append(dinslit, max(dinslit)+bin_size))
                es_slit[i][:, j] += temp
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
            print 'time',time.time()-pre_time
#    if chunk_index == 0:
#        break

dinslit -= p_start
for i in range(2):
    for j in range(len(res_targ)):
        for k, item in enumerate(n_dist[i][:, j]):
            if item != 0:
                es_sum[i][k, j] /= item
        for k, item in enumerate(ninslit[i][:, j]):
            if item != 0:
                es_slit[i][k, j] /= item
########## Output ##########
for i in range(2):
    np.savetxt('es_slit%d_x.txt' %i, np.column_stack((distances[i], es_sum[i])), fmt = ['%12.4f' for row in range(len(res_targ)+1)])
    np.savetxt('es_slit%d_z.txt' %i, np.column_stack((dinslit, es_slit[i])), fmt = ['%12.4f' for row in range(len(res_targ)+1)])

fig = plt.figure()
for i in range(len(res_targ)):
    plt.plot(distances[0], es_sum[0][:, i], label = res_name[i]+' es - slit 1')
    plt.plot(distances[1], es_sum[1][:, i], label = res_name[i]+' es - slit 2')
plt.ylabel('kJ/mol')
plt.legend(loc = 'best')
plt.savefig('Es_x.pdf')

fig = plt.figure()
for i in range(len(res_targ)):
    plt.plot(dinslit, es_slit[0][:, i], label = res_name[i]+' es - slit 1')
    plt.plot(dinslit, es_slit[1][:, i], label = res_name[i]+' es - slit 2')
#plt.xlim([0.3, 0.5])
plt.ylabel('kJ/mol')
plt.legend(loc = 'best')
plt.savefig('Es_z.pdf')
