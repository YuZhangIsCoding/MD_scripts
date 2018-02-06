#!/Users/yuzhang/anaconda/bin/python
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
temp = []
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_dir = comm.calc_xyz_direction(res_targ, frame, topology, para, direction)
        for item in xyz_dir:
            for i, subi in enumerate(item):
#                if subi < 1.4:
                if subi > 8.2:
                    temp.append(i)
        break
    break
mylist = [res_targ[0][i] for i in temp]
myfile = open('pairlist_neg.txt', 'w')
for i in mylist:
    count = 0
    for atom in topology.residue(i).atoms:
        if atom.name == 'N2' or atom.name == 'C3':
            count += 1
            myfile.write('%d\t' %(atom.index+1))
    if count == 2:
        count = 2
        myfile.write('\n')
myfile.close()
