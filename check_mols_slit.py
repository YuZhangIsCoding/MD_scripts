#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_density_planar_iter.py
# Description: This is a python script to check how many ions remained within the slit during the course of simulation
# Dates:    01-17-2017 Created

import mdtraj as md
import matplotlib.pyplot as plt
import calc_common as comm


########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
direction = comm.direction
bound = comm.load_bound(traj0, direction)

count = [[[] for i in range(len(res_targ))] for row in range(len(bound)/2)]
chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        mol_rem = [[[] for i in range(len(res_targ))] for row in range(len(bound)/2)]
        temp = [[] for row in range(len(res_targ))]
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        for res_type, res_list in enumerate(res_targ):
            for i in range(len(bound)/2):
                count[i][res_type].append(0)
            for res_seq, res_ind in enumerate(res_list):
                for i in range(len(bound)/2):
                    if xyz_com[res_type][res_seq, direction] >= bound[2*i] and xyz_com[res_type][res_seq, direction] <= bound[2*i+1]:
                        temp[res_type].append(res_ind)
                        count[i][res_type][-1] += 1
                        mol_rem[i][res_type].append(res_ind)
                        continue
        res_targ = temp
        if sub_ind%10 == 9:
            print [len(i) for i in res_targ]
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if chunk_index == 0:
#        break

frames_read = chunk_index*chunk_size+sub_ind+1
print frames_read

########## Outputs ##########

myfile = open('mol_remain.txt', 'w')
for slit in mol_rem:
    for mols in slit:
        for mol in mols:
            myfile.write('%d\t' %mol)
        myfile.write('\n')
myfile.close()

for slit, mol in enumerate(count):
    for i, item in enumerate(mol):
        plt.plot(range(frames_read), item, label = 'slit '+str(slit+1)+'- mol '+str(i))
plt.legend()
plt.savefig('mol_t.pdf')
