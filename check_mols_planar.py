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
#bound = comm.load_bound(traj0, direction)
#bound = [1.38, 1.61, 1.25, 4.98, 4.76, 5.11]
#bound_1 = [1.38, 1.61, 1.25, 4.98, 4.76, 5.11]
#bound_2 = [1.73, 1.96, 1.59, 4.63, 4.39, 4.76]
#bound_2 = [2.38, 2.61, 2.25, 3.98, 3.76, 4.11]
bound_1 = [1.75+0.682 for _ in range(6)]
bound_2 = [3.25+0.682 for _ in range(6)]

### OH functionalized
#bound = [1.39, 1.20, 1.24, 4.98, 5.15, 5.12]

count = [[] for _ in res_targ]
chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        mol_rem = [[] for _ in res_targ]
        temp = [[] for _ in res_targ]
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        for res_type, res_list in enumerate(res_targ):
            count[res_type].append(0)
            for res_seq, res_ind in enumerate(res_list):
#                if xyz_com[res_type][res_seq, direction] <= bound[res_type] or xyz_com[res_type][res_seq, direction] >= bound[3+res_type]:
                if xyz_com[res_type][res_seq, direction] >= bound_1[res_type] and xyz_com[res_type][res_seq, direction] <= bound_2[3+res_type]:
#                if (xyz_com[res_type][res_seq, direction] >= bound_1[res_type] and xyz_com[res_type][res_seq, direction] <= bound_2[res_type]) or (xyz_com[res_type][res_seq, direction] > bound_2[3+res_type] and xyz_com[res_type][res_seq, direction] < bound_1[3+res_type]):
#                if xyz_com[res_type][res_seq, direction] >= bound[res_type] and xyz_com[res_type][res_seq, direction] <= bound[3+res_type]:
                    temp[res_type].append(res_ind)
                    count[res_type][-1] += 1
                    mol_rem[res_type].append(res_ind)
                    continue
        res_targ = temp
        if sub_ind%10 == 9:
            print 'Molecules remained:', [len(i) for i in res_targ], 
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if chunk_index == 0:
#        break

frames_read = chunk_index*chunk_size+sub_ind+1
print frames_read

########## Outputs ##########

myfile = open('mol_remain.txt', 'w')
for mols in mol_rem:
    for mol in mols:
        myfile.write('%d\t' %mol)
    myfile.write('\n')
myfile.close()

for i, item in enumerate(count):
    plt.plot(range(frames_read), item, label = 'mol '+str(i+1))
plt.legend()
plt.savefig('mol_t.pdf')
