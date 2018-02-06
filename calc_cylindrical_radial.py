#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_cylindrical_distribution.py
# Description: This is a python script to calculate the cylindrical distribution between different spicies
# Date: 01-25-2016: Created

import pdb
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as calc_comm

############ main ##########
traj_name, traj0, topology = calc_comm.load_traj()
nonbond = calc_comm.load_nonbond()
para = calc_comm.load_para(nonbond)
direction = 2
bound = calc_comm.load_bound(traj0, direction)
bin_size = 0.01
dim_flat = []
for i in range(3):
    if i != direction:
        dim_flat.append(i)
d_plane = [traj0.unitcell_lengths[0, i] for i in dim_flat]
dimensions = [int(i/bin_size)+1 for i in d_plane]
cutoff = 3.0
distances = np.arange(0, cutoff+bin_size, bin_size)

# The following way creates a new function called radial_pairs in the calc_common.py file
#pair_targ = calc_comm.radial_pairs(traj0)

# Antoher way is to use the previous function select groups, and automatically make pairs based on the selected groups.
# This redueces the lines of the code and make the codes more consistent, but it has less freedom to tune.
# Current version is based on the followin code.
res_targ, res_name = calc_comm.select_groups(traj0)
#pair_targ = []
#for i in range(len(res_targ)):
#    for j in range(i, len(res_targ)):
#        pair_targ.append([res_targ[i], res_targ[j]])

n2d = [np.zeros(dimensions, dtype = np.int) for row in range(len(res_targ))]
chunk_size  = 1000
cyld = np.zeros((len(distances), sum(range(len(res_targ)+1))))
r2d = []
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    #if chunk_index == 1:
    #    break
    for sub_ind, frame in enumerate(traj):
        xyz_com = calc_comm.calc_xyz_com_con(res_targ, frame, topology, para, direction, bound)
        r2d = []
        count = []
        for res_type, xyz_res in enumerate(xyz_com):
            for ref_type, xyz_ref in enumerate(xyz_com[res_type:]):
                count.append(len(xyz_res))
                r2d.append([])
                for i in xyz_res[:, dim_flat]:
                    if ref_type == 0:
                        n2d[res_type][int(i[dim_flat[0]]/bin_size), int(i[dim_flat[1]]/bin_size)] += 1
                    for j in xyz_ref[:, dim_flat]:
                        if sum((i-j)**2) != 0:
                            temp = np.abs(i-j)
                            for subtemp in range(2):
                                if temp[subtemp] >= d_plane[subtemp]/2:
                                    temp[subtemp] = d_plane[subtemp]-temp[subtemp]
                            r2d[-1].append(np.sqrt(sum(temp**2)))
        for ind, item in enumerate(r2d):
            temp, _ = np.histogram(item, bins = np.append(distances, cutoff+2*bin_size))
            cyld[:, ind] += temp/float(count[ind])
        print 'Reading',chunk_index*chunk_size+sub_ind+1

frames_read = chunk_index*chunk_size+sub_ind+1
########## output ###########
label_name = []
for i in range(len(res_targ)):
    for j in range(i, len(res_targ)):
        label_name.append(res_name[i]+'-'+res_name[j])

for i, r in enumerate(distances):
    if r != 0:
        cyld[i, :] /= (2*np.pi*r*(bound[1]-bound[0])*frames_read*bin_size)
for i in range(len(cyld[0])):
    plt.plot(distances, cyld[:, i], label = label_name[i])
plt.legend()
plt.savefig('Cylindrical.pdf')

fmt = ['%12.4f' for row in range(1+len(cyld[0]))]
np.savetxt('Cylindrical_distri.txt', np.column_stack((distances, cyld)), fmt = fmt)
for i, item in enumerate(res_name):
    outfile = open('Plane_distri_%s.txt' %item, 'w')
    for row in n2d[i]:
        temp = [str(a) for a in row]
        outfile.write('%s\n' %(" ".join(temp)))
    outfile.close()
