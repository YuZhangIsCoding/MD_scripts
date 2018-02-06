#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_rdf.py
# Description: This is a python script to calculate the cylindrical distribution between different spicies
# Date: 01-28-2016: Created

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

############ main ##########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
direction = comm.direction   ## the direction where the slit are selected
v_dir = 2   ## the direction perpendicular to the slit surface
hght = 0.8
bin_size = 0.01
bound = comm.load_bound(traj0, direction)
#dim_flat = []
#for i in range(3):
#    if i != v_dir:
#        dim_flat.append(i)
#d_plane = [traj0.unitcell_lengths[0, i] for i in dim_flat]
box = traj0.unitcell_lengths[0, :]
cutoff = 3.0
distances = np.arange(0, cutoff+bin_size, bin_size)

# Antoher way is to use the previous function select groups, and automatically make pairs based on the selected groups.
# This redueces the lines of the code and make the codes more consistent, but it has less freedom to tune.
# Current version is based on the following code.
res_targ, res_name = comm.select_groups(traj0)
chunk_size  = 200
cyld = np.zeros((len(distances), sum(range(len(res_targ)+1))))
r2d = []
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        xyz_com = comm.calc_xyz_com_con(res_targ, frame, topology, para, direction, bound)
        r2d = []
        count = []
        for res_type, xyz_res in enumerate(xyz_com):
            for ref_type, xyz_ref in enumerate(xyz_com[res_type:]):
                count.append(len(xyz_res))
                r2d.append([])
                for i in xyz_res:
                    for j in xyz_ref:
                        if sum((i-j)**2) != 0:
                            temp = np.abs(i-j)
                            tt = 0
                            for subdim in range(3):
                                if subdim != v_dir:
                                    if subdim != direction:
                                        if temp[subdim] >= box[subdim]/2:
                                            temp[subdim] = box[subdim]-temp[subdim]
                                    tt += temp[subdim]**2
                            r2d[-1].append(np.sqrt(tt))
        for ind, item in enumerate(r2d):
            temp, _ = np.histogram(item, bins = np.append(distances, cutoff+2*bin_size))
            cyld[:, ind] += temp/float(count[ind])
        print 'Reading',chunk_index*chunk_size+sub_ind+1
    if comm.args.test == True:
        if chunk_index == 0:
            break

frames_read = chunk_index*chunk_size+sub_ind+1
########## output ###########
label_name = []
for i in range(len(res_targ)):
    for j in range(i, len(res_targ)):
        label_name.append(res_name[i]+'-'+res_name[j])

for i, r in enumerate(distances):
    if r != 0:
        cyld[i, :] /= (2*np.pi*r*hght*bin_size*frames_read)
for i in range(len(cyld[0])):
    plt.plot(distances, cyld[:, i], label = label_name[i])
plt.legend()
plt.savefig('RDF_2d.pdf')

fmt = ['%12.4f' for row in range(1+len(cyld[0]))]
np.savetxt('RDF_2d.txt', np.column_stack((distances, cyld)), fmt = fmt)
