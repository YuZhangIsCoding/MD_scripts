#!/Users/yuzhang/anaconda/bin/python
# Filename: 
# Date: 01-20-2015 Created
# Description: This is a file to compute the number density near the gold disk

import calc_common as calc_comm
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import pdb

########## main ##########
traj_name, traj0, topology = calc_comm.load_traj()
para = calc_comm.load_para()
res_targ, res_name = calc_comm.select_groups(traj0)
direction = calc_comm.direction
bin_size = calc_comm.bin_size
bound = calc_comm.load_bound(traj0, direction)
distances = np.arange(bound[0], bound[1], bin_size)
nd = np.zeros((len(distances), len(res_targ)))
chunk_size = 10
cyl = [4, 4, 1.5]
for chunk_index, traj in enumerate(md.iterload(traj_name,chunk = chunk_size, top = 'topol.pdb')):
    #if chunk_index == 1:
    #    break
    for i, frame in enumerate(traj):
        #xyz_com = calc_comm.calc_xyz_com(res_targ, frame, topology, para)
        xyz_com = calc_comm.calc_xyz_cyl(res_targ, frame, topology, para, cyl)
        for res_type, xyz in enumerate(xyz_com):
            temp, _ = np.histogram(xyz[:,direction], bins = np.append(distances, max(distances)+bin_size))
            nd[:,res_type] += temp
        print 'Reading',chunk_index*chunk_size+i+1
frames_read = chunk_index*chunk_size+i+1
#bin_volume = np.prod(traj.unitcell_lengths[0][:3])*bin_size/traj.unitcell_lengths[0][direction]
bin_volume = np.pi*cyl[2]**2*bin_size

plt.figure(figsize=(16, 12),dpi=1000)
for i in range(len(res_targ)):
    nd[:,i] /= frames_read*bin_volume 
    plt.plot(distances-0.682-0.07, nd[:,i])
plt.savefig('see.pdf')
fmt = ['%12.4f' for row in range(1+len(res_targ))]
np.savetxt('ND.txt', np.column_stack((distances,nd)), fmt=fmt)
