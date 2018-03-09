#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: check_mols_sphere.py
# Description:  A python script that compute the remaining ions in the cutoff
#               distance from a spherical surface
# Date:     03-09-2018  Created

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import calc_common as comm

traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)

if comm.args.center:
    center = comm.args.center
else:
    center = (4.48300466, 4.48300138, 4.48199593)

res_targ = [np.array(item) for item in res_targ]

if comm.args.v:
    print('The center is located at', center)
    print('Radius:', comm.args.exb)
    print('Cut off distance from surface:', comm.args.cutoff)
    print('Final molecules are writing to file %s...' %"mol_remain.txt")

chunk_size = comm.args.chunk
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk_size, top= 'begin.gro')):
    for sub_ind, frame in enumerate(traj):
        mol_rem = [[] for _ in res_targ]
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        temp = []
        for resid, res_list in enumerate(res_targ):
            if len(res_list):
                dist2 = np.sum((xyz_com[resid]-center)**2, axis = 1)
                temp.append(res_targ[resid][dist2 <= (comm.args.exb+comm.args.cutoff)**2])
            else:
                temp.append([])
        res_targ = temp
        if comm.args.v and sub_ind%10 == 9: 
            print('Molecules remained:', [len(i) for i in res_targ], end = ", ")
            print('Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1)
        if comm.args.test and comm.args.test == sub_ind+chunk_index*chunk_size:
            break
    else:
        continue
    break

np.savetxt('mol_remain.txt', res_targ, fmt = '%s')
outfile = open('mol_remain.txt', 'w')
for res in res_targ:
    for item in res:
        outfile.write('%d\t' %item)
    outfile.write('\n')
outfile.close()

