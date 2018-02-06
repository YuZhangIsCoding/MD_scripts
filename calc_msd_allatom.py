#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_msd.py
# Description: This is a python script to calculate the MSD and diffusion coefficients using Einstein equation
# Dates:    09-20-2016 Created
#           01-17-2017 Changed with nframe/2 origins and nframe-nframe/2+1 time steps

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

if comm.args.suffix != None:
    outname = 'msd_'+comm.args.suffix
else:
    outname = 'msd'
########## Main ###########

traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
bin_size = comm.bin_size
dt = comm.dt
nframe = comm.args.nframe

atom_targ = []
for res_list in res_targ:
    atom_targ.append([])
    for res in res_list:
        atom_targ[-1].extend([i.index for i in topology.residue(res).atoms])

def get_xyz(frame, atoms):
    '''Return a list of arrays of the coordinates of different types of molecules
    '''
    xyz = []
    for atom_list in atoms:
        xyz.append([frame.xyz[0, i, :] for i in atom_list])
        xyz[-1] = np.array(xyz[-1])
    return xyz
chunk_size = 100
xyz_pre = get_xyz(traj0, atom_targ)
xyz_ref = [xyz_pre]
per = [[ np.zeros((len(xyz_pre[i]), 3))] for i in range(len(xyz_pre))]
sd = [[[] for row in range(nframe-nframe/2+1)] for i in range(len(res_targ))]
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        frame_ind = chunk_index*chunk_size+sub_ind
        if frame_ind == 0:
            continue
        box = frame.unitcell_lengths[0, :]
        xyz_com = get_xyz(frame, atom_targ)
        if frame_ind < nframe/2:
            xyz_ref.append(xyz_com)
            p_start = 0
            p_end = frame_ind
        elif nframe-nframe/2-frame_ind >= 0:
            p_start = 0
            p_end = nframe/2
        else:
            p_start = frame_ind-nframe+nframe/2
            p_end = nframe/2
            xyz_ref[p_start-1] = 0
        for res_ind, xyz_res in enumerate(xyz_com):
            per[res_ind].append(per[res_ind][-1].copy())
            for dim in range(3):
                for i in range(len(xyz_res)):
                    if xyz_res[i, dim]-xyz_pre[res_ind][i, dim] > box[dim]/2:
                        per[res_ind][-1][i, dim] -= box[dim]
                    elif xyz_res[i, dim]-xyz_pre[res_ind][i, dim] < -box[dim]/2:
                        per[res_ind][-1][i, dim] += box[dim]
            for loop_ind in range(p_start, p_end):
                disp = xyz_res-xyz_ref[loop_ind][res_ind]+per[res_ind][-1]-per[res_ind][loop_ind]
                sd[res_ind][frame_ind-loop_ind].append(np.sum(disp**2))
        xyz_pre = xyz_com
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
    if comm.args.test and chunk_index == 9:
        break

time = np.arange(nframe-nframe/2+1)*dt
data = time.copy()
pt_start = len(time)/2
msd = []

print 'frames read:', chunk_index*chunk_size+sub_ind+1
for i, item in enumerate(sd):
    sd[i][0] = 0
    msd.append([])
    for j in item:
        msd[-1].append(np.mean(j))
    msd[i] = np.array(msd[i])/len(atom_targ[i])
    coef = np.polyfit(time[pt_start:], msd[i][pt_start:], 1)
    print coef/6/1e6
    fit = np.poly1d(coef)
    plt.plot(time, msd[i], label = str(i))
    plt.plot(time, fit(time), '--', label = 'fit'+str(i))
    data = np.column_stack((data, msd[i]))
plt.legend(loc = 'best')
plt.xlabel('time (ps)')
plt.ylabel('MSD (nm^2)')
plt.savefig(outname+'.pdf')
np.savetxt(outname+'.txt', data, fmt = ['%12.6f' for i in range(len(data[0]))])
