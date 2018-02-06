#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_msd.py
# Description: This is a python script to calculate the MSD and diffusion coefficients using Einstein equation
# Dates:    09-20-2016 Created
#           01-17-2017 Changed with nframe/2 origins

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

if comm.args.suffix != None:
    outname = 'msd_'+comm.args.suffix
else:
    outname = 'msd'

########## Main ###########

myfile = open('mol_remain.txt', 'r')
res_targ = []
for line in myfile:
    if line != '\n':
        res_targ.append([int(i) for i in line.split()])

traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
#res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
dt = comm.dt
nframe = comm.args.nframe
direction = comm.args.direction

chunk_size = 100
xyz_pre = comm.calc_xyz_com(res_targ, traj0, topology, para)
xyz_ref = []
xyz_ref.append(xyz_pre)
per = [[ np.zeros((len(xyz_pre[i]), 3))] for i in range(len(xyz_pre))]
sd = [[[] for row in range(nframe-nframe/2+1)] for i in range(len(res_targ))]
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        frame_ind = chunk_index*chunk_size+sub_ind
        if frame_ind == 0:
            continue
        box = frame.unitcell_lengths[0, :]
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
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
                disp = disp[:, :2]
                sd[res_ind][frame_ind-loop_ind].append(np.sum(disp**2))
        xyz_pre = xyz_com
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if chunk_index == 9:
#        break

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
    msd[i] = np.array(msd[i])/len(res_targ[i])
    coef = np.polyfit(time[pt_start:], msd[i][pt_start:], 1)
    print coef/4/1e6
    fit = np.poly1d(coef)
    plt.plot(time, msd[i], label = str(i))
    plt.plot(time, fit(time), '--', label = 'fit'+str(i))
    data = np.column_stack((data, msd[i]))
plt.legend(loc = 'best')
plt.xlabel('time (ps)')
plt.ylabel('MSD (nm^2)')
plt.savefig(outname+'.pdf')
np.savetxt(outname+'.txt', data, fmt = ['%12.4f' for i in range(len(data[0]))])
