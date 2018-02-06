#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_msd.py
# Description: This is a python script to calculate the MSD and diffusion coefficients using Einstein equation in select region
#              Since ions selected in specific dimension won't accross the boundary, some correction for periodicity can be negelected
# Dates:    09-29-2016 Created

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
bin_size = comm.bin_size
dt = comm.dt
nframe = comm.args.nframe

chunk_size = 100

def check_bound(xyz, bound):
    if xyz[0] > bound[0] and xyz[0] < bound[1] and xyz[1] > bound[2] and xyz[1] < bound[3] and xyz[2] > bound[4] and xyz[2] < bound[5]:
        return True
    else:
        return False


bound = [1, 4, 2, 4, 2, 4]

xyz_pre = comm.calc_xyz_com(res_targ, traj0, topology, para)
xyz_ref = []
xyz_ref.append(xyz_pre)

rem_stat = []
for xyz_res in xyz_pre:
    rem_stat.append([])
    for xyz in xyz_res:
        if check_bound(xyz, bound):
            rem_stat[-1].append(0)
        else:
            rem_stat[-1].append(-1)

per = [[ np.zeros((len(xyz_pre[i]), 3))] for i in range(len(xyz_pre))]
sd = [[[] for row in range(nframe-(nframe+1)/2+1)] for i in range(len(res_targ))]
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        frame_ind = chunk_index*chunk_size+sub_ind
        if frame_ind == 0:
            continue
        box = frame.unitcell_lengths[0, :]
        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)

        if frame_ind < (nframe+1)/2:
            xyz_ref.append(xyz_com)
            p_start = 0
            p_end = frame_ind
        else:
            p_start = frame_ind-nframe+(nframe+1)/2
            p_end = (nframe+1)/2

        for res_ind, xyz_res in enumerate(xyz_com):
            ## Correct for periodicity
            per[res_ind].append(per[res_ind][-1].copy())
            for dim in range(3):
                for i in range(len(xyz_res)):
                    if xyz_res[i, dim]-xyz_pre[res_ind][i, dim] > box[dim]/2:
                        per[res_ind][-1][i, dim] -= box[dim]
                    elif xyz_res[i, dim]-xyz_pre[res_ind][i, dim] < -box[dim]/2:
                        per[res_ind][-1][i, dim] += box[dim]
            for mol_ind, xyz in enumerate(xyz_res):
                if check_bound(xyz, bound):
                    rem_stat[res_ind][mol_ind] += 1
                else:
                    rem_stat[res_ind][mol_ind] = -1
                if rem_stat[res_ind][mol_ind] > 0:
                    if frame_ind-rem_stat[res_ind][mol_ind] >= (nframe+1)/2:
                        continue
                    else:
                        p_start = max(p_start, frame_ind-rem_stat[res_ind][mol_ind])
                    if frame_ind-rem_stat[res_ind][mol_ind] >= (nframe+1)/2:
                        print 'here'
                    for loop_ind in range(p_start, p_end):
                        disp = xyz - xyz_ref[loop_ind][res_ind][mol_ind]+per[res_ind][-1][mol_ind]-per[res_ind][loop_ind][mol_ind]
                        sd[res_ind][frame_ind-loop_ind].append(np.sum(disp**2))
        xyz_pre = xyz_com
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if chunk_index == 0:
#        break

#        box = frame.unitcell_lengths[0, :]
#        xyz_com = comm.calc_xyz_com(res_targ, frame, topology, para)
#        for idx_res, xyz_res in enumerate(xyz_com):
#            for i in range(3):
#                for j in range(len(xyz_res)):
#                    if xyz_res[j, i] - xyz_pre[idx_res][j, i] > box[i]/2:
#                        per[idx_res][j, i] -= 1
#                    elif xyz_res[j, i] - xyz_pre[idx_res][j, i] < -box[i]/2:
#                        per[idx_res][j, i] += 1
#            disp = xyz_res-xyz_ref[idx_res]+box*per[idx_res]
#            sd[idx_res].append(np.sum(disp**2))
#            msd[idx_res].append(np.mean(sd[idx_res]))
#        xyz_pre = xyz_com
time = np.arange(nframe-(nframe+1)/2+1)*dt
data = time.copy()
pt_start = int(len(time)*0.4)
pt_end = int(len(time)*0.8)
msd = []
print 'frames read:', chunk_index*chunk_size+sub_ind+1
for i, item in enumerate(sd):
    sd[i][0].append(0)
    msd.append([])
    count = min([len(j) for j in item[1:]])
    for j in item:
        msd[-1].append(np.mean(j[:count]))
    coef = np.polyfit(time[pt_start:pt_end], msd[i][pt_start:pt_end], 1)
    print coef/6
    fit = np.poly1d(coef)
    plt.plot(time, msd[i], label = str(i))
    plt.plot(time, fit(time), '--', label = 'fit'+str(i))
    data = np.column_stack((data, msd[i]))
plt.legend(loc = 'best')
plt.xlabel('time (ps)')
plt.ylabel('MSD (nm^2)')
plt.savefig('msd_sel.pdf')
np.savetxt('msd_sel.txt', data, fmt = ['%12.4f' for i in range(len(data[0]))])
