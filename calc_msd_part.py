#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_msd.py
# Description:  This is a python script to calculate the MSD and diffusion coefficients using Einstein equation
#               To calculate the diffusion coefficients, this script used the initial coordinate of a particle
#               to assign its box
# Dates:    09-20-2016 Created
#           03-30-2017 Modified the code so that different molecules can use different corresponding different boxes

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

if comm.args.suffix != None:
    outname = 'msd_'+comm.args.suffix
else:
    outname = 'msd'

def check_bound(xyz, bound):
    if xyz[0] > bound[0] and xyz[0] < bound[1] and xyz[1] > bound[2] and xyz[1] < bound[3] and xyz[2] > bound[4] and xyz[2] < bound[5]:
        return True
    else:
        return False

########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
#direction = comm.direction
bin_size = comm.bin_size
dt = comm.dt
nframe = comm.args.nframe
#bound = comm.load_bound(traj0, direction)

chunk_size = 100
xyz_pre = comm.calc_xyz_com(res_targ, traj0, topology, para)

### bound for testing bulk simulations ###
bound = [[0, 100, 0, 100, 0, 1.40],\
        [0, 100, 0, 100, 0, 1.22],\
        [0, 100, 0, 100, 0, 1.24]]

rem_stat = []
for i, xyz_res in enumerate(xyz_pre):
    rem_stat.append([])
    for j, xyz in enumerate(xyz_res):
        if check_bound(xyz, bound[i]):
            rem_stat[-1].append(res_targ[i][j])
    print len(rem_stat[-1]), 'of',len(res_targ[i]),'molecule', i+1, 'selected'

res_targ = rem_stat
xyz_pre = comm.calc_xyz_com(res_targ, traj0, topology, para)
xyz_ref = []
xyz_ref.append(xyz_pre)
#msd = [[0] for row in range(len(res_targ))]
per = [[ np.zeros((len(xyz_pre[i]), 3))] for i in range(len(xyz_pre))]
#for i in range(len(xyz_pre)):
#    per.append(np.zeros((len(xyz_pre[i]), 3)))
#sd = [[] for row in range(len(res_targ))]
#msd = [[] for row in range(len(res_targ))]
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
    if chunk_index == 0:
        if comm.args.test == True:
            break

time = np.arange(nframe-(nframe+1)/2+1)*dt
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
