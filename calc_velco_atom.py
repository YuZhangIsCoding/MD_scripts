#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: calc_msd.py
# Description:  This is a python script to calculate the velocity autocorrelation functions(VACF) and normalized VAF (NVACF)
#               VACF = <v(0)*v(t)>
#               NVACF = <v(0)*v(t)>/<v(0)*v(0)>
# Date:     03-22-2017 
#           03-23-2017  Output both velocity correlation functions and normalized velocity correlation functions

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm

if comm.args.suffix != None:
    outname = 'vcf_'+comm.args.suffix
else:
    outname = 'vcf'
########## Main ###########

traj_name, traj0, topology = comm.load_traj()
#para = comm.load_para()
#res_targ, res_name = comm.select_groups(traj0)
#massind, chargeind = comm.load_atomic(res_targ, topology, para)
bin_size = comm.bin_size
dt = comm.dt
nframe = comm.args.nframe-1

atom_targ = ['C8']
index = [topology.select('name %s' %item) for item in atom_targ]

xyz_pre = [traj0.xyz[0][item] for item in index]
vc = [[[] for row in range(nframe-int(nframe/2)+1)] for _ in atom_targ]
base = [[[] for row in range(nframe-int(nframe/2)+1)] for _ in atom_targ]
vel_ref = []
chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'begin.gro')):
    for sub_ind, frame in enumerate(traj):
        frame_ind = chunk_index*chunk_size+sub_ind
        if frame_ind == 0:
            continue
        box = frame.unitcell_lengths[0, :]
        xyz_com = [frame.xyz[0][item] for item in index]
        vel = [[] for _ in atom_targ]
        for res_ind, xyz_res in enumerate(xyz_com):
            vel[res_ind] = xyz_res-xyz_pre[res_ind]
            for dim in range(3):
                for i in range(len(xyz_res)):
                    if xyz_res[i, dim]-xyz_pre[res_ind][i, dim] > box[dim]/2:
                        vel[res_ind][i, dim] -= box[dim]
                    elif xyz_res[i, dim]-xyz_pre[res_ind][i, dim] < -box[dim]/2:
                        vel[res_ind][i, dim] += box[dim]
        xyz_pre = xyz_com
        if (frame_ind-1) < int(nframe/2):
            vel_ref.append(vel)
            p_start = 0
            p_end = frame_ind-1
            for res_ind, vel_res in enumerate(vel):
                base[res_ind][0].append(np.mean(np.sum(vel[res_ind]**2, axis = 1)))
                vc[res_ind][0].append(base[res_ind][0][-1])
        elif nframe-int(nframe/2)-(frame_ind-1) >= 0:
            p_start = 0
            p_end = int(nframe/2)
        else:
            p_start = (frame_ind-1)-nframe+int(nframe/2)
            p_end = int(nframe/2)
        for res_ind, vel_res in enumerate(vel):
            for loop_ind in range(p_start, p_end):
                vc[res_ind][frame_ind-1-loop_ind].append(np.mean(np.sum(vel[res_ind]*vel_ref[loop_ind][res_ind], axis = 1)))
                base[res_ind][frame_ind-1-loop_ind].append(np.mean(np.sum(vel_ref[loop_ind][res_ind]**2, axis = 1)))
        if sub_ind%10 == 9:
            print('Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1)
#    if chunk_index == 0:
#        break

time = np.arange(nframe-int(nframe/2)+1)*dt
data = time.copy()
data_raw = time.copy()
pt_start = len(time)/2
nvcf = []
vcf = []

print('frames read:', chunk_index*chunk_size+sub_ind+1)
fig = plt.figure()
for i, item in enumerate(vc):
    nvcf.append([])
    vcf.append([])
    for j in range(len(item)):
        nvcf[-1].append(np.mean(item[j])/np.mean(base[i][j]))
        vcf[-1].append(np.mean(item[j]))
    nvcf[i] = np.array(nvcf[i])
    vcf[i] = np.array(vcf[i])*1e6/dt**2  ### units are m^2/s^2
    data = np.column_stack((data, nvcf[i]))
    data_raw = np.column_stack((data_raw, vcf[i]))
    plt.plot(time, nvcf[i], label = str(i))
plt.legend(loc = 'best')
plt.xlabel('time (ps)')
plt.ylabel('Normalized Velociry Correlation Function')
plt.savefig('n'+outname+'.pdf')
np.savetxt('n'+outname+'.txt', data, fmt = ['%12.6f' for i in range(len(data[0]))])
np.savetxt(outname+'.txt', data_raw, fmt = ['%14.6E' for i in range(len(data_raw[0]))])
