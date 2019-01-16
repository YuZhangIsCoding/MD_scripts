#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_density_planar_iter.py
# Description:  This is a python script to get the time correlation function of 
#               ions/solvents with the planar surface
# Dates:    02-23-2017 Created

import sys
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm


########## Main ###########
traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj0)
massind, chargeind = comm.load_atomic(res_targ, topology, para)
direction = comm.direction
bin_size = comm.bin_size
bound = comm.load_bound(traj0, direction)
distances = np.arange(bound[0], bound[1], bin_size)

dt = comm.dt
nframe = comm.args.nframe
if nframe == 0:
    sys.exit('Please input the number of frames!')

if comm.args.suffix == None:
    outname = 'Time_correlation'
else:
    outname = 'Time_correlation_'+comm.args.suffix

## bound_conf is a list that specifies the boundary of the region
## Each species has upper limit (ul) and lower limit (ll)
## Example [ll1, ll2, ll3, ul1, ul2, ul3] for 3 molecule types.
#bound_conf = [1.12, 1.07, 1.02, 5.24, 5.27, 5.33]

## TEABF4 pristine
#bound_conf = [1.38, 1.61, 1.25, 4.98, 4.76, 5.11]
## TEABF4 OH
bound_conf = [1.39, 1.20, 1.24, 4.97, 5.17, 5.10]


#unknown
#bound_conf = [2.40, 2.42, 3.69, 3.66]

## This following code is for general analysis
## Similar to the way employed in the calculation of MSD, we use first half 
## of frames as the origins to effectively use the available data, and achieve
## good staticticle reliability.

chunk_size = 100
res_remain = [[] for _ in res_targ]
tc = [[[] for row in range(nframe-nframe/2+1)] for _ in res_targ]
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        frame_ind = chunk_index*chunk_size+sub_ind
        if frame_ind < nframe/2:
            p_start = 0
            p_end = frame_ind
        elif nframe-nframe/2-frame_ind >= 0:
            p_start = 0
            p_end = nframe/2
        else:
            p_start = frame_ind-nframe+nframe/2
            p_end = nframe/2
        xyz_dir = comm.calc_xyz_direction(res_targ, frame, topology, para, direction)
        for i in range(len(res_targ)):
            res_remain[i].append([])
            for j, item in enumerate(xyz_dir[i]):
                if item <= bound_conf[i] or item >= bound_conf[i+len(res_targ)]:
                    res_remain[i][-1].append(res_targ[i][j])
            temp = set(res_remain[i][-1])
            if frame_ind < nframe/2:
                tc[i][0].append(len(temp))
            for loop_ind in range(frame_ind-1, p_start-1, -1):
                temp = temp&set(res_remain[i][loop_ind])
                if loop_ind < p_end:
                    tc[i][frame_ind-loop_ind].append(len(temp))
        if sub_ind%10 == 9:
            print('Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1)
    if comm.args.test  and chunk_index == 0:
        break

frames_read = chunk_index*chunk_size+sub_ind+1
print(frames_read)

### The following block is for the cases when ions/solvents leaves the surface cage very rapidly
### It takes average of different trajectory chunks
#chunk_size = 100
#for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
#    res_remain = [item[:] for item in res_targ]
#    tc = [[] for i in range(len(res_targ))]
#    for sub_ind, frame in enumerate(traj):
#        xyz_dir = comm.calc_xyz_direction(res_remain, frame, topology, para, direction)
#        for i in range(len(res_remain)):
#            temp = res_remain[i][:]
#            #res_remain[i].append([])
#            for j, item in enumerate(xyz_dir[i]):
#                if item > bound_conf[i] and item < bound_conf[i+len(res_targ)]:
#                    temp.remove(res_remain[i][j])
#            res_remain[i] = temp[:]
#            tc[i].append(len(temp))
#        if sub_ind%10 == 9:
#            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
#    if 'tc_tot' not in globals():
#        tc_tot = np.array(tc).T
#    else:
#        if tc_tot.shape == np.array(tc).T.shape:
#            tc_tot += np.array(tc).T
##    if chunk_index == 0:
##        break
#frames_read = chunk_index*chunk_size+sub_ind+1
#print frames_read
#
#tc_tot = tc_tot*1.0/tc_tot[0]

########## Outputs ##########
fig = plt.figure(figsize = (16, 12), dpi = 300)
tc_avg = []
for i, item in enumerate(tc):
    tc_avg.append([])
    base = np.mean(item[0])
    for j in item:
        tc_avg[i].append(np.mean(j)*1.0/base)
    plt.plot(np.arange(len(item))*dt/1000.0, tc_avg[i], label = res_name[i])
plt.xlabel('Time (ns)')
plt.ylabel('Correlation function')
plt.legend()
plt.savefig(outname+'.pdf')

tc_avg = np.array(tc_avg).T
np.savetxt(outname+'.txt', np.column_stack((np.arange(len(tc_avg))*dt/1000.0, tc_avg)), fmt = ['%12.4f' for row in range(1+len(tc_avg[0]))])
