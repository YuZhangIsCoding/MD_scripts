#!/Users/yuzhang/anaconda3/bin/python
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import calc_common as comm
from numpy import linalg as LA

if comm.args.suffix:
    outname = 'amsd_'+comm.args.suffix
else:
    outname = 'amsd'

res_file = open('mol_remain.txt', 'r')
res_targ = []
for line in res_file:
    if line != '\n':
        res_targ.append([int(i) for i in line.split()])

traj_name, traj0, topology = comm.load_traj()
para = comm.load_para()

massind, chargeind = comm.load_atomic(res_targ, topology, para)
dt = comm.dt
nframe = comm.args.nframe

if comm.args.center:
    center = comm.args.center
else:
    center = (4.48300466, 4.48300138, 4.48199593)

#chunk_size = min(100, np.ceil(nframe/10).astype(int))
chunk_size = min(100, nframe//10)
total_chunks = np.ceil(nframe/chunk_size).astype(int)
period = [[0 for i in j] for j in res_targ]
r_ref = []
ad = [[[] for row in range(nframe-nframe//2+1)] for i in range(len(res_targ))]

print('Calculating msd')
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'begin.gro')):
    for sub_ind, frame in enumerate(traj):
        frame_ind = chunk_index*chunk_size+sub_ind
        r_com = comm.calc_xyz_com(res_targ, frame, topology, para)
        ## normalize the vector
        for i, item in enumerate(r_com):
            norm = LA.norm(item-center, axis = 1).reshape((-1, 1))
            r_com[i] = (item-center)/norm
        p_start = 0
        p_end = nframe//2
        if frame_ind < nframe//2:
            r_ref.append(r_com)
            p_end = frame_ind
        elif nframe-nframe//2-frame_ind < 0:
            p_start = frame_ind-nframe+nframe//2
            r_ref[p_start-1] = None # free up some space
        if not frame_ind:
            #r_pre = r_com # prepare for periodicity check
            continue
        period.append(period[-1][:])
        for res_ind, r in enumerate(r_com):
            ##
            ## block to check the periodic condition
            ## leave it for now because the molecules did not travel so far
            ##
            for loop_ind in range(p_start, p_end):
                theta = np.arccos(np.sum(r*r_ref[loop_ind][res_ind], axis = 1))
                ad[res_ind][frame_ind-loop_ind].append(np.mean(theta**2))
        #r_pre = r_com
        if comm.args.test and comm.args.test == frame_ind:
            break
    else:
        print('\r|'+u"\u2588"*(chunk_index+1)+' '*(total_chunks-chunk_index-1)+'| %2d%%' %((chunk_index+1)*100/total_chunks), end = '\r')
        continue
    break

amsd = []
for i, item in enumerate(ad):
    ad[i][0] = 0
    amsd.append([])
    for j in item:
        amsd[-1].append(np.mean(j))

data = np.column_stack((np.arange(nframe-nframe/2+1)*dt, np.array(amsd).T))
np.savetxt(outname+'.txt', data, fmt = ['%12.6f' for i in range(len(data[0]))])

for i in range(1, len(data[0])):
    plt.plot(data[:, 0], data[:, i])
plt.xlabel('ps')
plt.ylabel('$\mathregular{rad^2}$')
plt.savefig(outname+'.pdf')
