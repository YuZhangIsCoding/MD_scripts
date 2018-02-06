#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_density_planar_iter.py
# Description: This is a python script to get the density profiles of a channel simulation
# Dates:    10-08-2015 Created
#           02-05-2015 Modified to split the trajectory by chunks to save memory

import sys, os, pdb, time
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import calc_common as comm


########## Main ###########
traj_name, traj0, topology = comm.load_traj()
direction = comm.direction

atom_ind = comm.args.atom-1

chunk_size = 100
mycoord = []
pbc = 0
for chunk_index, traj in enumerate(md.iterload(traj_name, chunk = chunk_size, top = 'topol.pdb')):
    for sub_ind, frame in enumerate(traj):
        box = frame.unitcell_lengths[0, 0]
        if mycoord == []:
            pre = frame.xyz[0, atom_ind, 0]
        else:
            if frame.xyz[0, atom_ind, 0]-pre < -box/2:
                pbc += 1
            elif frame.xyz[0, atom_ind, 0]-pre > box/2:
                pbc -= 1
        mycoord.append(frame.xyz[0, atom_ind, 0]+pbc*box)
        pre = frame.xyz[0, atom_ind, 0]
        if sub_ind%10 == 9:
            print 'Reading chunk', chunk_index+1, 'and frame',chunk_index*chunk_size+sub_ind+1
fig, (ax1, ax2) = plt.subplots(2)

vel = []
time = np.arange(len(mycoord))*comm.args.dt
print 'time step (ps):', comm.args.dt
l = 10
fit_order = 1
print 'Fitting polynomial with order of', fit_order
for i in range(len(mycoord)):
    if i >= l and i < len(mycoord)-l:
        z = np.polyfit(time[i-l:i+l+1], mycoord[i-l:i+l+1], fit_order)
        df = np.poly1d(z[:-1]*np.arange(len(z)-1, 0, -1))
        vel.append(df(time[i])*1000)

ax1.plot(time, mycoord)
ax2.plot(time[l:l+len(vel)], vel)
ax2.set_xlim([time[0], time[-1]])
ax2.set_xlabel('Time (ps)')
ax1.set_ylabel('Displacement (nm)')
ax2.set_ylabel('Speed (nm/ns)')
plt.savefig('speed_atom.pdf')
