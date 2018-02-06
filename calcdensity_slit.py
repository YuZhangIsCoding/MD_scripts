#! /usr/bin/python
#  Filename : backup.py

import sys
import mdtraj as md
import numpy as np
import time
import matplotlib.pyplot as plt

time_start = time.time()

########### Functions ###########

def move_pbc(frame,resno):
    '''Make the molecule whole again according to PBC'''
    xyz_org = np.array([ frame.xyz[0,atom.index,:] for atom in frame.top.residue(resno).atoms])
    for i, item in enumerate(xyz_org.T[:2,:]):
        if np.ptp(item) > traj.unitcell_lengths[0, i]/2:
            item[item < traj.unitcell_lengths[0, i]/2] += traj.unitcell_lengths[0, i]
            xyz_org[:,i] = np.array(item).T
    return xyz_org



def calc_xyz_com(rtils_index,frame, bound_lim):
    '''Calculate the coordinate in z direction according to center of mass'''
    xyz_com = np.zeros((len(rtils_index),3))
    totalmass = np.zeros((len(rtils_index),1))
    for resno in rtils_index:
        xyz_pbc = move_pbc(frame, resno)
        for i, atom in enumerate(frame.top.residue(resno).atoms):
            for j in range(len(bound_lim)/2):
                if xyz_pbc[i,2]-bound_lim[2*j] > -2 and xyz_pbc[i, 2]-bound_lim[2*j+1] < 2:
                    xyz_com[resno-rtils_index[0], :] += pmass[atom.name]*xyz_pbc[i, :]
                    totalmass[resno-rtils_index[0]] += pmass[atom.name]
                    continue
        if totalmass[resno - rtils_index[0]] == 0:
            totalmass[resno-rtils_index[0]] = 1
    xyz_com /= totalmass
    return xyz_com

#################### Main ####################

# Load trajactory file
if '-i' in sys.argv:
    temp = sys.argv.index('-i')+1
    traj_name = sys.argv[temp]
else:
    traj_name = 'traj.trr'

traj = md.load(traj_name, top='topol.pdb')


# Read parameter file for molar mass and partial charge of each atom
para = np.loadtxt('para.dat',dtype=[('atom','S5'),('charge',float),('mass',float)])

# Create dictionary of parameters for each atom type
pcharge = {item[0] : item[1] for item in para}
pmass = {item[0] : item[2] for item in para}

# Inputs and initialization
frames = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+', 0 for default:\n')
if frames == 0:
    frames = [1, traj.n_frames]
frames_read = frames[1]-frames[0]+1

bin_size = input('Please input the bin size (nm)\n')
direction = input('Please input the direction of distribution:\n\
0   x\n\
1   y\n\
2   z\n')
direction_lim = input('Please input the direction that restrictions are applied:\n\
0   x\n\
1   y\n\
2   z\n')
bound = [1.75, 2.86]
bound_lim = [2.1, 4.2, 13.3+2.1, 13.3+4.2]
n_bins = int(np.ceil((bound[1]-bound[0])/bin_size))

# Create arrays to store mass density and charge density
mass_density = np.zeros(n_bins)
charge_density = np.zeros(n_bins)

# Lists of mass and charge for each atom
massind = np.array([pmass[atom.name] for atom in traj.top.atoms])
chargeind = np.array([pcharge[atom.name] for atom in traj.top.atoms])

# Create arrays to store number density of cation and anion
ncat = np.zeros(n_bins,dtype=np.int)
nanion = np.zeros(n_bins,dtype=np.int)

# Get the total number of cation and anion
rtils_index = [] # a list of residue index for RTILs
for res in traj.top.residues:
    if res.name == 'CAT' or res.name == 'Tf2N':
        rtils_index.append(res.index)

# Use histogram to get distribution over z direction
pre_time = time.time()
for i, frame in enumerate(traj[frames[0]-1:frames[1]]):
    xyz_com = calc_xyz_com(rtils_index, frame, bound_lim)
    xyz_com_cat = xyz_com[:len(rtils_index)/2, :]
    xyz_com_anion = xyz_com[len(rtils_index)/2:, :]
    xcount = 0
    targ_cat = []
    targ_anion = []
    for item in xyz_com_cat:
        if item[direction_lim] > bound_lim[0] and item[direction_lim] < bound_lim[1]:
            targ_cat.append(item)
    for item in xyz_com_anion:
        if item[direction_lim] > bound_lim[0] and item[direction_lim] < bound_lim[1]:
            targ_anion.append(item)
    targ_cat = np.array(targ_cat)
    targ_anion = np.array(targ_anion)
    temp, _ = np.histogram(targ_cat[:,direction], bins=n_bins, range=(bound[0], bound[1]))
    ncat += temp
    temp, _ = np.histogram(targ_anion[:,direction], bins=n_bins, range=(bound[0], bound[1]))
    nanion += temp
    if i%10 == 9:
        rem_time = (time.time()-pre_time)*(frames_read-i+1)/(i+1)
        if rem_time > 300:
            rem_time = '%.2f min' %(rem_time/60)
        else:
            rem_time = '%.2f s' %rem_time
        print 'Calculating frame', i+frames[0], 'of', traj.n_frames,'Estimated to finish in', rem_time

# Take average over all frames
bin_volume = traj.unitcell_lengths[0][1]*bin_size*(bound_lim[1]-bound_lim[0])

ncat = ncat/bin_volume/frames_read
nanion = nanion/bin_volume/frames_read

########## Output ##########

distances = np.arange(bound[0], bound[1], bin_size)-1.75
np.savetxt('NDinSlit.txt',np.column_stack((distances,ncat,nanion)),fmt=['%12.4f','%12.4f','%12.4f'])

fig2 = plt.figure()
plt.plot(distances,ncat,'r')
plt.plot(distances,nanion,'b')
plt.legend(['Cation','Anion'])
plt.xlabel('z (nm)')
plt.ylabel('Number density #/nm^3')
fig2.savefig('NDinSlit.pdf')

print time.time()-time_start
