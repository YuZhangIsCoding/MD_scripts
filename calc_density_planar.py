#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_template.py
# Description:  This is a python script that calculate the density profile in a channel simulation
#               It will reads the whole trajectory as a whole chunk, which is memory demanding
#               But it also estimates the ending time
# Date: 10-08-2015 Created
#       02-15-2015 Stop maintaining

import sys
import os
import mdtraj as md
import numpy as np
import time
import matplotlib.pyplot as plt

time_start = time.time()

########## Functions ##########

def load_para():
    '''Read necessary parameters from itp files'''
    para = {}
    temp = 0
    for itp_name in os.listdir('./'):
        if itp_name.endswith('.itp'):
            itpfile = open(itp_name, 'r')
            for line in itpfile:
                if '[' in line and ']' in line and 'atoms' in line:
                    temp = 1
                    continue
                elif '[' in line:
                    temp = 0
                elif line[0] == ';' or line == '\n':
                    continue
                if temp == 1:
                    entry = line.split()
                    para[entry[-4]] = [float(entry[-1]), float(entry[-2])]
    return para

def select_groups(topology, para):
    '''Select the groups of interest'''
    count = 0
    res_list = []
    res_dict = {}
    #res_dict = {} #necessary if want to recall the resname
    temp = '0\tALL\n'
    for res in topology.residues:
        if res.name not in res_dict:
            res_dict[res.name] = count
            temp += '%d\t%s\n' %(count+1, res.name)
            res_list.append([])
            count += 1
        res_list[res_dict[res.name]].append(res.index)
    choice = input('Please select the groups:\n'+temp)
    if choice == 0:
        res_targ = res_list
    elif type(choice) is tuple:
        res_targ = [res_list[item-1] for item in choice]
    else:
        res_targ = [res_list[choice-1]]
    massind = []
    chargeind = []
    for i, atom in enumerate(topology.atoms):
        if atom.residue.index in sum(res_targ, []):
            massind.append(para[atom.name][0])
            chargeind.append(para[atom.name][1])
        else:
            massind.append(0)
            chargeind.append(0)
    return res_targ, massind, chargeind

def calc_xyz_direction(res_targ, frame, topology, direction):
    ''' Calculate teh coordinates of a groups according to center of mass in just one direction'''
    xyz_dir = [[] for row in range(len(res_targ))]
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name])
        para_res[res_type] = np.array(para_res[res_type])
        xyz_dir[res_type] = np.zeros(len(res_list))
        for res_index in res_list:
            xyz_org = np.array([frame.xyz[0, atom.index, direction] for atom in topology.residue(res_index).atoms])
            xyz_dir[res_type][res_index-res_list[0]] = np.sum(np.multiply(xyz_org, para_res[res_type][:, 0]))/sum(para_res[res_type][:, 0])
    return xyz_dir

def calc_xyz_com(res_targ, frame, topology):
    '''Calculate the coordinate of a group according to center of mass'''
    box = frame.unitcell_lengths[0, :]
    xyz_com = [[] for row in range(len(res_targ))]
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name])
        para_res[res_type] = np.array(para_res[res_type])
        xyz_com[res_type] = np.zeros((len(res_list), 3))
        for res_index in res_list:
            xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
            for i, item in enumerate(xyz_org.T):
                if np.ptp(item) > box[i]/2:
                    item[item < box[i]/2] += box[i]
                    xyz_org[:,i] = np.array(item).T
            xyz_com[res_type][res_index-res_list[0], :] = np.sum(np.multiply(xyz_org.T, para_res[res_type][:,0]).T, axis=0)/sum(para_res[res_type][:,0])
            for item in range(3):
                if xyz_com[res_type][res_index-res_list[0], item] > box[i]:
                    xyz_com[res_type][res_index-res_list[0], item] -= box[i]
    return xyz_com

def sort_density_com():
    pre_time = time.time()
    for i, frame in enumerate(traj[frames[0]-1: frames[1]]):
        xyz_com = calc_xyz_com(res_targ, frame, topology)
        for j, item in enumerate(xyz_com):
            temp, _ = np.histogram(item[:, direction], bins = np.append(distances, max(distances)+bin_size))
            numd[:, j] += temp
        if i%10 == 9:
            rem_time = (time.time()-pre_time)*(frames_read-i+1)/(i+1)
            if rem_time > 300:
                rem_time = '%.2f min' %(rem_time/60)
            else:
                rem_time = '%.2f s' %rem_time
            print 'Calculating frame', i+frames[0], 'of', traj.n_frames,'Estimated to finish in', rem_time

########## Main ###########

# Load trajectory
if '-i' in sys.argv:
    temp = sys.argv.index('-i')+1
    traj_name = sys.argv[temp]
else:
    traj_name = 'traj.trr'

traj = md.load(traj_name, top='topol.pdb')
topology = traj.top
para = load_para()
#massind = np.array([para[atom.name][0] for atom in topology.atoms])
#chargeind = np.array([para[atom.name][1] for atom in topology.atoms])

# Select the direction
direction = input('Please input the direction of distribution:\n\
0   x\n\
1   y\n\
2   z\n')

# binning method
bin_size = input('Please input the bin size (nm):\n')
bound = input('Please select the boundary from %f to %f, 0 for default:\n' %(np.min(traj[0].xyz[0,:, direction]), np.max(traj[0].xyz[0,:, direction])))
if bound == 0:
    bound = [np.min(traj[0].xyz[0,:, direction]), np.max(traj[0].xyz[0,:, direction])]
distances = np.arange(bound[0], bound[1], bin_size)

# Define the frames wanted
frames = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+', 0 for default:\n')
if frames == 0:
    frames = [1, traj.n_frames]
frames_read = frames[1]-frames[0]+1

res_targ, massind, chargeind= select_groups(topology, para)

numd = np.zeros((len(distances), len(res_targ)))
mass_density = np.zeros(len(distances))
charge_density = np.zeros(len(distances))
pre_time = time.time()
for i, frame in enumerate(traj[frames[0]-1: frames[1]]):
    temp, _ = np.histogram(frame.xyz[0, :, direction], weights = massind, bins = np.append(distances, max(distances)+bin_size))
    mass_density += temp
    temp, _ = np.histogram(frame.xyz[0, :, direction], weights = chargeind, bins = np.append(distances, max(distances)+bin_size))
    charge_density += temp
    xyz_dir = calc_xyz_direction(res_targ, frame, topology, direction)
    for j, item in enumerate(xyz_dir):
        temp, _ = np.histogram(item, bins = np.append(distances, max(distances)+bin_size))
        numd[:, j] += temp
    if i%10 == 9:
        rem_time = (time.time()-pre_time)*(frames_read-i+1)/(i+1)
        if rem_time > 300:
            rem_time = '%.2f min' %(rem_time/60)
        else:
            rem_time = '%.2f s' %rem_time
        print 'Calculating frame', i+frames[0], 'of', traj.n_frames,'Estimated to finish in', rem_time

mass_density /= bin_size*np.prod(traj.unitcell_lengths[0])/traj.unitcell_lengths[0][direction]*frames_read/10*6.022
charge_density /= bin_size*np.prod(traj.unitcell_lengths[0])/traj.unitcell_lengths[0][direction]*frames_read
numd /= bin_size*np.prod(traj.unitcell_lengths[0])/traj.unitcell_lengths[0][direction]*frames_read

########## Outputs ##########
fmt = ['%12.4f' for row in range(1+numd.shape[1])]
np.savetxt('NumberDensity.txt', np.column_stack((distances,numd)), fmt=fmt)
np.savetxt('Density.txt', np.column_stack((distances, mass_density, charge_density)), fmt=['%12.4f','%12.4f','%12.4f'])

f,axarr=plt.subplots(2,sharex=True)
axarr[0].plot(distances, mass_density)
axarr[0].set_title('mass density profile')
axarr[0].set_xlim(bound)
axarr[0].set_ylim([0,5000])
axarr[0].set_ylabel('mass density ($kg/m^3$)')
axarr[1].plot(distances,charge_density)
axarr[1].set_title('charge density profile')
axarr[1].set_xlabel('z (nm)')
axarr[1].set_ylabel('charge density ($e/nm^3$)')
plt.savefig('Density.pdf')

fig = plt.figure()
plt.plot(distances, numd)
plt.xlabel('Separation (nm)')
plt.ylabel('Number density #/nm^3')
fig.savefig('NumberDensity.pdf')
