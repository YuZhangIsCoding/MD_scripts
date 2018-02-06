#!/usr/bin/python
# Filename: calc_template.py
# Date: 10-08-2015 Created
# Description: This is a template file for trajectory analysis

import sys
import os
import pdb
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

def select_groups():
    '''Select the groups of interest'''
    count = 0
    res_list = []
    res_dict = {}
    #res_dict = {} #necessary if want to recall the resname
    temp = '0\tALL\n'
    for res in traj.top.residues:
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
    return res_targ

def calc_xyz_com(res_targ, frame, topology):
    '''Calculate the coordinated of a group according to center of mass'''
    box = frame.unitcell_lengths[0, :]
    xyz_com = [[] for row in range(len(res_targ))]
    #totalmass = np.zeros(len(res_targ))
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
    pdb.set_trace()
    return xyz_com

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

# Define the frames wanted
frames = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+', 0 for default:\n')
if frames == 0:
    frames = [1, traj.n_frames]
frames_read = frames[1]-frames[0]+1

res_targ = select_groups()
pre_time = time.time()
for i, frame in enumerate(traj[frames[0]-1: frames[1]]):
    xyz_com = calc_xyz_com(res_targ, frame, topology)
    if i%10 == 9:
        rem_time = (time.time()-pre_time)*(frames_read-i+1)/(i+1)
        if rem_time > 300:
            rem_time = '%.2f min' %(rem_time/60)
        else:
            rem_time = '%.2f s' %rem_time
        print 'Calculating frame', i+frames[0], 'of', traj.n_frames,'Estimated to finish in', rem_time


########## Outputs ##########
