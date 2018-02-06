#!/Users/yuzhang/anaconda/bin/python
# Filename: calc_template.py
# Date: 10-08-2015 Created
# Description: This is a template file for trajectory analysis

import sys
import os
import mdtraj as md
import numpy as np
import time
import matplotlib.pyplot as plt

time_start = time.time()

########## Functions ##########

def remove_pbc(frame, res_index, topology):
    '''Remove the periodic boundary and make the molecule whole again'''
    box = frame.unitcell_lengths[0,:]
    xyz_org = np.array([ frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
    for i, item in enumerate(xyz_org.T):
        if np.ptp(item) > box[i]:
            item[item < box[i]/2] += box[i]
    return xyz_org

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
    #choice = input('Please select the groups:\n'+temp)
    choice = (2,3)
    if choice == 0:
        res_targ = res_list
    elif type(choice) is tuple:
        res_targ = [res_list[item-1] for item in choice]
    else:
        res_targ = [res_list[choice-1]]
    return res_targ

def calc_xyz_com(res_targ, frame, vec_index, topology):
    '''Calculate the coordinated of a group according to center of mass'''
    v_nor = [[], []]
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
                if xyz_com[res_type][res_index-res_list[0], item] > box[item]:
                    xyz_com[res_type][res_index-res_list[0], item] -= box[item]
            if res_type == 0:
                v1 = xyz_org[vec_index[0][1], :] - xyz_org[vec_index[0][0], :]
                v2 = xyz_org[vec_index[0][2], :] - xyz_org[vec_index[0][0], :]
                v_nor[0].append(np.cross(v1, v2))
            else:
                v_nor[1].append(xyz_org[vec_index[1][1], :] - xyz_org[vec_index[1][0], :])
    return xyz_com, v_nor

########## Main ###########

# Load trajectory
if '-i' in sys.argv:
    temp = sys.argv.index('-i')+1
    traj_name = sys.argv[temp]
else:
    traj_name = 'traj.trr'

traj = md.load_frame(traj_name, 0,top='topol.pdb')

para = load_para()
topology = traj.top

# Define the frames wanted
#frames = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+', 0 for default:\n')
#frames = 0
#if frames == 0:
#    frames = [1, traj.n_frames]
#frames_read = frames[1]-frames[0]+1

#b_dir = input('Please select the direction for boundary conditions:\n\
#0   x\n\
#1   y\n\
#2   z\n')
#p_dir = input('Please select the direction for number distribution profile:\n\
#0   x\n\
#1   y\n\
#2   z\n')
b_dir = 0
p_dir = 2

def calc_angle_emim(res_targ, frame, topology, vec_index, direction, bound):
    angle = []
    for i, res_index in enumerate(res_targ):
        xyz_pbc = remove_pbc(frame, res_index, topology)
        for j in range(len(bound)/2):
            if max(xyz_pbc[:,0]) > bound[2*j] and min(xyz_pbc[:,0]) <  bound[2*j+1]:
                v1 = xyz_pbc[vec_index[1], :] - xyz_pbc[vec_index[0], :]
                v2 = xyz_pbc[vec_index[2], :] - xyz_pbc[vec_index[0], :]
                v_nor = np.cross(v1, v2)
                angle.append(np.degrees(np.arccos(v_nor[direction]/np.sqrt(sum(v_nor**2)))))
                continue
    return angle

def calc_angle_tf2n(res_targ, frame, topology, vec_index, direction, bound):
    angle = []
    for i, res_index in enumerate(res_targ):
        xyz_pbc = remove_pbc(frame, res_index, topology)
        for j in range(len(bound)/2):
            if max(xyz_pbc[:,0]) > bound[2*j] and min(xyz_pbc[:,0]) <  bound[2*j+1]:
                v_nor = xyz_pbc[vec_index[1], :] - xyz_pbc[vec_index[0], :]
                angle.append(np.degrees(np.arccos(v_nor[direction]/np.sqrt(sum(v_nor**2)))))
                continue
    return angle



res_targ = select_groups()
vec_index = [[],[]]
name_type = ['OMIm','Tf2N']
for atom in topology.residue(res_targ[0][0]).atoms:
    if atom.name == 'C3' or atom.name == 'C5' or atom.name == 'C6':
        vec_index[0].append(atom.index-topology.residue(res_targ[0][0]).atom(0).index)
for atom in topology.residue(res_targ[1][0]).atoms:
    if atom.name == 'C1T' or atom.name == 'C5T':
        vec_index[1].append(atom.index-topology.residue(res_targ[1][0]).atom(0).index)

#bound = [1.75, 2.55]
#bound = [1.75, 4.35]
bound = [1.75, 2.85]
#bound = [1.75, 2.55]
#bound = [0 , 11]
#bound_lim = [3, 9, 19+3, 19+9]
bound_lim = [3, 7, 17+3, 17+7]
bin_size = 2
distances = np.arange(bound[0], bound[1], 0.005)-bound[0]
degrees = np.arange(0, 180, bin_size)
n_angle = [np.zeros((len(degrees),2)), np.zeros((len(degrees),2))]
nd = [np.zeros((len(distances),2)), np.zeros((len(distances),2))]
pre_time = time.time()
chunk_size = 100
for chunk_index, traj in enumerate(md.iterload(traj_name,chunk = chunk_size, top = 'topol.pdb')):
    for i, frame in enumerate(traj):
        xyz_com, v_nor = calc_xyz_com(res_targ, frame, vec_index, topology)
        for j in range(len(bound_lim)/2):
            for res_type, xyz in enumerate(xyz_com):
                xyz_fin = []
                angle = []
                for index, item in enumerate(xyz):
                    if item[b_dir] > bound_lim[2*j] and item[b_dir] < bound_lim[2*j+1]:
                        xyz_fin.append(item[p_dir])
                        angle.append(np.degrees(np.arccos(v_nor[res_type][index][p_dir]/np.sqrt(sum(v_nor[res_type][index]**2)))))
                n_temp, _ = np.histogram(angle, bins = 180/ bin_size, range = (0, 180))
                n_angle[j][:, res_type] += n_temp
                n_temp, _ = np.histogram(xyz_fin, bins = len(distances), range = (bound[0], bound[1]))
                nd[j][:, res_type]  += n_temp
        print 'Reading',chunk_index*chunk_size+i+1
#    if chunk_index == 0:
#        break

########## Outputs ##########
frames_read = chunk_index*chunk_size+i+1
print frames_read
n_angle = np.array(n_angle)
nd = np.array(nd)
bin_v = traj.unitcell_lengths[0][1]*(bound_lim[1]-bound_lim[0])*0.005
ele = ['POS', 'NEG']
for i in range(2):
    nd[i] = nd[i]/bin_v/frames_read
    np.savetxt('NDinSlit_%s.txt' %ele[i], np.column_stack((distances,nd[i][:,0], nd[i][:,1])),fmt=['%12.4f','%12.4f', '%12.4f'])
    np.savetxt('Angle_%s.txt' %ele[i], np.column_stack((degrees, n_angle[i][:,0], n_angle[i][:,1])), fmt=['%12d','%12d', '%12d'])
    plt.figure()
    plt.plot(distances, nd[i][:,0], 'red', label = 'Emim_%s' %ele[i])
    plt.plot(distances, nd[i][:,1], 'blue', label = 'Tf2N_%s' %ele[i])
    plt.xlabel('Distance (nm)')
    plt.ylabel('Number density (#/nm^3)')
    plt.legend()
    plt.savefig('NDinSlit_%s.pdf' %ele[i])
    plt.figure()
    plt.plot(degrees, n_angle[i][:,0], label = 'OMIm_%s' %ele[i])
    plt.plot(degrees, n_angle[i][:,1], label = 'Tf2N_%s' %ele[i])
    plt.xlabel('${\Theta}$')
    plt.ylabel('Probability density')
    plt.legend()
    plt.savefig('Angle_%s.pdf' %ele[i])
