#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename: calc_common.py
# Description: This is a common file used for post-processing
# Date: 12-25-2015  Created and added basic functions
#       02-15-2015  Added self part to enable command line input of multiple 
#                   parameters using python module argparse
#                   Combined load_nonbond() and load_para()
#       07-18-2016  For load_para(), fixed the problem for aqueous solution by
#                   changing the atom name to mdtraj defaults
#       09-14-2016  Added new function calc_pangle() to calculate the angle 
#                   distribution based on cos(theta) instead of theta
#       07-05-2017  Added new function load_gro() to load gromacs files
#       01-18-2018  Modified the load_atomic() function, speed up when the
#                   system is large.

import sys, os,  argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

########## Self ###########
## Use argparse module to handle multiple command line input for single option!
parser = argparse.ArgumentParser(description = 'User specified filenames, boundaries, etc.')
parser.add_argument('--aname', dest = 'aname', nargs = '*', help = 'Atom name for the angle vector')
parser.add_argument('--atom', dest = 'atom', type = int, help = 'Atom index, starting from 1')
parser.add_argument('-b', '--bound', dest = 'bound', nargs = '*', type = float, help = 'Boundary of interest')
parser.add_argument('-bs', '--bin_size', dest = 'bin_size', type = float, default = 0.01, help = 'Bin size, default is 0.01')
parser.add_argument('-c', '--cutoff', dest = 'cutoff', type = float, default = 1.1, help = 'Cut-off distance for nonbonded interaction, default is 1.1')
parser.add_argument('-et', '--energytype', dest = 'et', nargs = '*', help = 'Select the type of non-bonded energy of interest')
parser.add_argument('--exb', dest = 'exb', type = float, default = 0.682, help = 'The thickness of electrode')
parser.add_argument('-d', '--direction', dest = 'direction', type = int, choices = range(3), default = 2, help = 'Direction of profile, default in Z direction')
parser.add_argument('--delimiter', default = ' ', help = 'Delimiter for loading data')
parser.add_argument('-dim', '--dimennsion', dest = 'dim', nargs = '*', type = float, help = 'Dimension of the simulation box')
parser.add_argument('-g', '--group', dest = 'group', nargs = '*', help = 'Group of interest')
parser.add_argument('-i', '--input', dest = 'filename', default = 'traj.trr', help = 'Filename for the trajectory, default is traj.trr')
parser.add_argument('-mm', '--molmass', dest = 'mmass', nargs = '*', type = float, help = 'The molar mass of moleculars of interest')
parser.add_argument('--nframe', dest = 'nframe', type = int, default = 0, help = 'The total number of frames')
parser.add_argument('-q', dest = 'q', type = float, default = 2, help = 'The jump length considered mobile (default 2 angstrom)')
parser.add_argument('-rev', '--reverse', dest = 'rev', action = 'store_true', help = 'reverse charged states')
parser.add_argument('--showgroups', dest = 'showgroups', action = 'store_true', help = 'Show the groups of the system')
parser.add_argument('--suffix', dest = 'suffix', type = str, help = 'Suffix after each output file')
parser.add_argument('--scale', dest = 'scale', type = float, nargs = '*', help = 'Select the region you want for the fitting')
parser.add_argument('--scharge', '-sc', dest = 'sc', type = float, default = 0, help = 'surface charge density')
parser.add_argument('--skip', dest = 'skip', type = int, default = 1, help = 'Only calculate only nr-th frame')
parser.add_argument('-t', '--temperature', dest = 'temperature', type = float, default = 298, help = 'Temperature of the system')
parser.add_argument('--test', action = 'store_true', help = 'Test the code or not, default is false')
parser.add_argument('-wg', '--wallgroup', dest = 'wallgroup', nargs = '*', help = 'Group for the electrode')
parser.add_argument('-dt', '--dt', dest = 'dt', type = float, default = 2, help = 'The timestep for each trajectory frame')
parser.add_argument('-sub', '--subplots', dest = 'sub', action = 'store_true', help = 'subplots for each data set')
args = parser.parse_args()
direction = args.direction
bin_size = args.bin_size
dt = args.dt
########## Functions ##########
def load_traj():
    '''load_traj() -> traj_name , traj0, topology

    Load the trajectory file. Filename is specified after -i option
    If not filename specified, will read the default trajectory file traj.trr
    Return the filename, first frame and topology'''
    if args.filename != None:
        traj_name = args.filename
    else:
        traj_name = 'traj.trr'
    if not os.path.isfile(traj_name):
        sys.exit('Exit: not trajectory file named %s found in current directory!' %traj_name)
    try:
        traj0 = md.load_frame(traj_name, 0, top='begin.gro')
    except:
        sys.exit('Fail to load trajectory, check your gro files')
    topology = traj0.top
    return traj_name, traj0, topology

def load_nonbond():
    '''load_nonbond() -> nonbond

    Read the itp files from current directory and search for LJ parameters
    Returns a dictionary of atomtypes and its corresponding LJ parameters
    Incoporated into load_para() ...... 02-15-2016'''
    nonbond = {}
    temp = 0
    for itp_name in os.listdir('./'):
        if itp_name.endswith('.itp'):
            itpfile = open(itp_name, 'r')
            for line in itpfile:
                if '[' in line and 'atomtypes' in line and ']' in line:
                    temp = 1
                    continue
                elif '[' in line:
                    temp = 0
                elif line[0] == ';' or line == '\n':
                    continue
                if temp == 1:
                    entry = line.split()
                    nonbond[entry[0]] = [float(entry[-2]), float(entry[-1])]
    return nonbond

def load_para():
    '''loac_para() -> para

    Read itp files from current directories and search for force field parameters. 
    First obtain a dictionary of atomtypes and their corresponding LJ parameters
    Then obtain a dictionary of atom names and their correponding molar mass, 
    atomic charge and atomtype-based LJ parameters.
    Return a dictionary of atom names and corresponding force field parameters'''
    nonbond = {}
    temp = 0
    for itp_name in os.listdir('./'):
        if itp_name.endswith('.itp'):
            itpfile = open(itp_name, 'r')
            for line in itpfile:
                if '[' in line and 'atomtypes' in line and ']' in line:
                    temp = 1
                    continue
                elif '[' in line:
                    temp = 0
                elif line[0] == ';' or line == '\n':
                    continue
                if temp == 1:
                    entry = line.split()
                    #if entry[0] == 'OW':
                    #    nonbond['O'] = [float(entry[-2]), float(entry[-1])]
                    #elif entry[0][:2] == 'HW':
                    #    nonbond['H'] = [float(entry[-2]), float(entry[-1])]
                    #else:
                    #    nonbond[entry[0]] = [float(entry[-2]), float(entry[-1])]
                    nonbond[entry[0]] = [float(entry[-2]), float(entry[-1])]
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
                    if len(entry) != 0 and entry[0] != ';' and line[0] != '#':
                        if entry[4] == 'OW':
                            #para['O'] = [float(entry[-1]), float(entry[-2]), nonbond[entry[1]]]
                            para['O'] = [15.9994, -0.8476,nonbond[entry[1]]]
                        elif entry[4] == 'HW1':
                            para['H1'] = [1.008, 0.4238, nonbond[entry[1]]]
                        elif entry[4] == 'HW2':
                            para['H2'] = [1.008, 0.4238, nonbond[entry[1]]]
                        else:
                            para[entry[-4]] = [float(entry[-1]), float(entry[-2]), nonbond[entry[1]]]
    return para
def load_atomic(res_targ, topology, para):
    '''load(para, res_targ) -> massind, chargeind
    res_targ:   a list of lists of residue indexes
    topology:   mdtraj topology object
    para:       a dictionary of atom names and their corresponding ff parameters
    
    Creat lists of molar mass and atomic charge according to atom index. If residues
    not selected in res_targ, all atoms will be assign with 0 mass and atomic charge.
    Return a list of molar mass and a list of atomic charge.'''
    massind = np.zeros(topology.n_atoms)
    chargeind = np.zeros(topology.n_atoms)
    for res in sum(res_targ, []):
        for atom in topology.residue(res).atoms:
            massind[atom.index] = para[atom.name][0]
            chargeind[atom.index] = para[atom.name][1]

#    for res_list in res_targ:
#        for res in res_list:
#            for atom in topology.residue(res).atoms:
#                massind[atom.index] = para[atom.name][0]
#                chargeind[atom.index] = para[atom.name][1]

## First version slow when the system is large
#    for atom in topology.atoms:
#        if atom.residue.index in sum(res_targ, []):
#            massind.append(para[atom.name][0]) 
#            chargeind.append(para[atom.name][1])
#        else:
#            massind.append(0)
#            chargeind.append(0)
    return massind, chargeind

def load_bound(traj, direction):
    '''load_bound(traj, direction) -> bound
    traj:       mdtraj trajectory for only one frame
    direction:  int, indicating the direction of the profile

    Boundary of interest can be specified after -b option
    If not specified, use the boundary of the simulation system
    Return a list of 2 elements, lower and upper bounds respectively'''
    if args.bound != None:
        if len(args.bound) == 1 and args.bound[0] == 0:
            bound = [0, np.max(traj[0].xyz[0, :, direction])]
        elif len(args.bound) %2 == 0:
            bound = args.bound
        else:
            print('Boundary entered:', args.bound)
            sys.exit('Exit: boundary not recognized')
    else:
        temp = ['X', 'Y', 'Z']
        print('The system has a length of %f in %s direction' %(np.max(traj[0].xyz[0, :, direction]), temp[direction]))
        bound = input('Please select the boundary, 0 for default:\n')
        if bound == '0':
            bound = [0, np.max(traj[0].xyz[0, :, direction])]
        else:
            bound = [float(i) for i in bound.split()]
    return bound

def creat_atom_list(res_targ, topology):
    '''creat_atom_list(res_targ, topology) -> atom_targ, atom_list
    res_targ:   a list consists of sorted lists of residue index numbers
    topology:   mdtraj topology object
    
    Sort the atom index of selected residues by their atom names.
    Return a list of atom names and a list of lists of corresponding atom indexes'''
    atom_targ = []
    atom_dict = {}
    atom_list = []
    count = 0
    for i in res_targ:
        for item in i:
            for atom in topology.residue(item).atoms:
                if atom.name not in atom_targ:
                    atom_list.append([])
                    atom_targ.append(atom.name)
                    atom_dict[atom.name] = count
                    atom_list[atom_dict[atom.name]].append(atom.index)
                    count += 1
                else:
                    atom_list[atom_dict[atom.name]].append(atom.index)
    return atom_targ, atom_list

def select_groups(traj, direction = 0, bound_sep = [2, 10, 21, 29]):
    '''select_groups(traj) -> res_targ, res_name
    traj:   mdtraj trajectory object

    Select the groups of residues according to their names. If 0 is entered, all of the 
    residues will be selected.
    Return a list of lists of residues selected and a list of residue names selected'''
    count = 0
    res_list = []
    res_dict = {}
    name_list = []
    temp = '0\tALL\n'
    for res in traj.top.residues:
        if res.name not in res_dict:
            res_dict[res.name] = count
            name_list.append(res.name)
            temp += '%d\t%s\n' %(count+1, res.name)
            res_list.append([])
            count += 1
        res_list[res_dict[res.name]].append(res.index)
    if args.group != None:
        group = [int(i) for i in args.group]
    else:
        ## this might cause problems for python 2, where raw_input() should be used instead
        group = input('Please select the groups:\n'+temp)
        group = [int(i) for i in group.split()]
    if len(group) == 1 and group[0] == 0:
        res_targ = res_list
        res_name = name_list
    else:
        res_targ = [res_list[item-1] for item in group]
        res_name = [name_list[item-1] for item in group]
    if args.wallgroup != None:
        wallgroup = [res_list[int(item)-1] for item in args.wallgroup]
        wallname = [name_list[int(item)-1] for item in args.wallgroup]
        wallindex = [[], []]
        for wall_type, wall_list in enumerate(wallgroup):
            for wall_res in wall_list:
                for atom in traj.top.residue(wall_res).atoms:
                    if traj.xyz[0, atom.index, direction] >= bound_sep[0] and traj.xyz[0, atom.index, direction] <= bound_sep[1]:
                        wallindex[0].append(atom.index)
                    elif traj.xyz[0, atom.index, direction] >= bound_sep[2] and traj.xyz[0, atom.index, direction] <= bound_sep[3]:
                        wallindex[1].append(atom.index)
        return res_targ, res_name, wallindex, wallname
    else:
        return res_targ, res_name

def radial_pairs(traj):
    '''radial_pairs(traj) -> res_targ, res_name
    
    Select pairs of groups for the radial distribution. If 0 is entered,
    all of possible pairs will be selected.
    Return a list of pairs of residue selected and a list of names of residue pairs'''
    count = 0
    res_list = []
    res_dict = {}
    name_list = []
    res_targ = []
    res_name = []
    temp = '0\tALL\n'
    for res in traj.top.residues:
        if res.name not in res_dict:
            res_dict[res.name] = count
            name_list.append(res.name)
            temp += '%d\t%s\n' %(count+1, res.name)
            res_list.append([])
            count += 1
        res_list[res_dict[res.name]].append(res.index)
    if args.group != None:
        group = args.group
    else:
        group = input('Please select pairs of groups:\n'+temp).split()
        print(group)
    group = [int(i) for i in group]
    if len(group) == 1 and group[0] == 0:
        for i in range(len(res_list)):
            for j in range(i, len(res_list)):
                res_targ.append([res_list[i], res_list[j]])
                res_name.append([name_list[i], name_list[j]])
    elif len(group) %2 == 0:
        res_targ = [[] for row in range(len(group)/2)]
        res_name = [[] for row in range(len(group)/2)]
        for i, item in enumerate(group):
            res_targ[i/2].append(res_list[item-1])
            res_name[i/2].append(name_list[item-1])
    else:
        print('Groups selected:', group)
        sys.exit('Exit: please enter pairs of groups!')
    return res_targ, res_name

def select_vec(res_targ, topology):
    '''select_vect(res_targ, topology) -> vec_ind
    res_targ:   a list of lists of residue indexes
    topology:   mdtraj topology object

    Select the vectors that can represent the orientations. The index is base on
    the sequence in that molecule, not the real index in the whole system.
    Return the atom index of the selected atoms in that molecule'''
    if args.aname != None:
        aname = args.aname
    else:
        aname = input('Please input the atoms for the angle vector:\n').split()
    vec_ind = []
    for item in res_targ:
        vec_ind.append([])
        for atom in topology.residue(item[0]).atoms:
            if atom.name in aname:
                vec_ind[-1].append(atom.index - topology.residue(item[0]).atom(0).index)
            elif topology.residue(item[0]).name == 'HOH':
                if atom.name in ['H1', 'O']:
                    vec_ind[-1].append(atom.index - topology.residue(item[0]).atom(0).index)
    return vec_ind

def calc_interaction(xyz1, xyz2, nb1, nb2, box, cutoff):
    '''calc_dist(xyz1, xyz2, box) -> vdw, es
    xyz1, xyz2: array of coordinates of an atom
    nb1, nb2:   list of nonbond parameters
    box:        array of box lengths
    cutoff:     cutoff distance for the nonbonded interaction, int

    Return the interaction energy between 2 atoms, van der waals and electrostatic respectively'''
    temp = np.abs(xyz1-xyz2)
    r2 = 0
    for i in range(3):
        if temp[i] < box[i]-temp[i]:
            r2 += temp[i]*temp[i]
        else:
            r2 += (box[i]-temp[i])*(box[i]-temp[i])
    r = np.sqrt(r2)
    es = 0
    vdw = 0
    if r <= cutoff:
        if args.et == None or 'es' in args.et:
            es = nb1[0]*nb2[0]/(4*np.pi*r)*1.74587*1000
        if args.et == None or 'vdw' in args.et:
            if len(nb1[1]) == 2:
                sig6 = (nb1[1][0]+nb2[1][0])**6/2**6
                sig12 = sig6**2
                r6 = r2*r2*r2
                epsilon = np.sqrt(nb1[1][1]*nb2[1][1])
                vdw = 4*epsilon*(sig12/r6/r6-sig6/r6)
            elif len(nb1[1]) == 3:
                sys.exit('Exit: Buckingham potential not available right now!')
    return vdw, es


def calc_nb(res_targ, wallindex, frame, topology, para, direction, cutoff):
    '''
    Return a list of lists of nonbond energy'''
    box = frame.unitcell_lengths[0, :]
    xyz_com = [[] for row in range(len(res_targ))]
    vdw = [[] for row in range(len(res_targ))]
    es = [[] for row in range(len(res_targ))]
    para_me_res = [[] for row in range(len(res_targ))]
    para_nb_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_me_res[res_type].append(para[atom.name][:2])
            para_nb_res[res_type].append(para[atom.name][2])
        para_me_res[res_type] = np.array(para_me_res[res_type])
        para_nb_res[res_type] = np.array(para_nb_res[res_type])
        for res_index in res_list:
            xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
            if min(xyz_org[:, direction]) < cutoff+args.exb:
                vdw_sum = 0
                es_sum = 0
                xyz_com[res_type].append(np.sum(np.multiply(xyz_org[:, direction], para_me_res[res_type][:, 0]))/sum(para_me_res[res_type][:, 0]))
                for atom in topology.residue(res_index).atoms:
                    for ref_type, ref_list in enumerate(wallindex):
                        for ref_index in ref_list:
                            vdw_sub, es_sub = calc_interaction(frame.xyz[0, atom.index, :], frame.xyz[0, ref_index, :], para[atom.name][1:], para[topology.atom(ref_index).name][1:], box, cutoff)
                            vdw_sum += vdw_sub
                            es_sum += es_sub
                vdw[res_type].append(vdw_sum)
                es[res_type].append(es_sum)
    return xyz_com, vdw, es

def calc_nb_slit(res_targ, wallindex, frame, topology, para, direction, cutoff, bound_slit):
    '''
    Return a list of lists of nonbond energy for ions inside slits'''
    box = frame.unitcell_lengths[0, :]
    xyz_com = [[[] for row in range(len(res_targ))] for slit in range(2)]
    vdw = [[[] for row in range(len(res_targ))] for slit in range(2)]
    es = [[[] for row in range(len(res_targ))] for slit in range(2)]
    para_me_res = [[] for row in range(len(res_targ))]
    para_nb_res = [[] for row in range(len(res_targ))]
    bound_slit = [4, 8, 23, 27]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_me_res[res_type].append(para[atom.name][:2])
            para_nb_res[res_type].append(para[atom.name][2])
        para_me_res[res_type] = np.array(para_me_res[res_type])
        para_nb_res[res_type] = np.array(para_nb_res[res_type])
        for res_index in res_list:
            xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
            for slit_ind in range(2):
                if min(xyz_org[:, direction]) <= bound_slit[2*slit_ind+1] and max(xyz_org[:, direction]) >= bound_slit[2*slit_ind]:
                    temp = np.sum(np.multiply(xyz_org.T, para_me_res[res_type][:, 0]), axis = 1)/sum(para_me_res[res_type][:, 0])
                    if temp[direction] >= bound_slit[2*slit_ind] and temp[direction] <= bound_slit[2*slit_ind+1]:
                        xyz_com[slit_ind][res_type].append(temp)
                        vdw_sum = 0
                        es_sum = 0
                        for atom in topology.residue(res_index).atoms:
                            for ref_index in wallindex[slit_ind]:
                                vdw_sub, es_sub = calc_interaction(frame.xyz[0, atom.index, :], frame.xyz[0, ref_index, :], para[atom.name][1:], para[topology.atom(ref_index).name][1:], box, cutoff)
                                vdw_sum += vdw_sub
                                es_sum += es_sub
                        vdw[slit_ind][res_type].append(vdw_sum)
                        es[slit_ind][res_type].append(es_sum)
    return xyz_com, vdw, es

def calc_angle(res_targ, vec_ind, frame, topology, para, direction, bound):
    '''calc_angle(res_targ, vec_ind, frame, topology, para, direction) -> angle
    res_targ:   a list of lists of residue indexes
    vec_ind:    a list of lists of atoms for the angle vector
    frame:      mdtraj single trajectory object
    topology:   mdtraj topology object
    para:       a dictionary of atom names and their corresponding ff parameters

    Return a list of lists of angle vector within the range of selected boundaries.
    '''
    angle = [[[] for row in range(len(res_targ))] for row in range(len(bound)/2)]
    box = frame.unitcell_lengths[0, :]
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name][:2])
        para_res[res_type] = np.array(para_res[res_type])
        for res_index in res_list:
            xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
            for i, item in enumerate(xyz_org.T):
                if np.ptp(item) > box[i]/2:
                    item[item < box[i]/2] += box[i]
                    xyz_org[:,i] = np.array(item).T
            xyz_temp = np.sum(np.multiply(xyz_org.T, para_res[res_type][:,0]).T, axis=0)/sum(para_res[res_type][:,0])
            for j in range(len(bound)/2):
                if xyz_temp[direction] >= bound[2*j] and xyz_temp[direction] <= bound[2*j+1]:
                    if len(vec_ind[res_type]) == 2:
                        v_nor = xyz_org[vec_ind[res_type][1], :] - xyz_org[vec_ind[res_type][0], :]
                    elif len(vec_ind[res_type]) == 3:
                        v1 = xyz_org[vec_ind[res_type][1], :] - xyz_org[vec_ind[res_type][0], :]
                        v2 = xyz_org[vec_ind[res_type][2], :] - xyz_org[vec_ind[res_type][0], :]
                        v_nor = np.cross(v1, v2)
                    temp = v_nor[direction]/np.sqrt(sum(v_nor**2))
                    if temp <= -1:
                        print('Warning product less than -1:', temp)
                        angle[j][res_type].append(180)
                    elif temp >= 1:
                        print('Warning product higher than 1:', temp)
                        angle[j][res_type].append(0)
                    else:
                        angle[j][res_type].append(np.degrees(np.arccos(temp)))
    return angle

def calc_cosangle(res_targ, vec_ind, frame, topology, para, direction, bound):
    '''calc_angle(res_targ, vec_ind, frame, topology, para, direction) -> cosa
    res_targ:   a list of lists of residue indexes
    vec_ind:    a list of lists of atoms for the angle vector
    frame:      mdtraj single trajectory object
    topology:   mdtraj topology object
    para:       a dictionary of atom names and their corresponding ff parameters

    Return a list of lists of p2cos within the range of selected boundaries.
    '''
    cosa = [[[] for row in range(len(res_targ))] for row in range(len(bound)/2)]
    box = frame.unitcell_lengths[0, :]
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name][:2])
        para_res[res_type] = np.array(para_res[res_type])
        for res_index in res_list:
            xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
            for i, item in enumerate(xyz_org.T):
                if np.ptp(item) > box[i]/2:
                    item[item < box[i]/2] += box[i]
                    xyz_org[:,i] = np.array(item).T
            xyz_temp = np.sum(np.multiply(xyz_org.T, para_res[res_type][:,0]).T, axis=0)/sum(para_res[res_type][:,0])
            for j in range(len(bound)/2):
                if xyz_temp[direction] >= bound[2*j] and xyz_temp[direction] <= bound[2*j+1]:
                    if len(vec_ind[res_type]) == 2:
                        v_nor = xyz_org[vec_ind[res_type][1], :] - xyz_org[vec_ind[res_type][0], :]
                    elif len(vec_ind[res_type]) == 3:
                        v1 = xyz_org[vec_ind[res_type][1], :] - xyz_org[vec_ind[res_type][0], :]
                        v2 = xyz_org[vec_ind[res_type][2], :] - xyz_org[vec_ind[res_type][0], :]
                        v_nor = np.cross(v1, v2)
                    cosa[j][res_type].append(v_nor[direction]/np.sqrt(sum(v_nor**2)))
    return cosa

def calc_pangle(res_targ, vec_ind, frame, topology, para, direction, bound):
    '''calc_angle(res_targ, vec_ind, frame, topology, para, direction) -> cosa
    res_targ:   a list of lists of residue indexes
    vec_ind:    a list of lists of atoms for the angle vector
    frame:      mdtraj single trajectory object
    topology:   mdtraj topology object
    para:       a dictionary of atom names and their corresponding ff parameters

    Return a list of lists of p2cos within the range of selected boundaries.
    '''
    cosa = [[[] for row in range(len(res_targ))] for row in range(len(bound)/2)]
    box = frame.unitcell_lengths[0, :]
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name][:2])
        para_res[res_type] = np.array(para_res[res_type])
        for res_index in res_list:
            xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
            for i, item in enumerate(xyz_org.T):
                if np.ptp(item) > box[i]/2:
                    item[item < box[i]/2] += box[i]
                    xyz_org[:,i] = np.array(item).T
            xyz_temp = np.sum(np.multiply(xyz_org.T, para_res[res_type][:,0]).T, axis=0)/sum(para_res[res_type][:,0])
            for j in range(len(bound)/2):
                if xyz_temp[direction] >= bound[2*j] and xyz_temp[direction] <= bound[2*j+1]:
                    if len(vec_ind[res_type]) == 2:
                        v_nor = xyz_org[vec_ind[res_type][1], :] - xyz_org[vec_ind[res_type][0], :]
                    elif len(vec_ind[res_type]) == 3:
                        v1 = xyz_org[vec_ind[res_type][1], :] - xyz_org[vec_ind[res_type][0], :]
                        v2 = xyz_org[vec_ind[res_type][2], :] - xyz_org[vec_ind[res_type][0], :]
                        v_nor = np.cross(v1, v2)
                    cosa[j][res_type].append(1.5*v_nor[direction]**2/sum(v_nor**2)-0.5)
    return cosa

def calc_xyz_direction(res_targ, frame, topology, para, direction):
    '''calc_xyz_direction(res_targ, frame, topology, direction) -> xyz_dir
    res_targ:   a list of lists of residue indexes
    frame:      mdtraj single trajectory object
    topology:   mdtraj topology object
    para:       a dictionary of atom names and their corresponding ff parameters

    Calculate the coordinates of molecules in selected direction based on their center of mass. 
    Assume molecules did not move accross the boundary in selected direction.
    Return a list of arrays of coordinate in one direction, sorted by their residue types.'''
    xyz_dir = [[] for row in range(len(res_targ))]
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        if res_list == []:
            continue
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name][:2])
        para_res[res_type] = np.array(para_res[res_type])
        xyz_dir[res_type] = np.zeros(len(res_list))
        #for res_index in res_list:
        for res_seq, res_index in enumerate(res_list):
            xyz_org = np.array([frame.xyz[0, atom.index, direction] for atom in topology.residue(res_index).atoms])
            xyz_dir[res_type][res_seq] = np.sum(np.multiply(xyz_org, para_res[res_type][:, 0]))/sum(para_res[res_type][:, 0])
    return xyz_dir

def calc_xyz_com(res_targ, frame, topology, para):
    '''calc_xyz_com(res_targ, frame, topology, para) -> xyz_com
    res_targ:   a list of lists of residue indexes
    frame:      mdtraj single trajectory object
    topology:   mdtraj topology object
    para:       a dictionary of atom names and their corresponding ff parameters
    
    Calculate the xyz coordinates of molecules in all 3 directions based on 
    their center of mass. All COM coordinates are inside the system box.
    Return a list of arrays of xyz, sorted by their residue types.'''
    box = frame.unitcell_lengths[0, :]
    xyz_com = [[] for row in range(len(res_targ))]
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        if res_list == []:
            continue
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name][:2])
        para_res[res_type] = np.array(para_res[res_type])
        xyz_com[res_type] = np.zeros((len(res_list), 3))
        for res_seq, res_index in enumerate(res_list):
            xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
            for i, item in enumerate(xyz_org.T):
                if np.ptp(item) > box[i]/2:
                    item[item < box[i]/2] += box[i]
                    xyz_org[:,i] = np.array(item).T
            xyz_com[res_type][res_seq, :] = np.sum(np.multiply(xyz_org.T, para_res[res_type][:,0]).T, axis=0)/sum(para_res[res_type][:,0])
            for item in range(3):
                if xyz_com[res_type][res_seq, item] > box[item]:
                    xyz_com[res_type][res_seq, item] -= box[item]
    return xyz_com

def calc_xyz_cyl(res_targ, frame, topology, para, cyl):
    '''calc_xyz_cyl(res_targ, frame, topology, para, cyl) -> xyz_sel
    res_targ:   a list of lists of residue indexes
    frame:      mdtraj single trajectory object
    topology:   mdtraj topology object
    para:       a dictionary of atom names and their corresponding ff parameters
    cyl:        a list parameters for cylinder, including coordinates for the center and radius

    Calculate the xyz coordinates of molecules in all 3 directions based on 
    their center of mass. All COM coordinates are inside the cylinder specified by cyl.
    This is useful to get the atomic coordinates in selected cylinderical space, for example to get
    the number density profile near AFM tips.
    Return a list of arrays of xyz inside cylinder, sorted by their residue types'''
    box = frame.unitcell_lengths[0, :]
    xyz_com = [[] for row in range(len(res_targ))]
    xyz_sel = [[] for row in range(len(res_targ))]
    #totalmass = np.zeros(len(res_targ))
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name][:2])
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
        for item in xyz_com[res_type]:
            if (item[0]-cyl[0])**2+(item[1]-cyl[1])**2 <= cyl[2]**2:
                xyz_sel[res_type].append(item)
        xyz_sel[res_type] = np.array(xyz_sel[res_type])
    return xyz_sel

def calc_xyz_com_con(res_targ, frame, topology, para, direction, bound):
    '''calc_xyz_com_con(res_targ, frame, topology, para, direction, bound) -> xyz_com
    res_targ:   a list of lists of residue indexes
    frame:      mdtraj single trajectory object
    topology:   mdtraj topolpgy object
    para:       a dictionary of atom names and their corresponding ff parameters
    direction:  the direction of boundary restrictions
    bound:      the selected boundary limits

    Calculate the xyz coordinates of molecule in selected region.
    Return a list of arrays of xyz in selected region, sorted by their residue types'''
    box = frame.unitcell_lengths[0, :]
    xyz_com = [[] for row in range(len(res_targ))]
    para_res = [[] for row in range(len(res_targ))]
    for res_type, res_list in enumerate(res_targ):
        for atom in topology.residue(res_list[0]).atoms:
            para_res[res_type].append(para[atom.name][:2])
        para_res[res_type] = np.array(para_res[res_type])
        #xyz_com[res_type] = np.zeros((len(res_list), 3))
        for res_index in res_list:
            xyz_org = np.array([frame.xyz[0, atom.index, :] for atom in topology.residue(res_index).atoms])
            # prejudge if it is a suitable ion 
            if np.max(xyz_org[:, direction]) > bound[0] and np.min(xyz_org[:, direction]) < bound[1]:
                for i, item in enumerate(xyz_org.T):
                    if np.ptp(item) > box[i]/2:
                        item[item < box[i]/2] += box[i]
                        xyz_org[:,i] = np.array(item).T
                temp = np.sum(np.multiply(xyz_org.T, para_res[res_type][:,0]).T, axis=0)/sum(para_res[res_type][:, 0])
                if temp[direction] >= bound[0] and temp[direction] < bound[1]:
                    for dim in range(3):
                        if temp[dim] > box[dim]:
                            temp[dim] -= box[dim]
                    xyz_com[res_type].append(temp)
        xyz_com[res_type] = np.array(xyz_com[res_type])
    return xyz_com

def plot_axis(fig, out_name, fsize):
    '''plot_axis(fig, outname)
    fig:        matplotlib figure object
    outname:    name for output, string

    Adjust axis and some common setup and save a fig'''
    for axis in fig.axes:
        axis.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
        for direction in ['top', 'bottom', 'left', 'right']:
            axis.spines[direction].set_linewidth(2)
    plt.tight_layout()
    plt.savefig(out_name+'.pdf')

def out_nd(distances, nd, res_name, out_name = None, fsize = 36, lwidth = 4, exb = 0.682):
    '''out_nd(nd) -> NumberDensity_suffix.txt, NumberDensity.pdf
    distances:  array of distances
    nd:         array of number densities

    Output a txt file and plot of the number density'''
    if not out_name:
        if args.suffix != None:
            out_name = 'NumberDensity_'+args.suffix
        else:
            out_name = 'NumberDensity'
    fmt = ['%12.4f' for row in range(1+nd.shape[1])]
    np.savetxt(out_name+'.txt', np.column_stack((distances, nd)), fmt=fmt)
    fig = plt.figure(figsize = (16, 12), dpi = 1000)
    for i in range(len(nd[0])):
        plt.plot(distances-exb, nd[:, i], linewidth = lwidth, label = res_name[i-1])
    plt.xlim((0, max(distances)-2*exb))
    plt.legend(loc = 'upper center', fontsize = fsize, frameon = False)
    plt.xlabel('Distance (nm)', fontsize = fsize)
    plt.ylabel('Number density $\mathregular{\#/nm^3}$', fontsize = fsize)
    plot_axis(fig, out_name, fsize)

def out_density(distances, mass_density, charge_density, out_name = None, fsize = 36, lwidth = 4, exb = 0.682):
    '''Output a txt file and plot of the mass density and charge density'''
    
    if not out_name:
        if args.suffix != None:
            out_name = 'Density_'+args.suffix
        else:
            out_name = 'Density'
    np.savetxt(out_name+'.txt', np.column_stack((distances, mass_density, charge_density)), fmt=['%12.4f','%12.4f','%12.4f'])
    fig, (ax1, ax2) = plt.subplots(2, figsize = (16, 12), dpi = 1000, sharex = True)
    ax1.plot(distances-exb, mass_density, linewidth = lwidth)
    ax1.set_title('Mass density profile', fontsize = fsize)
    ax1.set_xlim([0, max(distances)-2*exb])
    ax1.set_ylabel('Mass density $\mathregular{(kg/m^3)}$', fontsize = fsize)
    ax2.plot(distances-exb, charge_density, linewidth = lwidth)
    ax2.set_title('Charge density profile', fontsize = fsize)
    ax2.set_xlabel('Distance (nm)', fontsize = fsize)
    ax2.set_ylabel('Charge density $\mathregular{(e/nm^3)}$', fontsize = fsize)
    plot_axis(fig, out_name, fsize)

def load_gro(filename):
    '''Read gro file and return atom info and coordinate

    Return type: a list of string and numpy array.
    '''
    myfile = open(filename, 'r')
    coord = []
    for i, line in enumerate(myfile):
        temp = line.split()
        if i > 1:
            coord.append([line[:20], [float(item) for item in temp[-3:]]])
    return coord 

def load_xvg(filename):
    '''load xvg file into numpy array'''

    myfile = open(filename, 'r')
    data = [] 
    for line in myfile:
        if line[0] != '#' and line[0] != '@':
            data.append([float(i) for i in line.split()])
    myfile.close()
    return np.array(data)

def sort_data(data):
    '''Sort the data for charged system from negative charge to positive charge.
    The initial data should starting from pzc to systems with higher surface charges'''
    hl = len(data[0])/2
    charges = np.transpose([np.arange(-len(data)+1, len(data))])
    temp = np.concatenate((data[1:, hl:][::-1], [(data[0, :hl]+data[0, hl:])/2], data[1:, :hl]))
    data = np.column_stack((charges, temp))
    return data

########## Self_2 ############
if args.showgroups == True:
    traj_name, traj0, topology = load_traj()
    temp = '0\tALL'
    print('Groups in the system:\n',temp)
    count = 1
    for res in topology.residues:
        if res.name != temp:
            print(count,'\t',res.name)
            temp = res.name
            count += 1
    sys.exit()
