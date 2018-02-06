#!/Users/yuzhang/anaconda/bin/python
import pdb
import mdtraj as md
import numpy as np
import time
import matplotlib.pyplot as plt
time_start = time.time()

########## Functions ##########

def move_pbc(frame,resno):
    '''Make the molecule whole again according to PBC'''
    xyz_org = np.array([ frame.xyz[0,atom.index,:] for atom in frame.top.residue(resno).atoms])
    for i, item in enumerate(xyz_org.T[:2,:]):
        if np.ptp(item) > traj.unitcell_lengths[0, i]/2:
            item[item < traj.unitcell_lengths[0, i]/2] += traj.unitcell_lengths[0, i]
            xyz_org[:,i] = np.array(item).T
    return xyz_org

def calc_com(resid_index, frame, pmass):
    ''' General function to calculate the coordinate of a molecule based on Center of Mass'''
    j = len(resid_index)
    xyz_com = np.zeros((j,3))
    totalmass = np.zeros(j)
    for resno in resid_index:
        xyz_pbc = move_pbc(frame, resno)
        for i, atom in enumerate(frame.top.residue(resno).atoms):
            xyz_com[resno, : ] += pmass[atom.name]*xyz_pbc[i,:]
            totalmass[resno] += pmass[atom.name]
    xyz_com /= np.array([totalmass]).T # some of the coordinates may extend the boundary limit
    return xyz_com

def calc_dist(xy,spot_coord,bc):
    '''Calculate the smallest projected distance between any cation and the functionalized spot'''
    x = abs(xy[0]-spot_coord[0])
    y = abs(xy[1]-spot_coord[1])
    x = min(x, abs(bc[0]-x))
    y = min(y, abs(bc[1]-y))
    return (x**2+y**2)**0.5

def calc_flist(xyz, bot_list, top_list, botzone, topzone, frange, bc):
    '''Create a list of distances to the functional spot'''
    bot_list_new = []
    top_list_new = []
    for spot_coord in botzone:
        for item in bot_list:
            if calc_dist(xyz[item], spot_coord, bc) <= frange:
                bot_list_new.append(item)
    for spot_coord in topzone:
        for item in top_list:
            if calc_dist(xyz[item], spot_coord, bc) <= frange:
                top_list_new.append(item)
    return bot_list_new, top_list_new

############ Main #############

# Load trajactory file
traj = md.load('traj.trr', top='topol.pdb')

# Load the mole mass and partial charge from 'para.dat', which is manually made from topology file
para = np.loadtxt('para.dat',dtype=[('atom','S5'),('charge',float),('mass',float)])

# Create dictionaries of partial charge and molar mass for each atom type
pcharge = {item[0] : item[1] for item in para} # Not used in this situation
pmass = {item[0] : item[2] for item in para}

# Select the system
#choice = input('Please select your system:\n\
#                1. Prestine graphene\n\
#                2. Hydrogenated graphene\n\
#                3. Epoxy functionalized graphene\n\
#                4. Hydroxul functionalized graphene\n')

choice = 2

# Some necessary inputs
#frame_start, frame_end = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+':\n')

frame_start = 1
frame_end = 751

#cutoff = input('Please input the cutoff distance away from the electrode surface (nm) :\n')
#binsize = input('Please input the binsize (nm) :\n')
#frange = input('Please input the effective range of functional group')
frange = 1
bondlength = 0.139
frange *= bondlength
cutoff = 1
binsize = 0.01
distances = np.arange(0, cutoff, binsize)

# count the number of cations
cat_index = []
anion_index = []
for res in traj.top.residues:
    if res.name == 'CAT':
        cat_index.append(res.index)
    elif res.name == 'BF4':
        anion_index.append(res.index)

atom_all = [] # Create a list of all atom names
for i in traj.top.residue(cat_index[0]).atoms:
    atom_all.append(i.name)
for i in traj.top.residue(anion_index[0]).atoms:
    atom_all.append(i.name)

# Calculate the index for functional group
ztop = max(traj.xyz[0,:,2]) # the z value of top electrode surface
if choice != 1:
    fungroup = ['HG','OG','OG'] # atom name of functional group
    botspot = [] # a list of index for functional groups at bottom electrode
    topspot = [] # a list of index for functional groups at top electrode
    selection = traj.top.select('name '+fungroup[choice-2])
    for i in selection:
        if traj.xyz[0, i, 2] > 0.68 and traj.xyz[0, i, 2] < 2:
            botspot.append(i)
        elif traj.xyz[0, i, 2] > ztop-0.682-2 and traj.xyz[0, i, 2] < ztop-0.68:
            topspot.append(i)

ndist_bot = np.zeros((len(atom_all),len(distances)), dtype = np.int)
ndist_top = np.zeros((len(atom_all),len(distances)), dtype = np.int)

for i,frame in enumerate(traj[frame_start-1:frame_end]):
    xyz_com = calc_com( cat_index+anion_index, frame, pmass)
    xyz_cat = xyz_com[cat_index,:]
    xyz_anion = xyz_com[anion_index,:]
    bot_list = []
    top_list = []
    for res_ind, item in enumerate(xyz_com):
        if item[2] < 0.682+cutoff:
            bot_list.append(res_ind)
        elif item[2] > ztop-0.682-cutoff:
            top_list.append(res_ind)
    if choice != 1:
        botzone = frame.xyz[0, botspot, :2]
        topzone = frame.xyz[0, topspot, :2]
        bot_list, top_list = calc_flist(xyz_com, bot_list, top_list, botzone, topzone, frange, traj.unitcell_lengths[0,:2])
    rlist_bot = [[] for row in range(len(atom_all))]
    rlist_top = [[] for row in range(len(atom_all))]
    for atom_ind, atom_name in enumerate(atom_all):
        for res_ind in bot_list:
            for atom in frame.top.residue(res_ind).atoms:
                if atom.name == atom_name:
                    rlist_bot[atom_ind].append(frame.xyz[0,atom.index,2]-0.682)
        for res_ind in top_list:
            for atom in frame.top.residue(res_ind).atoms:
                if atom.name == atom_name:
                    rlist_top[atom_ind].append(ztop-0.682-frame.xyz[0,atom.index,2])
    for atom_ind, nlist in enumerate(rlist_bot):
        n_temp, _ = np.histogram(nlist, bins = len(distances), range = (0, cutoff))
        ndist_bot[atom_ind] += n_temp
    for atom_ind, nlist in enumerate(rlist_top):
        n_temp, _ = np.histogram(nlist, bins = len(distances), range = (0, cutoff))
        ndist_top[atom_ind] += n_temp
    print 'calculating frame',i+frame_start,'of',traj.n_frames

frames_read = frame_end-frame_start+1
ndist_bot = ndist_bot.astype(float)
ndist_top = ndist_top.astype(float)
ndist_bot /= frames_read*np.prod(traj.unitcell_lengths[0][:2])*binsize
ndist_top /= frames_read*np.prod(traj.unitcell_lengths[0][:2])*binsize
########## Output ##########

np.savetxt('Atom_z_bot.txt',np.concatenate((np.array([distances]).T,ndist_bot.T), axis = 1),fmt=['%12.4f']*(len(ndist_bot)+1))
np.savetxt('Atom_z_top.txt',np.concatenate((np.array([distances]).T,ndist_top.T), axis = 1),fmt=['%12.4f']*(len(ndist_bot)+1))
#plt.plot(distances, ndist_bot_cat,'r')
#plt.plot(distances, ndist_bot_anion,'b')
#plt.plot(distances, ndist_top_cat,'r',ls='--')
#plt.plot(distances, ndist_top_anion,'b',ls='--')
#plt.legend(['bottom-cation','bottom-anion','top-cation','top-anion'])
#plt.xlabel('Distance (nm)')
#plt.ylabel('Number density ($\#/nm^2$)')
#plt.savefig('DistanceDistributionProfile_cation.pdf')

print time.time()-time_start
