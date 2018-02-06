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

def calc_rlist(xyz, botzone, topzone, bc, zmin, zmax):
    '''Create a list of distances to the functional spot'''
    rlist_bot = [] # create a list to store the distance
    rlist_top = []
    for spot_coord in botzone:
        for item in xyz:
            if item[2] <= zmin:
                rlist_bot.append(calc_dist(item[:2], spot_coord, bc))
    for spot_coord in topzone:
        for item in xyz:
            if item[2] >= zmax:
                rlist_top.append(calc_dist(item[:2], spot_coord, bc))
    return rlist_bot, rlist_top

############ Main #############

# Load trajactory file
traj = md.load('traj.trr', top='topol.pdb')

# Load the mole mass and partial charge from 'para.dat', which is manually made from topology file
para = np.loadtxt('para.dat',dtype=[('atom','S5'),('charge',float),('mass',float)])

# Create dictionaries of partial charge and molar mass for each atom type
pcharge = {item[0] : item[1] for item in para} # Not used in this situation
pmass = {item[0] : item[2] for item in para}

# Select the system
choice = input('Please select your system:\n\
                1. Prestine graphene\n\
                2. Hydrogenated graphene\n\
                3. Epoxy functionalized graphene\n\
                4. Hydroxul functionalized graphene\n')

# Some necessary inputs
frame_start, frame_end = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+':\n')
#zmin, zmax = input('Please enter zmin and zmax (nm):\n')
#binsize = input('Please enter the bin size (nm):\n')
#cutoff = input('Please enter the cutoff (nm):\n')
zmin = 1.3
zmax = 9.42
binsize = 0.01
cutoff = 3
distances = np.arange(0, cutoff, binsize)
ndist_bot_cat =  np.zeros(int(cutoff/binsize), dtype = np.int)
ndist_bot_anion =  np.zeros(int(cutoff/binsize), dtype = np.int)
ndist_top_cat =  np.zeros(int(cutoff/binsize), dtype = np.int)
ndist_top_anion =  np.zeros(int(cutoff/binsize), dtype = np.int)

# count the number of cations
graphene_index = []
cat_index = []
anion_index = []
for res in traj.top.residues:
    if res.name == 'CAT':
        cat_index.append(res.index)
    if res.name == 'BF4':
        anion_index.append(res.index)
    if res.name == 'GPH':
        graphene_index.append(res.index)

# Calculate the index for functional group
if choice != 1:
    fungroup = ['HG','OG','OG'] # atom name of functional group
    ztop = [10.7, 10.7, 9.78]
    botspot = [] # a list of index for functional groups at bottom electrode
    topspot = [] # a list of index for functional groups at top electrode
    selection = traj.top.select('name '+fungroup[choice-2])
    for i in selection:
        if traj.xyz[0, i, 2] > 0.68 and traj.xyz[0, i, 2] < 2:
            botspot.append(i)
        elif traj.xyz[0, i, 2] > ztop[choice-2]-0.682-2 and traj.xyz[0, i, 2] < ztop[choice-2]-0.68:
            topspot.append(i)

# Get the relative indexes of the atoms in the ring to calculate the comparing vector of the cation
ring_index = [] # a list of index for the atoms of the ring
for atom in traj.top.residue(0).atoms:
    if atom.name == 'C2' or\
        atom.name == 'N3' or\
        atom.name == 'C4' or\
        atom.name == 'C5':
        ring_index.append(atom.index)
    elif atom.name == 'N1':
        ring_index.insert(0,atom.index)

for i,frame in enumerate(traj[frame_start-1:frame_end]):
    xyz_com = calc_com( cat_index+anion_index, frame, pmass)
    xyz_cat = xyz_com[cat_index,:]
    xyz_anion = xyz_com[anion_index,:]
    if choice == 1:
        break
    else:
        botzone = frame.xyz[0, botspot, :2]
        topzone = frame.xyz[0, topspot, :2]
        rlist_bot_cat, rlist_top_cat = calc_rlist(xyz_cat, botzone, topzone, traj.unitcell_lengths[0,:2], zmin, zmax)
        rlist_bot_anion, rlist_top_anion = calc_rlist(xyz_anion, botzone, topzone, traj.unitcell_lengths[0,:2], zmin, zmax)
        n_temp, _ = np.histogram(rlist_bot_cat, bins = cutoff/binsize, range = (0, cutoff))
        ndist_bot_cat += n_temp
        n_temp, _ = np.histogram(rlist_top_cat, bins = cutoff/binsize, range = (0, cutoff))
        ndist_top_cat += n_temp
        n_temp, _ = np.histogram(rlist_bot_anion, bins = cutoff/binsize, range = (0, cutoff))
        ndist_bot_anion += n_temp
        n_temp, _ = np.histogram(rlist_top_anion, bins = cutoff/binsize, range = (0, cutoff))
        ndist_top_anion += n_temp
    print 'calculating frame',i+frame_start,'of',traj.n_frames

ndist_bot_cat = ndist_bot_cat.astype(float)
ndist_top_cat = ndist_top_cat.astype(float)
ndist_bot_anion = ndist_bot_anion.astype(float)
ndist_top_anion = ndist_top_anion.astype(float)
frames_read = frame_end-frame_start+1
for i in range(1, len(ndist_bot_cat)):
    area = np.pi*(distances[i]**2-distances[i-1]**2)
    ndist_bot_cat[i] /= area*frames_read*len(botspot)
    ndist_top_cat[i] /= area*frames_read*len(botspot)
    ndist_bot_anion[i] /= area*frames_read*len(botspot)
    ndist_top_anion[i] /= area*frames_read*len(botspot)

########## Output ##########

np.savetxt('DistanceDistribution_cation.txt',np.column_stack((distances, ndist_bot_cat, ndist_bot_anion, ndist_top_cat, ndist_top_anion)),fmt=['%12.4f','%12.4f','%12.4f', '%12.4f', '%12.4f'])
plt.plot(distances, ndist_bot_cat,'r')
plt.plot(distances, ndist_bot_anion,'b')
plt.plot(distances, ndist_top_cat,'r',ls='--')
plt.plot(distances, ndist_top_anion,'b',ls='--')
plt.legend(['bottom-cation','bottom-anion','top-cation','top-anion'])
plt.xlabel('Distance (nm)')
plt.ylabel('Number density ($\#/nm^2$)')
plt.savefig('DistanceDistributionProfile_cation.pdf')

print time.time()-time_start
