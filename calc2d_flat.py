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
    for i in range(len(xyz_com)):
        for j in range(3):
            if xyz_com[i, j ] > frame.unitcell_lengths[0, j]:
                xyz_com[i, j] -= frame.unitcell_lengths[0, j]
            elif xyz_com[i, j] < 0:
                xyz_com[i, j] += frame.unitcell_lengths[0, j]
    return xyz_com

def calc_dist(xy,spot_coord,bc):
    '''Calculate the smallest projected distance between any cation and the functionalized spot'''
    x = abs(xy[0]-spot_coord[0])
    y = abs(xy[1]-spot_coord[1])
    x = min(x, abs(bc[0]-x))
    y = min(y, abs(bc[1]-y))
    return (x**2+y**2)**0.5

def grid_match(xy, mapping, binsize):
    '''Distribute ion to grid coordinates'''
    mapping[int(xy[0]/binsize),int(xy[1]/binsize)] += 1
    return mapping

############ Main #############

# Load trajactory file
traj = md.load('traj.trr', top='topol.pdb')

# Load the mole mass and partial charge from 'para.dat', which is manually made from topology file
para = np.loadtxt('para.dat',dtype=[('atom','S5'),('charge',float),('mass',float)])

# Create dictionaries of partial charge and molar mass for each atom type
pcharge = {item[0] : item[1] for item in para} # Not used in this situation
pmass = {item[0] : item[2] for item in para}

## Select the system
#choice = input('Please select your system:\n\
#                1. Prestine graphene\n\
#                2. Hydrogenated graphene\n\
#                3. Epoxy functionalized graphene\n\
#                4. Hydroxul functionalized graphene\n')
choice = 3

# Some necessary inputs
#frame_start, frame_end = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+':\n')
#zmin, zmax = input('Please enter zmin and zmax (nm):\n')
#binsize = input('Please enter the bin size (nm):\n')
#cutoff = input('Please enter the cutoff (nm):\n')
frame_start = 1
frame_end = traj.n_frames
binsize = 0.01
zmin = 1.3
zmax = 9.32
dimensions = [int(traj.unitcell_lengths[0, 0]/binsize+1), int(traj.unitcell_lengths[0, 1]/binsize)+1]
#distances = np.arange(0, cutoff, binsize)
#ndist_bot_cat =  np.zeros(int(cutoff/binsize), dtype = np.int)
#ndist_bot_anion =  np.zeros(int(cutoff/binsize), dtype = np.int)
#ndist_top_cat =  np.zeros(int(cutoff/binsize), dtype = np.int)
#ndist_top_anion =  np.zeros(int(cutoff/binsize), dtype = np.int)

# count the number of cations
graphene_index = []
cat_index = []
anion_index = []
for res in traj.top.residues:
    if res.name == 'CAT':
        cat_index.append(res.index)
    elif res.name == 'BF4':
        anion_index.append(res.index)
    elif res.name == 'GPH':
        graphene_index.append(res.index)

#Calculate the index for functional group
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

cat_bot = np.zeros(dimensions,dtype = np.int)
cat_top = np.zeros(dimensions,dtype = np.int)
anion_bot = np.zeros(dimensions,dtype = np.int)
anion_top = np.zeros(dimensions,dtype = np.int)
fun_bot = np.zeros(dimensions,dtype = np.int)
fun_top = np.zeros(dimensions,dtype = np.int)

for i,frame in enumerate(traj[frame_start-1:frame_end]):
    xyz_com = calc_com( cat_index+anion_index, frame, pmass)
    xyz_cat = xyz_com[cat_index,:]
    xyz_anion = xyz_com[anion_index,:]
    for item in xyz_cat:
        if item[2] <= zmin:
            cat_bot = grid_match(item[:2], cat_bot, binsize)
        elif item[2] >= zmax:
            cat_top = grid_match(item[:2], cat_top, binsize)
    for item in xyz_anion:
        if item[2] <= zmin:
            anion_bot = grid_match(item[:2], anion_bot, binsize)
        elif item[2] >= zmax:
            anion_top = grid_match(item[:2], anion_top, binsize)
    if choice != 1:
        botzone = frame.xyz[0, botspot, :2]
        topzone = frame.xyz[0, topspot, :2]
        for item in botzone:
            fun_bot = grid_match(item, fun_bot, binsize)
        for item in topzone:
            fun_top = grid_match(item, fun_top, binsize)
    print 'calculating frame',i+frame_start,'of',traj.n_frames

########## Output ##########

#np.savetxt('DistanceDistribution_cation.txt',np.column_stack((distances, ndist_bot_cat, ndist_bot_anion, ndist_top_cat, ndist_top_anion)),fmt=['%12.4f','%12.4f','%12.4f', '%12.4f', '%12.4f'])
#plt.plot(distances, ndist_bot_cat,'r')
outname = ['cat_bot', 'cat_top', 'anion_bot', 'anion_top', 'fun_bot', 'fun_top']
for i, item in enumerate([cat_bot, cat_top, anion_bot, anion_top, fun_bot, fun_top]):
    outfile = open('2d_'+outname[i]+'.txt', 'w')
    for row in item:
        temp = [str(a) for a in row]
        outfile.write('%s\n' %(" ".join(temp)))
    outfile.close()
#np.savetxt('2dmapping_bot.txt', mapping_bot)
#np.savetxt('2dmapping_top.txt'. mapping_top)
#plt.savefig('DistanceDistributionProfile_cation.pdf')

print time.time()-time_start
