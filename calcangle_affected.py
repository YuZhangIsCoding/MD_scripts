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

def calc_cat_com(cat_index,frame,ring_index,pmass):
    '''Calculate coordinates of RTILs using center of mass'''
    j = len(cat_index)
    xyz_com = np.zeros((j,3))
    totalmass = np.zeros(j)
    cog_ring = np.zeros((j,3))
    vhead = np.zeros((j,3))
    for resno in cat_index:
        xyz_pbc = move_pbc(frame, resno)
        for i, atom in enumerate(frame.top.residue(resno).atoms): 
            xyz_com[resno, : ] += pmass[atom.name]*xyz_pbc[i,:]
            totalmass[resno] += pmass[atom.name]
        vhead[resno,:] = xyz_pbc[ring_index[0],:]
        cog_ring[resno,:] = np.sum(xyz_pbc[ring_index,:],axis=0)/5
    xyz_com /= np.array([totalmass]).T # some of the coordinates may extend the boundary limit
    g2h = vhead-cog_ring
    return xyz_com, g2h

def calc_dist(xy,tzone,bc,funrange):
    '''Calculate the smallest projected distance between any cation and the functionalized spot'''
    dist=[]
    for i, item in enumerate(tzone):
        x = abs(xy[0]-item[0])
        y = abs(xy[1]-item[1])
        x = min(x, abs(bc[0]-x))
        y = min(y, abs(bc[1]-y))
        dist.append((x**2+y**2)**0.5)
    if min(dist) < funrange:
        return 1
    else:
        return 0

def calc_range(xyz,g2h,tzone,bc,zmin,zmax,funrange):
    '''Judge if the coordinate is within the effective range of functional group'''
    funindx=[]
    for i, item in enumerate(xyz):
        if item[2] <= zmin or item[2] >= zmax:
            if calc_dist(item[:2],tzone,bc,funrange):
                funindx.append(i)
    return g2h[funindx]

############ Main #############

# Load trajactory file
traj = md.load('test.trr', top='topol.pdb')

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
zmin, zmax = input('Please enter zmin and zmax (nm):\n')
funrange = input('Please enter the effect range of functional spot by the number of bond:\n')
funrange *= 0.139 # Multiply by the bond length of C-C in graphene
binsize = 2
nangle = np.zeros((180/binsize,3),dtype=np.int) # Summary of the frequecy of each bin

# count the number of cations
graphene_index = []
cat_index = []
for res in traj.top.residues:
    if res.name == 'CAT':
        cat_index.append(res.index)
    if res.name == 'GPH':
        graphene_index.append(res.index)

# Calculate the index for functional group
if choice != 1:
    fungroup=['HG','OG','OG'] # atom name of functional group
    fspot = [] # a list of index for functional groups
    selection = traj.top.select('name '+fungroup[choice-2])
    for i in selection:
        if traj.xyz[0, i, 2] > 0.68 and traj.xyz[0, i, 2] < 2:
            fspot.append(i)



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
    [xyz_cat, g2h] = calc_cat_com(cat_index, frame, ring_index, pmass)
    if choice == 1:
        zcount=0
        for indx, z_com in enumerate(xyz_cat[:,2]):
            if z_com > zmin and z_com < zmax:
                zcount += 1
                g2h = np.delete(g2h, indx-zcount,0)
    else:
        tzone = frame.xyz[0, fspot, :2]
        g2h = calc_range(xyz_cat, g2h, tzone, traj.unitcell_lengths[0,:2], zmin,zmax,funrange)
    cosangle=g2h/(np.array([np.sqrt(np.sum(g2h**2,axis=1))]).T)
    angle=np.degrees(np.arccos(cosangle))
    for j, item in enumerate(angle.T):
        n_temp, _ = np.histogram(item, bins=180/binsize, range=(0,180))
        nangle[:,j] += n_temp
    print 'calculating frame',i+frame_start,'of',traj.n_frames

########## Output ##########

degrees = np.arange(0,180,binsize)
np.savetxt('AngleDistribution.txt',np.column_stack((degrees,nangle[:,0],nangle[:,1],nangle[:,2])),fmt=['%12.4f','%12.4f','%12.4f','%12.4f'])
plt.plot(degrees,nangle[:,0],'r')
plt.plot(degrees,nangle[:,1],'g')
plt.plot(degrees,nangle[:,2],'b')
plt.legend(['x','y','z'])
plt.xlabel('${\Theta}$')
plt.ylabel('Frequency')
plt.savefig('AngleDistributionProfile.pdf')

print time.time()-time_start
