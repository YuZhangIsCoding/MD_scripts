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

def calc_dist(xyz,spot_coord,bc):
    '''Calculate the smallest projected distance between any cation and the functionalized spot'''
    x = abs(xyz[0]-spot_coord[0])
    y = abs(xyz[1]-spot_coord[1])
    z = xyz[2]-spot_coord[2]
    x = min(x, abs(bc[0]-x))
    y = min(y, abs(bc[1]-y))
    return (x**2+y**2+z**2)**0.5

def calc_buckingham(coeff, dist, cutoff):
    '''Calculate the buckingham interaction between 2 atoms'''
    if dist > cutoff:
        E = 0
    else:
        E = coeff[0]*np.exp(-coeff[1]*dist)-coeff[2]/(dist**6)
    return E

def calc_lj69(coeff, dist, cutoff):
    '''Calculate the Lennard-Lones(6-9) potential of 2 atoms'''
    if dist > cutoff:
        E = 0
    else:
        E = -coeff[0]/dist**6+coeff[1]/dist**9
    return E

def calc_elect(coeff, dist, cutoff):
    '''Calculate the electrostatic potential between 2 atoms'''
    if dist > cutoff:
        E = 0
    else:
        E = 40.46*coeff[0]*coeff[1]/dist # 40.46 is the constant
    return E


def calc_potential(all_bot, spot, frame, cutoff, pcharge):
    ''' Calculate the interaction potential between molecules and electrode'''
    potential = np.zeros(len(all_bot)) 
    for i, resind in enumerate(all_bot):
        for atom in frame.top.residue(resind).atoms:
            for spotind in spot:
                spotname = frame.top.atom(spotind).name
                coeff = [pcharge[spotname],pcharge[atom.name]]
                dist = calc_dist(frame.xyz[0,atom.index,:], frame.xyz[0, spotind, :], frame.unitcell_lengths[0,:2])
                potential[i] += calc_elect(coeff,dist, cutoff)           
    return potential

############ Main #############

# Load trajactory file
traj = md.load('test.trr', top='topol.pdb')

# Load the mole mass and partial charge from 'para.dat', which is manually made from topology file
para = np.loadtxt('para.dat',dtype=[('atom','S5'),('charge',float),('mass',float)])

# Create dictionaries of partial charge and molar mass for each atom type
pcharge = {item[0] : item[1] for item in para} # Not used in this situation
pmass = {item[0] : item[2] for item in para}

# Some necessary inputs
#frame_start, frame_end = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+':\n')

frame_start = 1500
frame_end = 1600

#cutoff = input('Please input the cutoff distance away from the electrode surface (nm) :\n')
#binsize = input('Please input the binsize (nm) :\n')
#frange = input('Please input the effective range of functional group')
cutoff = 1.1
binsize = 0.01
distances = np.arange(0, cutoff, binsize)

# count the number of cations
cat_index = []
bf4_index = []
for res in traj.top.residues:
    if res.name == 'CAT':
        cat_index.append(res.index)
    elif res.name == 'BF4':
        bf4_index.append(res.index)

# Specify the inner electrodes
ztop = max(traj.xyz[0,:,2]) # the z value of top electrode surface
selection = traj.top.select('resname GPH and not name CG')
botspot = []
topspot = []
for i in selection:
    if traj.xyz[0, i, 2] > 0.4 and traj.xyz[0, i ,2] < 1:
        botspot.append(i)
    elif traj.xyz[0, i ,2]  > ztop-0.682-0.3 and traj.xyz[0, i, 2] < ztop-0.682+0.3:
        topspot.append(i)

ndist_bot = np.zeros((len(distances),3), dtype = np.int)
ndist_top = np.zeros((len(distances),3), dtype = np.int)
npotential = np.zeros((len(distances),3))

for i,frame in enumerate(traj[frame_start-1:frame_end]):
    xyz_com = calc_com( cat_index+bf4_index, frame, pmass)
    #xyz_cat = xyz_com[cat_index,:]
    #xyz_bf4 = xyz_com[bf4_index,:]
    #xyz_tf2n = xyz_com[tf2n_index,:]
    all_bot = [[] for row in range(2)]
    all_top = [[] for row in range(2)]
    for res_ind, item in enumerate(xyz_com):
        if item[2] < 0.682+cutoff:
            if res_ind <= max(cat_index):
                all_bot[0].append(res_ind)
            elif res_ind <= max(bf4_index):
                all_bot[1].append(res_ind)
        elif item[2] > ztop-0.682-cutoff:
            if res_ind <= max(cat_index):
                all_top[0].append(res_ind)
            elif res_ind <= max(bf4_index):
                all_top[1].append(res_ind)
    potential = [[] for row in range(2)]
    for j, item in enumerate(all_bot):
        potential[j] = calc_potential( item, botspot, frame, cutoff, pcharge)
        xyz_bot = xyz_com[item]-0.682
        n_temp, _ = np.histogram(xyz_bot[:,2], bins = len(distances), range = (0, cutoff))
        ndist_bot[:,j] += n_temp
        n_temp, _ = np.histogram(xyz_bot[:,2], weights = potential[j], bins = len(distances), range = (0, cutoff))
        npotential[:,j] += n_temp
    print 'calculating frame',i+frame_start,'of',traj.n_frames,
    print 'time', time.time()-time_start
for i in range(2):
    for j, item in enumerate(ndist_bot[:,i]):
        if item != 0:
            npotential[j,i] /= item
########## Output ##########

np.savetxt('elec_number.txt',np.column_stack((distances, ndist_bot[:,0], ndist_bot[:,1])), fmt=['%12.4f','%12.4f','%12.4f'])
np.savetxt('electrostat.txt', np.column_stack((distances, npotential[:,0], npotential[:,1])), fmt=['%12.4f','%12.4f','%12.4f'])
plt.plot(distances, npotential[:,0])
plt.plot(distances, npotential[:,1])
#plt.plot(distances, ndist_bot_cat,'r')
#plt.plot(distances, ndist_bot_anion,'b')
#plt.plot(distances, ndist_top_cat,'r',ls='--')
#plt.plot(distances, ndist_top_anion,'b',ls='--')
#plt.legend(['bottom-cation','bottom-anion','top-cation','top-anion'])
#plt.xlabel('Distance (nm)')
#plt.ylabel('Number density ($\#/nm^2$)')
plt.savefig('electrostat.pdf')

print time.time()-time_start
