#! /usr/bin/python
#  Filename : backup.py

import pdb
import mdtraj as md
import numpy as np
import time
import matplotlib.pyplot as plt
import sys

time_start = time.time()

########### Functions ###########

def calc_z_com(rtils_index,frame):
    '''Calculate the coordinate in z direction according to center of mass'''
    z_com = np.zeros(len(rtils_index))
    totalmass = np.zeros(len(rtils_index))
    for resno in rtils_index:
        for atom in frame.top.residue(resno).atoms:
            z_com[resno] += pmass[atom.name]*frame.xyz[0, atom.index, 2]
            totalmass[resno] += pmass[atom.name]
            #print atom.name, pmass[atom.name], frame.xyz[0, atom.index, 2], z_com[resno], totalmass[resno]
    z_com /= totalmass
    return z_com

#################### Main ####################

# Load trajactory file
if '-i' in sys.argv:
    temp = sys.argv.index('-i')+1
    traj = md.load(sys.argv[temp], top='topol.pdb')
else:
    traj = md.load('traj.trr', top='topol.pdb')

# Read parameter file for molar mass and partial charge of each atom
para = np.loadtxt('para.dat',dtype=[('atom','S5'),('charge',float),('mass',float)])

# Create dictionary of parameters for each atom type
pcharge = {item[0] : item[1] for item in para}
pmass = {item[0] : item[2] for item in para}

# Inputs and initialization
frame_start, frame_end = input('Please choose the starting frame and ending frame from 1 to '+str(traj.n_frames)+':\n')
bin_size = input('Please input the bin size (nm)\n')
bound = np.max(traj[0].xyz[0,:,2])
n_bins = int(np.ceil(bound/bin_size))

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
for i, frame in enumerate(traj[frame_start-1:frame_end]):
    z_org = frame.xyz[0, :, 2]
    z_com = calc_z_com(rtils_index, frame)
    z_com_cat = z_com[:len(rtils_index)/2]
    z_com_anion = z_com[len(rtils_index)/2:]
    temp, _ = np.histogram(z_org, weights = massind, bins = n_bins, range = (0.0, bound))
    mass_density += temp
    temp, _ = np.histogram(z_org, weights = chargeind, bins = n_bins, range = (0.0, bound))
    charge_density += temp
    temp, _ = np.histogram(z_com_cat, bins=n_bins, range=(0.0, bound))
    ncat += temp
    temp, _ = np.histogram(z_com_anion, bins=n_bins, range=(0.0, bound))
    nanion += temp
    print 'calculating frame',i+frame_start,'of',traj.n_frames

# Take average over all frames
frames_read = frame_end-frame_start+1
bin_volume = np.prod(traj.unitcell_lengths[0][:2])*bin_size

mass_density /= bin_volume*frames_read/10*6.022
charge_density /= bin_volume*frames_read
ncat = ncat/bin_volume/frames_read
nanion = nanion/bin_volume/frames_read

########## Output ##########

distances = np.arange(0, bound, bin_size)
np.savetxt('Density.txt', np.column_stack((distances, mass_density, charge_density)), fmt=['%12.4f','%12.4f','%12.4f'])
np.savetxt('NumberDensity.txt',np.column_stack((distances,ncat,nanion)),fmt=['%12.4f','%12.4f','%12.4f'])

f,axarr=plt.subplots(2,sharex=True)
axarr[0].plot(distances, mass_density)
axarr[0].set_title('mass density profile')
axarr[0].set_xlim([0,bound])
axarr[0].set_ylim([0,5000])
axarr[0].set_ylabel('mass density ($kg/m^3$)')
axarr[1].plot(distances,charge_density)
axarr[1].set_title('charge density profile')
axarr[1].set_xlabel('z (nm)')
axarr[1].set_ylabel('charge density ($e/nm^3$)')
plt.savefig('DensityProfile.pdf')

fig2 = plt.figure()
plt.plot(distances,ncat,'r')
plt.plot(distances,nanion,'b')
plt.legend(['Cation','Anion'])
plt.xlabel('z (nm)')
plt.ylabel('Number density')
fig2.savefig('NumberDensityProfile.pdf')

print time.time()-time_start
