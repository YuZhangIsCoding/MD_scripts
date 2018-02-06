#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Description:      Uses mdtraj.compute_distances() and mdtraj.compute_angles()
#                   to compute distance and angles. And then use numpy.where to
#                   quickly locate the atom pairs that greater than the cutoff
#                   angle and within cutoff distance.
import mdtraj as md
import numpy as np
import time
import argparse

outfile = open('hbond_index.txt', 'w')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'filename', help = 'input filename')
args = parser.parse_args()

if args.filename:
    filename = args.filename
else:
    filename = input('Please input filename:\n')

start = time.time()
temp = np.array([0])
## This loads the whole trajectory just once into the memory.
## May want to modify to use md.iterload to load a chunk of traj each time.
traj = md.load(filename, top = 'topol.pdb')
print(len(traj))
topology = traj.top

cutoff = 0.4
acri = 150*np.pi/180
donor_name = ['C3', 'C5', 'C6']
H_name = ['H3', 'H5', 'H6']

donor = []
for i in range(3):
    donor.append([topology.select("name %s" %donor_name[i]), topology.select("name %s" %H_name[i])])
acceptor = topology.select("name =~ 'F.*'")

apairs = []
for k in range(3):
    for i in range(len(donor[k][0])):
        for j in acceptor:
            apairs.append([donor[k][0][i], donor[k][1][i], j])
apairs = np.array(apairs)
indexes = np.array([np.arange(len(apairs)) for _ in range(len(traj))])

dists = md.compute_distances(traj, apairs[:, [0, 2]])
ha_dists = md.compute_distances(traj, apairs[:, [1, 2]])
angles = md.compute_angles(traj, apairs)
# Using this way we could also obtain the info for the distance distribution
#remaining = angles[np.where(dists <= cutoff)]
#final = remaining[np.where(remaining > acri)]
locs = np.where(angles > acri)
remaining = dists[locs]
ha = ha_dists[locs]
indexes = indexes[locs]
locs = np.where(remaining <= cutoff)
final = remaining[locs]
ha = ha[locs]
indexes = indexes[locs]

## Compare the sequence pair indexes to determine each frame.
## Most cases should be fine, but there could be cases that the starting index
## of next frame is higher than the ending index of previous frame.
myfile = open('hbond_index.txt', 'w')
pre = -1
for item in indexes:
    if item > pre:
        myfile.write('%d ' %item)
    else:
        myfile.write('\n%d ' %item)
    pre = item
myfile.close()
print('Hydrogen bond per frame:', len(final)/float(len(traj)))
np.savetxt('hbond_da_dist.txt', final)
np.savetxt('hbond_ha_dist.txt', ha)
print('Total time (s):', time.time()-start)
