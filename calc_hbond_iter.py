#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Description:      Uses the same idea as calc_hbond.py, but load trajectory
#                   in iterative way. As in the case of large trajectory,
#                   the memory load and total cpu time could be reduced.
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
## This loads the whole trajectory just once into the memory.
## May want to modify to use md.iterload to load a chunk of traj each time.
topology = md.load_frame(filename, 0, top = 'begin.gro').top

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

chunk_size = 100
indexes = np.array([np.arange(len(apairs)) for _ in range(chunk_size)])

index_all =  np.array([], dtype = int)
ha_all = np.array([], dtype = float)
da_all = np.array([], dtype = float)

# variable to sum up total frames
tt = 0
for chunk_index, traj in enumerate(md.iterload(filename, chunk = chunk_size, top = 'begin.gro')):
    print('Reading chunk %d of %d frames' %(chunk_index, chunk_size))
    tt += len(traj)
    dists = md.compute_distances(traj, apairs[:, [0, 2]])
    ha_dists = md.compute_distances(traj, apairs[:, [1, 2]])
    angles = md.compute_angles(traj, apairs)
    # Using this way we could also obtain the info for the distance distribution
    #remaining = angles[np.where(dists <= cutoff)]
    #final = remaining[np.where(remaining > acri)]
    locs = np.where(angles > acri)
    da = dists[locs]
    ha = ha_dists[locs]
    ind_temp = indexes[locs]
    locs = np.where(da <= cutoff)
    da_all = np.append(da_all, da[locs])
    ha_all = np.append(ha_all, ha[locs])
    index_all = np.append(index_all, ind_temp[locs])

## Compare the sequence pair indexes to determine each frame.
## Most cases should be fine, but there could be cases that the starting index
## of next frame is higher than the ending index of previous frame.
myfile = open('hbond_index.txt', 'w')
pre = -1
for item in index_all:
    if item > pre:
        myfile.write('%d ' %item)
    else:
        myfile.write('\n%d ' %item)
    pre = item
myfile.close()
print('Hydrogen bond per frame:', len(da_all)/float(tt))
np.savetxt('hbond_da_dist.txt', da_all)
np.savetxt('hbond_ha_dist.txt', ha_all)
print('Total time (s):', time.time()-start)
