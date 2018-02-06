#!/Users/yuzhang/anaconda/envs/myenv/bin/python
import mdtraj as md
import numpy as np
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'filename', help = 'input filename')
args = parser.parse_args()

if args.filename:
    filename = args.filename
else:
    filename = input('Please input filename:\n')

start = time.time()
temp = np.array([0])
traj = md.load(filename, top = 'topol.pdb')
print(len(traj))
topology = traj.top

cutoff = 0.4
acri = 150*np.pi/180
donor_name = ['C3', 'C5', 'C6']
H_name = ['H3', 'H5', 'H6']
#donor_name = ['OG']
#H_name = ['HG']

donor = []
for i in range(1):
    donor.append([topology.select("name %s" %donor_name[i]), topology.select("name %s" %H_name[i])])
#acceptor = topology.select("name =~ 'F.*'")
#acceptor = topology.select("name N1S")
acceptor = topology.select("name OG")

apairs = []
for k in range(1):
    for i in range(len(donor[k][0])):
        for j in acceptor:
            apairs.append([donor[k][0][i], donor[k][1][i], j])
apairs = np.array(apairs)

#
#temp = []
#for i in range(3):
#    temp.append(topology.select_pairs(topology.select("name %s" %donor_name[i]), topology.select("name =~ 'F.*'")))
#pdb.set_trace()
dists = md.compute_distances(traj, apairs[:, [0, 2]])
dist_all = []
hdist_all = []
print(time.time()-start)
for i, frame in enumerate(traj):
    if i%10 == 0:
        print(i)
    remaining = []
    for j, item in enumerate(dists[i]):
        if item <= cutoff:
            remaining.append(j)
    temp = apairs[remaining, :]
    angles = md.compute_angles(frame, temp)[0]
    final = []
    ind = []
    for j, item in enumerate(angles):
        if item >= acri:
            ind.append(remaining[j])
    dist_all.extend(dists[i][ind])
    hdist_all.extend(md.compute_distances(frame, apairs[ind, :][:, [1, 2]])[0])
print('Average distance:', np.mean(dist_all))
print('Hydrogen bonds per frame:', len(dist_all)/float(len(traj)))
print(time.time()-start)
np.savetxt('hbond_da_dist_ho.txt', dist_all)
np.savetxt('hbond_ha_dist_ho.txt', hdist_all)
