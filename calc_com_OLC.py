#!/Users/yuzhang/anaconda/envs/py3/bin/python 
# Filename:     calc_com_OLC.py
# Description:  This is a python script that calculates the center of mass
#               coordinates of onion like carbon in an supercapacitor
#               system
# Date:         11-22-2017  Created

import calc_common as comm
import mdtraj as md

print('Reading input file', comm.args.filename)
traj = md.load(comm.args.filename)
para = comm.load_para()
res_targ, res_name = comm.select_groups(traj)
xyz_com = comm.calc_xyz_com(res_targ, traj, traj.topology, para)
print('The center of the group:', )
for i in xyz_com:
    print(i)

i = 0
avg = xyz_com[0][0]
for i in range(1, len(xyz_com)):
    avg += xyz_com[i][0]
print(avg/(i+1))
