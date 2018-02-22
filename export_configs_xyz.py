#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Use python to automatically generate config files for pmf

import argparse, os
import numpy as np

parser = argparse.ArgumentParser(description = 'specify the beginning frame and ending frame')
parser.add_argument('-b', '--begin', type = float, help = 'beginning frame (ps)')
parser.add_argument('-e', '--end', type = float, help = 'ending frame (ps)')
parser.add_argument('-dt', type = float, help = 'time step between each frame')
args = parser.parse_args()

for i, time in enumerate(np.arange(args.begin, args.end+args.dt, args.dt)):
    script = "[atomselect top \"x>210 and x<290\"] writexyz %03d.xyz\nquit\n" %(i+1)
    myfile = open('writexyz.tcl', 'w')
    myfile.write(script)
    myfile.close()
    command = "echo 0 |gmx_mpi trjconv -f ../traj_comp.xtc -s ../topol.tpr -b %s -e %s -o %d.gro" %(time, time, i) 
    os.system(command)
    os.system("vmd %d.gro -dispdev text -e writexyz.tcl" %i)
