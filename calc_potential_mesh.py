#!/Users/yuzhang/anaconda/envs/py3/bin/python
# Filename:     calc_potential_mesh.py
# Description:  This is a python script that calculates the potential on carbon
#               onion on a two-electrode simulation, based on the mesh method.
# Date:         01-30-2018

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'Specify some parameters')
parser.add_argument('-b', '--begin', type = int, default = 0, help = 'The starting frame')
parser.add_argument('-v', action = 'store_true', help = 'Show the details')
args = parser.parse_args()

grid = pd.read_csv('gridElectrode1.dat', skiprows = 1, names = ['x', 'y', 'z'], sep = '\s+')

tot = len(grid)
half = int(tot/2)

pot_raw = pd.read_csv('potential_surf_kr.dat', names = ['k', 'r'], sep = '\s+', na_values = 'gfeng', )
pot_raw = pot_raw.dropna(how = 'any')

pot_com = pd.DataFrame()
for i in range(args.begin, int(len(pot_raw)/tot)):
    temp = pot_raw[i*tot: (i+1)*tot]
    temp.index = range(tot)
    pot_com = pd.concat([pot_com, temp['k']+temp['r']], axis = 1)

pot_com.columns = range(args.begin, i+1)
pot_com['avg'] = pot_com.mean(axis = 1)
pot_com['std'] = pot_com.iloc[:, :-1].std(axis = 1)

pot_avg = pot_com[['avg', 'std']]*0.010364272
if args.v:
    print('Total frames:', int(len(pot_raw)/tot))
    print('The potentials for anode and cathode are:')
print(pot_avg[:half]['avg'].mean(), pot_avg[half:]['avg'].mean())

bin_size = 0.05
pot_avg['hist_x'] = np.floor(grid.x/bin_size)*bin_size
pot_avg['hist_y'] = np.floor(grid.y/bin_size)*bin_size
pot_avg['hist_z'] = np.floor(grid.z/bin_size)*bin_size

groupz = pot_avg.groupby('hist_z')['avg']
temp = pd.concat((groupz.mean(), groupz.std()), axis = 1)
temp.to_csv('pot_z.txt', float_format = '%12.6f', header = False)

groupx1 = pot_avg[:half].groupby('hist_x')['avg']
temp = pd.concat((groupx1.mean(), groupx1.std()), axis = 1)
temp.to_csv('pot_x_pos.txt', float_format = '%12.6f', header = False)
groupx2 = pot_avg[half:].groupby('hist_x')['avg']
temp = pd.concat((groupx2.mean(), groupx2.std()), axis = 1)
temp.to_csv('pot_x_neg.txt', float_format = '%12.6f', header = False)

groupy1 = pot_avg[:half].groupby('hist_y')['avg']
temp = pd.concat((groupy1.mean(), groupy1.std()), axis = 1)
temp.to_csv('pot_y_pos.txt', float_format = '%12.6f', header = False)
groupy2 = pot_avg[half:].groupby('hist_y')['avg']
temp = pd.concat((groupy2.mean(), groupy2.std()), axis = 1)
temp.to_csv('pot_y_neg.txt', float_format = '%12.6f', header = False)
