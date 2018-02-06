#!/usr/bin/python
# Filename: plot_smd.py
# Description: This is a python script to generate plots to compare force from SMD to AFM measurements
# Date: 11-19-2015 Created

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('f_x.txt')

geom = input('Please specify the geometry of the tip\n\
1   ball\n\
2   disk\n')

if geom == 1:
    r = input('Please input the radius of the ball (nm):\n')
    data[:, 0] = data[:, 0]-0.682-r-0.07
elif geom == 2:
    h = input('Please input the thickness of the disk (nm):\n')
    data[:, 0] = data[:, 0]-0.682-h
else:
    print 'Geometry not recognized\n'
    exit()

file_in = open('50pct.txt','r')
xp = []
for i, line in enumerate(file_in):
    if i > 2:
        xp.append([float(j) for j in line.split()])
xp =np.array(xp)

#plt.figure(figsize=(16, 12),dpi=1000)
fsize = 36
lwidth = 4.0
msize = 16.0

fig, (ax1, ax3) = plt.subplots(2,sharex=True,figsize=(16, 12),dpi=1000)
ax1.errorbar(data[:,0], data[:,1],yerr = data[:,2], fmt='-ro',linewidth=2, markersize= 8,color='blue',markeredgecolor='none', label='MD ')
plt.xlabel('Distance (nm)', fontsize = fsize)
plt.ylabel('Force (kJ/mol/nm)', fontsize = fsize)
ax1.set_ylabel('Force (kJ/mol/nm)', fontsize = fsize)
ax1.tick_params(axis='both',labelsize=fsize, width =2, length=6)
ax1.legend(bbox_to_anchor=(0.98,1), fontsize=fsize,frameon=False)
ax1.set_ylim([0,2500])

ax2 = ax1.twinx()
ax2.plot(xp[:,0], xp[:,1], linewidth=lwidth, color='black', label='AFM')
ax2.set_ylabel('Frequency', fontsize = fsize)
ax2.tick_params(axis='both',labelsize=fsize, width =2, length=6)
ax2.legend(bbox_to_anchor=(1,0.9), fontsize=fsize,frameon=False)

nd = np.loadtxt('numberDensity.txt')
ax3.plot(nd[:,0]-0.682-0.07,nd[:,1], linewidth=lwidth, color='blue', markeredgecolor='none' ,markersize= msize,label='Emim')
ax3.plot(nd[:,0]-0.682-0.07,nd[:,2], linewidth=lwidth, color='red', markeredgecolor='none' ,markersize= msize,label='$\mathregular{BF_4}$')
ax3.plot(nd[:,0]-0.682-0.07,nd[:,3], linewidth=lwidth, color='green', markeredgecolor='none' ,markersize= msize,label='$\mathregular{Tf_2N}$')
ax3.legend(fontsize=fsize,frameon=False)
ax3.set_xlim([-0.5,4])
xticks = np.arange(0,4.1,1)
ax3.xaxis.set_ticks(xticks)
ax3.tick_params(axis='both',labelsize=fsize, width =2, length=6)
ax3.set_xlabel('Distance (nm)', fontsize = fsize)
ax3.set_ylabel('Number density $\mathregular{\#/nm^3}$', fontsize = fsize)

for axis in fig.axes:
    for direction in ['top', 'bottom', 'left', 'right']:
        axis.spines[direction].set_linewidth(2)

plt.tight_layout()
plt.savefig('f_x.tiff', format='tiff')
