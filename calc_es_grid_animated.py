#!/Users/yuzhang/anaconda/bin/python
#Filename: calc_es_grid.py
#Description: This is a python script to generate the potential grids for slit pore models
#Date: 10-09-2017 Created


import numpy as np
import matplotlib.pyplot as plt
import argparse, sys, pdb
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation

parser = argparse.ArgumentParser(description = 'User specified filenames, boundaries, etc.')
parser.add_argument('-b', '--begin', dest = 'begin', default = 0, type = int, help = 'input the first frame')
parser.add_argument('-e', '--end', dest = 'end', default = -1, type = int, help = 'input the last frame')
parser.add_argument('-d', '--dimension', dest = 'dim', nargs = '*', type = float, help = 'input the dimension of the simulaition box')
parser.add_argument('-t', '--timestep', dest = 'timestep', default = 0.1, type = float, help = 'The timestep for each frame')
parser.add_argument('--suffix', type = str, help = 'The suffix for output')
parser.add_argument('-bs', '--bin_size', dest = 'bin_size', type = float, default = 0.1, help = 'the step size for the grid')
args = parser.parse_args()

size = args.bin_size

if args.dim == None:
    print 'No boxsize specified, gonna use the default value: (38, 3.19, 6.1)'
    dim = [38.0, 3.19, 6.1]
else:
    dim = args.dim
if args.suffix == None:
    suffix = ''
else:
    suffix = '_'+args.suffix

dim = [int(np.ceil(i*np.round(1/size))) for i in dim]

fsize = 28
lwidth = 4.0
msize = 16.0

myfile = open('potential_surf_kr.dat', 'r')
chunk = 0
mark = 0
grid_temp = []
grid_2d = []
for line in myfile:
    if 'gfeng' in line:
        chunk += 1
        if grid_temp != []:
            grid_2d.append([[] for row in range(dim[0])])
            for i in range(dim[0]):
                for k in range(dim[2]):
                    temp = 0
                    for j in range(dim[1]):
                        label = i*dim[1]*dim[2]+j*dim[2]+k
                        temp += grid_temp[label]
                    grid_2d[-1][i].append(temp/dim[1])
            grid_2d[-1] = np.array(grid_2d[-1])
            grid_temp = []
        if chunk == args.end:
            break
    elif line == '\n':
        continue
    elif chunk >= args.begin:
        temp = [float(i) for i in line.split()]
        grid_temp.append(sum(temp)*0.010364272)


fig = plt.figure(figsize = (40, 6), dpi = 300)
plt.clf()
im = plt.imshow(grid_2d[0].T, cmap = 'jet', extent = [0, 38, 0, 4.3])
plt.yticks(np.arange(0, 4.3, 1))
plt.xlabel('x (nm)', fontsize = fsize)
plt.ylabel('z (nm)', fontsize = fsize)
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes('right', size='5%', pad=0.5)
cbar = plt.colorbar(im, cax = cax)
cbar.set_label('potential (V)', fontsize = fsize, rotation = 90)
for axes in fig.axes:
    axes.tick_params(labelsize = fsize, width = 2, length = 6, pad = 10)
    for direction in ['top', 'bottom', 'left', 'right']:
        axes.spines[direction].set_linewidth(2)
plt.tight_layout()

def animate(i):
    im.set_data(grid_2d[i].T)
    plt.clim([-10, 10])
    return im

ani = animation.FuncAnimation(fig, animate, range(chunk-args.begin), interval = 25)
writer = animation.writers['ffmpeg'](fps = 1)
ani.save('demo.mp4', writer = writer)
